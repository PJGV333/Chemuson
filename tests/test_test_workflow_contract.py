"""Static contract for the complete, unmasked and bounded pytest campaign."""

from __future__ import annotations

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
WORKFLOW_PATH = ROOT / ".github/workflows/test.yml"
SHARD_TOOL = ROOT / "tools/ci_pytest_shards.py"
EXPECTED_NODEIDS = ROOT / "tests/ci/expected_pytest_nodeids.txt"


def _workflow():
    workflow_text = WORKFLOW_PATH.read_text(encoding="utf-8")
    return workflow_text, yaml.load(workflow_text, Loader=yaml.BaseLoader)


def test_python_ci_installs_runtime_and_development_requirements_before_collection_and_shards() -> None:
    _workflow_text, workflow = _workflow()
    for job_name in ("pytest-plan", "pytest-shard"):
        job = workflow["jobs"][job_name]
        install_step = next(step for step in job["steps"] if step["name"] == "Install dependencies")
        install_commands = install_step["run"]
        assert "python -m pip install -r requirements.txt -r requirements-dev.txt" in install_commands
        assert "python -m pip install -e ." in install_commands
        setup_step = next(step for step in job["steps"] if step["name"] == "Setup Python")
        assert setup_step["with"]["python-version"] == "3.11"
        assert setup_step["with"]["cache"] == "pip"

    dev_requirements = (ROOT / "requirements-dev.txt").read_text(encoding="utf-8")
    runtime_requirements = (ROOT / "requirements.txt").read_text(encoding="utf-8")
    assert any(line.strip().lower() == "pyyaml" for line in dev_requirements.splitlines())
    assert not any(line.strip().lower() == "pyyaml" for line in runtime_requirements.splitlines())


def test_pytest_matrix_executes_every_manifest_id_with_bounded_fail_closed_shards() -> None:
    workflow_text, workflow = _workflow()
    plan_job = workflow["jobs"]["pytest-plan"]
    shard_job = workflow["jobs"]["pytest-shard"]
    summary_job = workflow["jobs"]["pytest-summary"]

    plan_step = next(
        step for step in plan_job["steps"]
        if step["name"] == "Collect and plan the complete pytest suite"
    )
    assert "timeout" in plan_step["run"] and "5m" in plan_step["run"]
    assert "tools/ci_pytest_shards.py plan" in plan_step["run"]
    assert "--shards 8" in plan_step["run"]
    assert "tests/ci/expected_pytest_nodeids.txt" in plan_step["run"]
    assert plan_job["outputs"]["shard-matrix"].endswith(
        "steps.create-plan.outputs.shard-matrix }}"
    )

    matrix = shard_job["strategy"]["matrix"]
    assert "fromJSON(needs.pytest-plan.outputs.shard-matrix)" in matrix
    assert shard_job["strategy"]["fail-fast"] == "false"
    assert shard_job["strategy"]["max-parallel"] == "8"
    run_step = next(
        step for step in shard_job["steps"]
        if step["name"] == "Execute one externally bounded shard"
    )
    assert "tools/ci_pytest_shards.py run" in run_step["run"]
    assert "--timeout-seconds 300" in run_step["run"]
    assert "matrix.shard" in run_step["run"]
    assert run_step["timeout-minutes"] == "6"
    assert any(
        step.get("if") == "always()"
        and "upload-artifact@v4" in step.get("uses", "")
        for step in shard_job["steps"]
    )
    assert summary_job["if"] == "always()"
    assert "pytest-plan" in summary_job["needs"]
    assert "pytest-shard" in summary_job["needs"]
    assert any("Verify full collection" in step["name"] for step in summary_job["steps"])

    tool_text = SHARD_TOOL.read_text(encoding="utf-8")
    assert "pytest_collection_finish" in tool_text
    assert "Counter(actual_nodeids) != Counter(assigned)" in tool_text
    assert "faulthandler_timeout=120" in tool_text
    assert "os.killpg" in tool_text
    assert "pytest collection does not match" in tool_text
    assert "missing shard result reports" in tool_text
    assert "external timeout exceeds 300 seconds" in tool_text
    assert "--ignore" not in workflow_text + tool_text
    assert "--deselect" not in workflow_text + tool_text
    assert "continue-on-error" not in workflow_text
    assert "|| true" not in workflow_text
    assert "python -m pytest -q" not in workflow_text

    assert {"windows-smoke", "flatpak-smoke"}.issubset(workflow["jobs"])
    assert workflow["jobs"]["windows-smoke"]["timeout-minutes"] == "20"
    assert workflow["jobs"]["flatpak-smoke"]["timeout-minutes"] == "10"

    release_gate = (ROOT / ".github/workflows/release.yml").read_text(encoding="utf-8")
    assert release_gate.count("tests/test_test_workflow_contract.py") == 2


def test_junit_nodeid_mapping_keeps_ipv6_parameter_colons_inside_the_test_id() -> None:
    module_spec = spec_from_file_location("ci_pytest_shards_probe", SHARD_TOOL)
    assert module_spec is not None and module_spec.loader is not None
    module = module_from_spec(module_spec)
    module_spec.loader.exec_module(module)
    split_nodeid = module._split_pytest_nodeid

    assert split_nodeid(
        "tests/test_network.py::test_ipv6[http://[::1]:8080/v1]"
    ) == ["tests/test_network.py", "test_ipv6[http://[::1]:8080/v1]"]
    assert split_nodeid(
        "tests/test_network.py::TestNetwork::test_ipv6[http://[2001:db8::1]/v1]"
    ) == [
        "tests/test_network.py",
        "TestNetwork",
        "test_ipv6[http://[2001:db8::1]/v1]",
    ]


def test_pytest_nodeid_manifest_is_sorted_unique_and_covers_the_existing_suite() -> None:
    nodeids = EXPECTED_NODEIDS.read_text(encoding="utf-8").splitlines()

    assert len(nodeids) > 2000
    assert nodeids == sorted(nodeids)
    assert len(nodeids) == len(set(nodeids))
    assert all(nodeid.startswith("tests/") and "::" in nodeid for nodeid in nodeids)
