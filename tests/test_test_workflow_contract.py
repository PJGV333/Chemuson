"""Static contract for CI test dependencies and unmasked test execution."""

from __future__ import annotations

from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
WORKFLOW_PATH = ROOT / ".github/workflows/test.yml"


def test_python_ci_installs_runtime_and_development_requirements_before_pytest() -> None:
    workflow_text = WORKFLOW_PATH.read_text(encoding="utf-8")
    workflow = yaml.load(workflow_text, Loader=yaml.BaseLoader)
    python_job = workflow["jobs"]["pytest"]
    install_step = next(step for step in python_job["steps"] if step["name"] == "Install dependencies")
    install_commands = install_step["run"]
    assert "python -m pip install -r requirements.txt -r requirements-dev.txt" in install_commands
    assert "python -m pip install -e ." in install_commands

    dev_requirements = (ROOT / "requirements-dev.txt").read_text(encoding="utf-8")
    runtime_requirements = (ROOT / "requirements.txt").read_text(encoding="utf-8")
    assert any(line.strip().lower() == "pyyaml" for line in dev_requirements.splitlines())
    assert not any(line.strip().lower() == "pyyaml" for line in runtime_requirements.splitlines())


def test_python_ci_runs_real_pytest_without_failure_masking_or_broad_exclusions() -> None:
    workflow_text = WORKFLOW_PATH.read_text(encoding="utf-8")
    workflow = yaml.load(workflow_text, Loader=yaml.BaseLoader)
    python_job = workflow["jobs"]["pytest"]
    run_step = next(step for step in python_job["steps"] if step["name"] == "Run tests")
    assert run_step["run"].strip() == "python -m pytest -q"
    assert "continue-on-error" not in workflow_text
    assert "--ignore" not in workflow_text
    assert "--deselect" not in workflow_text
    assert "|| true" not in workflow_text

    release_gate = (ROOT / ".github/workflows/release.yml").read_text(encoding="utf-8")
    assert release_gate.count("tests/test_test_workflow_contract.py") == 2
