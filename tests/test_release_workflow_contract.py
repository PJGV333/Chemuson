"""Static contract tests for the official tag-only release workflow."""

from __future__ import annotations

from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
WORKFLOW_PATH = ROOT / ".github/workflows/release.yml"


def _workflow() -> tuple[str, dict]:
    text = WORKFLOW_PATH.read_text(encoding="utf-8")
    workflow = yaml.load(text, Loader=yaml.BaseLoader)
    assert isinstance(workflow, dict)
    return text, workflow


def test_official_release_is_tag_only_and_defaults_to_read_permissions() -> None:
    text, workflow = _workflow()
    assert workflow["on"] == {"push": {"tags": ["v*"]}}
    assert workflow["permissions"] == {"contents": "read"}
    gate = workflow["jobs"]["release_gate"]
    assert gate["permissions"] == {"contents": "read"}
    preflight = next(step for step in gate["steps"] if step.get("id") == "preflight")
    assert preflight["env"]["RELEASE_PREFLIGHT_TOKEN"] == "${{ github.token }}"
    assert "workflow_dispatch" not in text
    assert "inputs.version" not in text
    assert "inputs.channel" not in text


def test_all_platform_builds_depend_on_gate_and_checkout_its_exact_sha() -> None:
    _text, workflow = _workflow()
    jobs = workflow["jobs"]
    for name in ("build_windows", "build_linux", "build_flatpak"):
        job = jobs[name]
        assert job["needs"] == "release_gate"
        checkout = next(step for step in job["steps"] if step.get("uses", "").startswith("actions/checkout@"))
        assert checkout["with"]["ref"] == "${{ needs.release_gate.outputs.sha }}"
        assert job["permissions"] == {"contents": "read"}
    assert jobs["release"]["needs"] == ["release_gate", "build_windows", "build_linux", "build_flatpak"]
    assert "release_gate" in jobs["publish_flatpak_remote"]["needs"]


def test_only_publication_jobs_have_write_permission_and_are_explicit() -> None:
    _text, workflow = _workflow()
    jobs = workflow["jobs"]
    writable = {
        name
        for name, job in jobs.items()
        if job.get("permissions", {}).get("contents") == "write"
    }
    assert writable == {"release", "publish_flatpak_remote"}
    assert jobs["release"]["environment"]["name"] == "${{ needs.release_gate.outputs.channel }}"
    assert jobs["publish_flatpak_remote"]["environment"]["name"] == (
        "flatpak-${{ needs.release_gate.outputs.channel }}"
    )


def test_existing_release_build_and_flatpak_publication_checks_remain() -> None:
    text, _workflow_data = _workflow()
    for required in (
        "validate_release_tag.py",
        "validate_release_artifacts.py",
        "generate_checksums.py",
        "sign_artifacts.ps1",
        "validate_flatpak_remote_artifacts.py",
        "Validate local Flatpak commit metadata",
        "Verify GitHub Pages publication",
        "softprops/action-gh-release@v2",
        "overwrite_files: false",
        "fail_on_unmatched_files: true",
    ):
        assert required in text
    assert "set_version.py" not in text
    assert "target_commitish:" not in text


def test_inno_requires_an_explicit_version_and_existing_smoke_supplies_one() -> None:
    setup_script = (ROOT / "packaging/windows/Chemuson.iss").read_text(encoding="utf-8")
    test_workflow = (ROOT / ".github/workflows/test.yml").read_text(encoding="utf-8")
    assert "#error" in setup_script
    assert 'MyAppName "ChemUSON"' in setup_script
    assert 'MyAppPublisher "ChemUSON"' in setup_script
    assert 'MyAppExeName "Chemuson.exe"' in setup_script
    assert "0.0.0-dev" not in setup_script
    assert '$env:CHEMUSON_VERSION = "0.0.0-ci"' in test_workflow


def test_gate_uses_bounded_non_monolithic_release_checks() -> None:
    text, workflow = _workflow()
    gate_text = text.split("  build_windows:", 1)[0]
    assert "timeout 5m python -m pytest -q tests/architecture" in gate_text
    assert "timeout 8m python -m pytest -q" in gate_text
    assert "pytest -q\n          tests/" in gate_text
    assert "pytest -q" in workflow["jobs"]["release_gate"]["steps"][-1]["run"]
