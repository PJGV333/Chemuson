"""Fail-closed frozen executable tests for ChemName assets and output."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "packaging" / "release"))

from chemuson.chemname.packaged_smoke import run_smoke  # noqa: E402
import validate_packaged_chemname as validator  # noqa: E402


def _valid_report(executable: Path) -> dict[str, Any]:
    report = run_smoke(require_frozen=False)
    bundle = executable.parent / "_MEI123"
    report.update(
        {
            "frozen": True,
            "executable": str(executable.resolve()),
            "meipass": str(bundle.resolve()),
            "gui_modules": [],
        }
    )
    for resource in report["template_resources"]:
        resource["path"] = str(bundle / "chemuson" / "chemname" / resource["relative_path"])
    return report


def test_frozen_chemname_validator_accepts_only_complete_in_bundle_report(
    tmp_path: Path,
) -> None:
    executable = tmp_path / "Chemuson.exe"
    executable.touch()
    report = _valid_report(executable)

    validator._validate_report(report, executable)

    report["template_resources"].pop()
    with pytest.raises(ValueError, match="inventory is incomplete"):
        validator._validate_report(report, executable)


def test_validator_executes_the_exact_binary_and_compares_against_source(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    executable = tmp_path / "Chemuson"
    executable.touch()
    observed: dict[str, Any] = {"commands": []}
    original_run = subprocess.run

    def fake_run(command: list[str], **kwargs: Any) -> subprocess.CompletedProcess[str]:
        observed["commands"].append(command)
        if command[0] != str(executable.resolve()):
            return original_run(command, **kwargs)
        observed["environment"] = kwargs["env"]
        report_path = Path(command[-1])
        report_path.write_text(json.dumps(_valid_report(executable)), encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(validator.subprocess, "run", fake_run)
    report = validator.validate(executable)

    frozen_command = next(
        command for command in observed["commands"] if command[0] == str(executable.resolve())
    )
    assert frozen_command[1] == "--chemname-packaged-smoke-test"
    assert observed["environment"]["CHEMUSON_CHEMNAME_PACKAGED_SMOKE"] == "1"
    assert report["source_comparison"]["ok"] is True
    assert report["source_comparison"]["template_names_match"] is True
    assert report["source_comparison"]["molecule_names_match"] is True
