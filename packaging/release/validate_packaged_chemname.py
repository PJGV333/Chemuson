"""Fail-closed ChemName validation against the actual packaged executable."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from chemuson.chemname.packaged_smoke import (  # noqa: E402
    EXPECTED_MOLECULE_NAMES,
    EXPECTED_TEMPLATE_NAMES,
    run_smoke,
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _source_comparison(frozen_report: dict[str, Any]) -> dict[str, Any]:
    source_report = run_smoke(require_frozen=False)
    _require(source_report.get("ok") is True, f"Source ChemName smoke failed: {source_report.get('failures')}")

    def names(report: dict[str, Any], key: str, id_key: str) -> dict[str, str]:
        entries = report.get(key)
        _require(isinstance(entries, list), f"ChemName smoke lacks {key} results.")
        return {str(item.get(id_key, "")): str(item.get("name", "")) for item in entries if isinstance(item, dict)}

    source_templates = names(source_report, "template_resources", "relative_path")
    frozen_templates = names(frozen_report, "template_resources", "relative_path")
    source_molecules = names(source_report, "molecule_results", "id")
    frozen_molecules = names(frozen_report, "molecule_results", "id")
    _require(source_templates == EXPECTED_TEMPLATE_NAMES, "Source names differ from the approved template expectations.")
    _require(frozen_templates == source_templates, "Frozen template names differ from Python source results.")
    _require(source_molecules == EXPECTED_MOLECULE_NAMES, "Source names differ from the approved molecule expectations.")
    _require(frozen_molecules == source_molecules, "Frozen molecule names differ from Python source results.")
    return {
        "ok": True,
        "template_names_match": True,
        "molecule_names_match": True,
        "source_template_names": source_templates,
        "frozen_template_names": frozen_templates,
        "source_molecule_names": source_molecules,
        "frozen_molecule_names": frozen_molecules,
        "source_smoke": source_report,
    }


def _validate_report(report: dict[str, Any], executable: Path) -> None:
    _require(report.get("ok") is True, f"Frozen ChemName smoke reported failure: {report.get('failures') or report.get('error')}")
    _require(report.get("frozen") is True, "ChemName smoke process did not report a frozen executable.")
    _require(not report.get("gui_modules"), "ChemName smoke imported the ChemUSON GUI.")
    _require(
        Path(str(report.get("executable", ""))).resolve() == executable.resolve(),
        "ChemName smoke ran from a different executable.",
    )
    meipass_text = str(report.get("meipass", ""))
    _require(bool(meipass_text), "Frozen ChemName smoke did not report sys._MEIPASS.")
    meipass = Path(meipass_text).resolve()
    resources = report.get("template_resources")
    _require(isinstance(resources, list), "Frozen ChemName smoke did not report template resources.")
    resource_map = {
        str(item.get("relative_path", "")): item
        for item in resources
        if isinstance(item, dict)
    }
    _require(set(resource_map) == set(EXPECTED_TEMPLATE_NAMES), "Frozen ChemName resource inventory is incomplete or unexpected.")
    for relative_path, expected_name in EXPECTED_TEMPLATE_NAMES.items():
        item = resource_map[relative_path]
        _require(item.get("status") == "pass", f"Template failed in frozen executable: {relative_path}.")
        _require(item.get("exists") is True and int(item.get("size_bytes", 0)) > 0, f"Frozen template is missing or empty: {relative_path}.")
        _require(item.get("name") == expected_name, f"Frozen template name mismatch for {relative_path}.")
        _require(
            Path(str(item.get("path", ""))).resolve().is_relative_to(meipass),
            f"Frozen template {relative_path} resolved outside sys._MEIPASS.",
        )
        digest = str(item.get("sha256", ""))
        _require(len(digest) == 64 and all(char in "0123456789abcdef" for char in digest.lower()), f"Frozen template hash is invalid: {relative_path}.")

    molecules = report.get("molecule_results")
    _require(isinstance(molecules, list), "Frozen ChemName smoke did not report molecule results.")
    molecule_map = {
        str(item.get("id", "")): item
        for item in molecules
        if isinstance(item, dict)
    }
    _require(set(molecule_map) == set(EXPECTED_MOLECULE_NAMES), "Frozen molecule smoke inventory is incomplete or unexpected.")
    for case_id, expected_name in EXPECTED_MOLECULE_NAMES.items():
        item = molecule_map[case_id]
        _require(item.get("status") == "pass", f"Frozen molecule naming failed for {case_id}.")
        _require(item.get("name") == expected_name, f"Frozen molecule name mismatch for {case_id}.")


def validate(executable: Path, *, timeout: int = 120) -> dict[str, Any]:
    executable = executable.resolve(strict=True)
    if not executable.is_file():
        raise ValueError(f"Frozen executable is not a file: {executable}")
    source_report = run_smoke(require_frozen=False)
    _require(source_report.get("ok") is True, f"Python source ChemName smoke failed: {source_report.get('failures')}")

    with tempfile.TemporaryDirectory(prefix="chemuson-frozen-chemname-smoke-") as temporary:
        scratch = Path(temporary)
        report_path = scratch / "chemname-smoke-report.json"
        environment = os.environ.copy()
        environment["CHEMUSON_CHEMNAME_PACKAGED_SMOKE"] = "1"
        for variable in (
            "CHEMUSON_INTERNAL_RDKIT_WORKER",
            "CHEMUSON_RDKIT_PACKAGED_SMOKE",
            "CHEMUSON_ICON_SMOKE_TEST",
        ):
            environment.pop(variable, None)
        command = [
            str(executable),
            "--chemname-packaged-smoke-test",
            "--chemname-smoke-report",
            str(report_path),
        ]
        try:
            process = subprocess.run(
                command,
                cwd=scratch,
                env=environment,
                stdin=subprocess.DEVNULL,
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False,
            )
        except subprocess.TimeoutExpired as exc:
            raise ValueError(f"Frozen ChemName smoke timed out after {timeout}s.") from exc
        except OSError as exc:
            raise ValueError(f"Could not start frozen executable {executable}: {exc}") from exc
        if not report_path.is_file():
            raise ValueError(
                "Frozen executable did not write its ChemName smoke report "
                f"(exit={process.returncode}, stdout={process.stdout[-2000:]!r}, stderr={process.stderr[-2000:]!r})."
            )
        try:
            report = json.loads(report_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError("Frozen executable wrote an invalid ChemName smoke report.") from exc
        if not isinstance(report, dict):
            raise ValueError("Frozen ChemName smoke report must be a JSON object.")
        if process.returncode != 0:
            raise ValueError(
                f"Frozen executable ChemName smoke exited {process.returncode}: "
                f"{report.get('failures') or report.get('error') or process.stderr[-2000:]}"
            )
        _validate_report(report, executable)
        report["source_comparison"] = _source_comparison(report)
        report["executable_sha256"] = hashlib.sha256(executable.read_bytes()).hexdigest()
        report["ok"] = True
        return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    report = validate(args.executable, timeout=args.timeout)
    serialized = json.dumps(report, sort_keys=True)
    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(serialized + "\n", encoding="utf-8")
    print(serialized)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"Packaged ChemName validation failed: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
