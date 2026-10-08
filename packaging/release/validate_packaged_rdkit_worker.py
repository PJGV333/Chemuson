"""Fail-closed smoke for RDKit in the actual PyInstaller executable."""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _validate_report(report: dict[str, Any], executable: Path) -> None:
    _require(report.get("ok") is True, f"Frozen RDKit smoke reported failure: {report.get('failures') or report.get('error')}")
    _require(report.get("frozen") is True, "Smoke process did not report a frozen executable.")
    actual_executable = Path(str(report.get("executable", ""))).resolve()
    _require(actual_executable == executable.resolve(), "Smoke ran from a different executable.")
    _require(not report.get("parent_rdkit_before"), "RDKit was imported in the parent before worker calls.")
    _require(not report.get("parent_rdkit_after"), "RDKit was imported in the parent during worker calls.")
    _require(not report.get("parent_gui_before"), "GUI was imported before packaged smoke dispatch.")
    _require(not report.get("parent_gui_after"), "GUI was imported by the worker smoke parent.")

    diagnostics = report.get("worker_diagnostics")
    _require(isinstance(diagnostics, dict) and diagnostics.get("ok") is True, "RDKit worker import diagnostics failed.")
    worker_executable = Path(str(diagnostics.get("python_executable", ""))).resolve()
    _require(worker_executable == executable.resolve(), "RDKit worker used a different Python executable.")
    _require(diagnostics.get("worker_frozen") is True, "RDKit worker did not run in frozen mode.")
    worker_meipass = str(diagnostics.get("worker_meipass", ""))
    _require(bool(worker_meipass), "RDKit worker did not report its bundle extraction root.")
    bundle_root = Path(worker_meipass).resolve()
    extensions = diagnostics.get("native_extensions")
    _require(isinstance(extensions, dict), "RDKit worker did not report its native extensions.")
    for name in (
        "rdkit.Chem.rdchem",
        "rdkit.Chem.rdMolDescriptors",
        "rdkit.Chem.rdDistGeom",
    ):
        extension = extensions.get(name)
        _require(isinstance(extension, dict), f"RDKit extension {name} was not reported.")
        _require(extension.get("exists") is True, f"RDKit extension {name} did not load from a file.")
        _require(extension.get("native") is True, f"RDKit extension {name} is not compiled native code.")
        _require(extension.get("inside_bundle") is True, f"RDKit extension {name} loaded outside the executable bundle.")
        extension_path = Path(str(extension.get("path", ""))).resolve()
        _require(
            extension_path.is_relative_to(bundle_root),
            f"RDKit extension {name} path is outside the worker bundle root.",
        )

    descriptors = report.get("descriptors")
    _require(isinstance(descriptors, dict), "Ethanol descriptor calculation returned no values.")
    for key, expected, tolerance in (
        ("logp", -0.0014, 1e-4),
        ("tpsa", 20.23, 1e-6),
        ("molecular_weight", 46.069, 1e-6),
    ):
        value = float(descriptors.get(key, float("nan")))
        _require(math.isfinite(value) and abs(value - expected) <= tolerance, f"Ethanol descriptor {key} was {value!r}, expected {expected}.")
    _require(descriptors.get("hbd") == 1, "Ethanol HBD must equal 1.")
    _require(descriptors.get("hba") == 1, "Ethanol HBA must equal 1.")
    _require(report.get("canonical_smiles") == "CCO", "Isolated SMILES worker did not return canonical ethanol CCO.")
    smiles_input = report.get("smiles_input")
    _require(
        isinstance(smiles_input, dict)
        and smiles_input.get("ok") is True
        and int(smiles_input.get("molblock_chars", 0)) > 0,
        "Isolated SMILES parser worker did not return an ethanol molecule.",
    )

    conformer = report.get("conformer3d")
    _require(isinstance(conformer, dict) and conformer.get("atom_count") == 3, "Isolated 3D worker did not return coordinates for all ethanol atoms.")
    positions = conformer.get("positions")
    _require(isinstance(positions, dict) and len(positions) == 3, "Isolated 3D worker returned malformed coordinates.")
    for atom_id, position in positions.items():
        _require(
            isinstance(position, list)
            and len(position) == 3
            and all(isinstance(value, (int, float)) and math.isfinite(value) for value in position),
            f"Isolated 3D coordinates for atom {atom_id} are invalid.",
        )


def validate(executable: Path, *, timeout: int = 120) -> dict[str, Any]:
    executable = executable.resolve(strict=True)
    if not executable.is_file():
        raise ValueError(f"Frozen executable is not a file: {executable}")
    with tempfile.TemporaryDirectory(prefix="chemuson-frozen-rdkit-smoke-") as temporary:
        scratch = Path(temporary)
        report_path = scratch / "rdkit-smoke-report.json"
        environment = os.environ.copy()
        environment["CHEMUSON_RDKIT_PACKAGED_SMOKE"] = "1"
        environment.pop("CHEMUSON_INTERNAL_RDKIT_WORKER", None)
        environment.pop("CHEMUSON_ENABLE_DIRECT_RDKIT", None)
        command = [
            str(executable),
            "--rdkit-packaged-smoke-test",
            "--rdkit-smoke-report",
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
            raise ValueError(f"Frozen RDKit smoke timed out after {timeout}s.") from exc
        except OSError as exc:
            raise ValueError(f"Could not start frozen executable {executable}: {exc}") from exc
        if not report_path.is_file():
            raise ValueError(
                "Frozen executable did not write its RDKit smoke report "
                f"(exit={process.returncode}, stdout={process.stdout[-2000:]!r}, stderr={process.stderr[-2000:]!r})."
            )
        try:
            report = json.loads(report_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError("Frozen executable wrote an invalid RDKit smoke report.") from exc
        if not isinstance(report, dict):
            raise ValueError("Frozen executable RDKit smoke report must be a JSON object.")
        if process.returncode != 0:
            raise ValueError(
                f"Frozen executable RDKit smoke exited {process.returncode}: "
                f"{report.get('failures') or report.get('error') or process.stderr[-2000:]}"
            )
        _validate_report(report, executable)
        return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=120)
    args = parser.parse_args()
    report = validate(args.executable, timeout=args.timeout)
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"Packaged RDKit worker validation failed: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
