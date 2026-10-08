"""Fail-closed RDKit smoke executed by the frozen packaging gate."""

from __future__ import annotations

import math
import sys
import time
from typing import Any

from chemuson.chemio.rdkit_safe import (
    conformer_3d_isolated,
    molecular_descriptors_isolated,
    molgraph_to_smiles_isolated,
    rdkit_worker_diagnostics,
    text_to_molblock,
)
from chemuson.core.model import MolGraph


def _loaded_rdkit_modules() -> list[str]:
    return sorted(
        name for name in sys.modules if name == "rdkit" or name.startswith("rdkit.")
    )


def _loaded_gui_modules() -> list[str]:
    return sorted(
        name
        for name in sys.modules
        if name == "PyQt6"
        or name.startswith("PyQt6.")
        or name == "chemuson.gui"
        or name.startswith("chemuson.gui.")
    )


def _ethanol_graph() -> MolGraph:
    graph = MolGraph()
    carbon_1 = graph.add_atom("C", 0.0, 0.0)
    carbon_2 = graph.add_atom("C", 40.0, 0.0)
    oxygen = graph.add_atom("O", 80.0, 0.0)
    graph.add_bond(carbon_1.id, carbon_2.id, order=1)
    graph.add_bond(carbon_2.id, oxygen.id, order=1)
    return graph


def run_smoke() -> dict[str, Any]:
    """Exercise bundled native imports and isolated RDKit operations only."""
    started = time.monotonic()
    report: dict[str, Any] = {
        "ok": False,
        "executable": sys.executable,
        "frozen": bool(getattr(sys, "frozen", False)),
        "meipass": str(getattr(sys, "_MEIPASS", "")),
        "expected_executable": sys.executable,
        "parent_rdkit_before": _loaded_rdkit_modules(),
        "parent_gui_before": _loaded_gui_modules(),
    }
    failures: list[str] = []
    if not report["frozen"]:
        failures.append("not_frozen")
    if report["parent_rdkit_before"]:
        failures.append("parent_loaded_rdkit_before_smoke")
    if report["parent_gui_before"]:
        failures.append("parent_loaded_gui_before_smoke")

    diagnostics = rdkit_worker_diagnostics(timeout_s=15.0)
    report["worker_diagnostics"] = diagnostics
    if diagnostics.get("ok") is not True:
        failures.append(f"worker_diagnostics:{diagnostics.get('error', 'failed')}")
    if diagnostics.get("python_executable") != sys.executable:
        failures.append("worker_executable_mismatch")
    if diagnostics.get("worker_frozen") is not True:
        failures.append("worker_not_frozen")
    worker_meipass = str(diagnostics.get("worker_meipass", ""))
    if not worker_meipass:
        failures.append("worker_meipass_missing")
    extensions = diagnostics.get("native_extensions", {})
    required_extensions = (
        "rdkit.Chem.rdchem",
        "rdkit.Chem.rdMolDescriptors",
        "rdkit.Chem.rdDistGeom",
    )
    for name in required_extensions:
        extension = extensions.get(name, {}) if isinstance(extensions, dict) else {}
        if extension.get("exists") is not True:
            failures.append(f"native_extension_missing:{name}")
        if extension.get("native") is not True:
            failures.append(f"native_extension_not_compiled:{name}")
        if extension.get("inside_bundle") is not True:
            failures.append(f"native_extension_outside_bundle:{name}")

    graph = _ethanol_graph()
    descriptors, descriptor_error = molecular_descriptors_isolated(graph, timeout_s=15.0)
    report["descriptors"] = descriptors
    if descriptor_error:
        failures.append(f"descriptors:{descriptor_error}")
    elif descriptors is None:
        failures.append("descriptors:empty")
    else:
        expected = {
            "logp": (-0.0014, 1e-4),
            "tpsa": (20.23, 1e-6),
            "molecular_weight": (46.069, 1e-6),
        }
        for key, (value, tolerance) in expected.items():
            try:
                actual = float(descriptors[key])
            except (KeyError, TypeError, ValueError):
                failures.append(f"descriptor_missing:{key}")
                continue
            if not math.isfinite(actual) or abs(actual - value) > tolerance:
                failures.append(f"descriptor_mismatch:{key}:{actual}")
        for key in ("hbd", "hba"):
            if descriptors.get(key) != 1:
                failures.append(f"descriptor_mismatch:{key}:{descriptors.get(key)}")

    smiles, smiles_error = molgraph_to_smiles_isolated(graph, timeout_s=15.0)
    report["canonical_smiles"] = smiles
    if smiles_error or smiles != "CCO":
        failures.append(f"smiles:{smiles_error or smiles or 'empty'}")
    smiles_import = text_to_molblock("CCO", fmt="smiles", timeout_s=15.0)
    report["smiles_input"] = {
        "ok": smiles_import.get("ok") is True,
        "molblock_chars": len(str(smiles_import.get("molblock", ""))),
        "error": smiles_import.get("error", ""),
    }
    if report["smiles_input"]["ok"] is not True or report["smiles_input"]["molblock_chars"] == 0:
        failures.append(f"smiles_input:{smiles_import.get('error', 'empty_molblock')}")

    coordinates, metadata, conformer_error = conformer_3d_isolated(graph, timeout_s=20.0)
    report["conformer3d"] = {
        "metadata": metadata,
        "atom_count": len(coordinates or {}),
        "positions": {
            str(atom_id): [float(component) for component in position]
            for atom_id, position in (coordinates or {}).items()
        },
    }
    if conformer_error or coordinates is None:
        failures.append(f"conformer3d:{conformer_error or 'empty'}")
    elif set(coordinates) != {atom.id for atom in graph.atoms.values()}:
        failures.append("conformer3d:atom_ids_mismatch")
    elif any(
        len(position) != 3 or not all(math.isfinite(float(value)) for value in position)
        for position in coordinates.values()
    ):
        failures.append("conformer3d:invalid_coordinates")

    report["parent_rdkit_after"] = _loaded_rdkit_modules()
    report["parent_gui_after"] = _loaded_gui_modules()
    if report["parent_rdkit_after"]:
        failures.append("parent_loaded_rdkit_after_smoke")
    if report["parent_gui_after"]:
        failures.append("parent_loaded_gui_after_smoke")
    report["elapsed_seconds"] = round(time.monotonic() - started, 3)
    report["failures"] = failures
    report["ok"] = not failures
    return report
