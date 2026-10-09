"""Runtime acceptance smoke for ChemName package resources and names."""

from __future__ import annotations

import hashlib
import sys
import time
from pathlib import Path
from typing import Any

from chemuson.chemio.rdkit_io import molfile_to_molgraph
from chemuson.chemname import NameOptions, iupac_name
from chemuson.core.model import MolGraph
from chemuson.platform.resources import open_resource_path

EXPECTED_TEMPLATE_NAMES = {
    "templates/fused/pyrene_cas.mol": "pyrene",
    "templates/fused/pyrene_iupac2004.mol": "pyrene",
    "templates/simple/benzene.mol": "benzene",
    "templates/special/alpha_d_glucopyranose.mol": "alpha-d-glucopyranose",
    "templates/special/androstane_core.mol": "androstane",
    "templates/special/beta_d_fructofuranose.mol": "beta-d-fructofuranose",
    "templates/special/beta_d_glucopyranose.mol": "beta-d-glucopyranose",
    "templates/special/cholestane_core.mol": "cholestane",
    "templates/special/d_ribose.mol": "d-ribose",
}

EXPECTED_MOLECULE_NAMES = {
    "ethanol": "ethan-1-ol",
    "benzene": "benzene",
    "acetamide": "1-aminoethanamide",
    "ethane": "ethane",
    "cyclohexane": "cyclohexane",
    "unsupported_element": "N/D",
}


def _graph(elements: tuple[str, ...], bonds: tuple[tuple[int, int, int], ...]) -> MolGraph:
    graph = MolGraph()
    atoms = [graph.add_atom(element, float(index), 0.0) for index, element in enumerate(elements)]
    for first, second, order in bonds:
        graph.add_bond(atoms[first].id, atoms[second].id, order=order)
    return graph


def _molecule_cases() -> dict[str, MolGraph]:
    benzene = MolGraph()
    ring = [benzene.add_atom("C", float(index), 0.0) for index in range(6)]
    for index in range(6):
        benzene.add_bond(
            ring[index].id,
            ring[(index + 1) % 6].id,
            order=1,
            is_aromatic=True,
        )
    return {
        "ethanol": _graph(("C", "C", "O"), ((0, 1, 1), (1, 2, 1))),
        "benzene": benzene,
        "acetamide": _graph(("C", "C", "O", "N"), ((0, 1, 1), (1, 2, 2), (1, 3, 1))),
        "ethane": _graph(("C", "C"), ((0, 1, 1),)),
        "cyclohexane": _graph(
            ("C", "C", "C", "C", "C", "C"),
            ((0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 0, 1)),
        ),
        "unsupported_element": _graph(("Xx",), ()),
    }


def _exception_fields(exc: Exception) -> dict[str, str]:
    cause = exc.__cause__
    return {
        "exception_type": type(exc).__name__,
        "exception": str(exc),
        "cause_type": type(cause).__name__ if cause is not None else "",
        "cause": str(cause) if cause is not None else "",
    }


def run_smoke(
    *,
    require_frozen: bool = True,
    expected_package_root: str | Path | None = None,
) -> dict[str, Any]:
    """Check every ChemName MOL template and representative molecular names.

    ``expected_package_root`` is used by Flatpak's build sandbox to prove that
    imports and resources resolve from the installed ``/app`` tree.
    """
    started = time.perf_counter()
    frozen = bool(getattr(sys, "frozen", False))
    meipass_value = getattr(sys, "_MEIPASS", None)
    meipass = Path(str(meipass_value)).resolve() if frozen and meipass_value else None
    package_file = Path(__file__).resolve()
    expected_root = Path(expected_package_root).resolve() if expected_package_root else None
    gui_modules = sorted(
        name
        for name in sys.modules
        if name == "chemuson.gui" or name.startswith("chemuson.gui.")
    )
    report: dict[str, Any] = {
        "ok": False,
        "frozen": frozen,
        "gui_modules": gui_modules,
        "executable": str(Path(sys.executable).resolve()),
        "meipass": str(meipass or ""),
        "package_file": str(package_file),
        "expected_package_root": str(expected_root or ""),
        "template_resources": [],
        "molecule_results": [],
        "failures": [],
    }
    failures: list[str] = report["failures"]

    if require_frozen and not frozen:
        failures.append("not_frozen")
    if require_frozen and meipass is None:
        failures.append("meipass_missing")
    if require_frozen and gui_modules:
        failures.append("chemuson_gui_imported")
    if expected_root is not None and not package_file.is_relative_to(expected_root):
        failures.append("package_module_outside_expected_root")

    for relative_path, expected_name in sorted(EXPECTED_TEMPLATE_NAMES.items()):
        resource: dict[str, Any] = {
            "relative_path": relative_path,
            "expected_name": expected_name,
            "exists": False,
            "size_bytes": 0,
            "sha256": "",
            "path": "",
            "name": "",
            "status": "error",
        }
        try:
            parts = ("chemname", *Path(relative_path).parts)
            with open_resource_path(*parts) as resource_path:
                path = Path(resource_path).resolve()
                content = path.read_bytes()
                graph = molfile_to_molgraph(content.decode("utf-8", errors="replace"))
                name = iupac_name(
                    graph,
                    NameOptions(
                        return_nd_on_fail=False,
                        fused_numbering_scheme=("cas" if relative_path.endswith("pyrene_cas.mol") else "iupac2004"),
                    ),
                )
                resource.update(
                    {
                        "exists": path.is_file(),
                        "size_bytes": len(content),
                        "sha256": hashlib.sha256(content).hexdigest(),
                        "path": str(path),
                        "name": name,
                        "status": "pass" if name == expected_name else "fail",
                    }
                )
                if not path.is_file() or not content:
                    failures.append(f"template_missing_or_empty:{relative_path}")
                if name != expected_name:
                    failures.append(f"template_name_mismatch:{relative_path}:{name!r}")
                if meipass is not None and not path.is_relative_to(meipass):
                    failures.append(f"template_outside_bundle:{relative_path}")
                if expected_root is not None and not path.is_relative_to(expected_root):
                    failures.append(f"template_outside_expected_root:{relative_path}")
        except Exception as exc:
            resource.update(_exception_fields(exc))
            failures.append(f"template_error:{relative_path}:{type(exc).__name__}:{exc}")
        report["template_resources"].append(resource)

    for case_id, graph in _molecule_cases().items():
        expected_name = EXPECTED_MOLECULE_NAMES[case_id]
        result: dict[str, Any] = {
            "id": case_id,
            "expected_name": expected_name,
            "name": "",
            "status": "error",
        }
        try:
            options = NameOptions(return_nd_on_fail=(case_id == "unsupported_element"))
            name = iupac_name(graph, options)
            result.update({"name": name, "status": "pass" if name == expected_name else "fail"})
            if name != expected_name:
                failures.append(f"molecule_name_mismatch:{case_id}:{name!r}")
        except Exception as exc:
            result.update(_exception_fields(exc))
            failures.append(f"molecule_error:{case_id}:{type(exc).__name__}:{exc}")
        report["molecule_results"].append(result)

    report["elapsed_seconds"] = round(time.perf_counter() - started, 3)
    report["ok"] = not failures
    return report
