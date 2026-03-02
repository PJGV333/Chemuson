"""Ejecución aislada de tareas RDKit propensas a fallos nativos."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Any


def run_rdkit_stereo_extract(
    smiles_or_molblock: str,
    fmt: str = "smiles",
    timeout_s: float = 8.0,
) -> dict[str, Any]:
    """Extrae estereoquímica en subproceso a partir de SMILES o molblock."""
    request = {
        "mode": "text",
        "format": str(fmt),
        "value": str(smiles_or_molblock or ""),
    }
    return _run_worker(request, timeout_s=timeout_s)


def stereo_descriptors_for_chain(
    graph,
    chain_atom_ids: list[int],
    timeout_s: float = 8.0,
) -> list[str]:
    """Obtiene designadores estereo para una cadena orientada usando worker aislado."""
    request = _graph_request_payload(
        graph=graph,
        chain_atom_ids=chain_atom_ids,
        mode="graph",
    )
    response = _run_worker(request, timeout_s=timeout_s)
    if not response.get("ok"):
        return []
    descriptors = response.get("descriptors", [])
    if not isinstance(descriptors, list):
        return []
    return [str(item) for item in descriptors if str(item).strip()]


def advanced_stereo_descriptors_for_chain(
    graph,
    chain_atom_ids: list[int],
    timeout_s: float = 8.0,
) -> list[str]:
    """Obtiene designadores estereo avanzados (M/P, R_a/S_a, endo/exo, si/re)."""
    request = _graph_request_payload(
        graph=graph,
        chain_atom_ids=chain_atom_ids,
        mode="advanced_graph",
    )
    response = _run_worker(request, timeout_s=timeout_s)
    if not response.get("ok"):
        return []
    descriptors = response.get("descriptors", [])
    if not isinstance(descriptors, list):
        return []
    return [str(item) for item in descriptors if str(item).strip()]


def _graph_request_payload(
    graph,
    chain_atom_ids: list[int],
    mode: str,
) -> dict[str, Any]:
    """Serializa un grafo interno para el worker aislado de RDKit."""
    return {
        "mode": str(mode),
        "chain": [int(atom_id) for atom_id in chain_atom_ids],
        "atoms": [
            {
                "id": int(atom.id),
                "element": str(atom.element),
                "formal_charge": int(getattr(atom, "formal_charge", getattr(atom, "charge", 0)) or 0),
                "isotope": getattr(atom, "isotope", None),
                "radical_electrons": int(getattr(atom, "radical_electrons", 0) or 0),
                "stereo_cip": getattr(atom, "stereo_cip", None),
                "stereo_axial": getattr(atom, "stereo_axial", None),
                "stereo_helical": getattr(atom, "stereo_helical", None),
                "stereo_si_re": getattr(atom, "stereo_si_re", None),
            }
            for atom in sorted(getattr(graph, "atoms", {}).values(), key=lambda a: a.id)
        ],
        "bonds": [
            {
                "a1_id": int(bond.a1_id),
                "a2_id": int(bond.a2_id),
                "order": int(getattr(bond, "order", 1) or 1),
                "is_aromatic": bool(getattr(bond, "is_aromatic", False)),
                "style": str(getattr(getattr(bond, "style", None), "value", getattr(bond, "style", ""))),
                "donor_atom_id": getattr(bond, "donor_atom_id", None),
                "stereo_ez": getattr(bond, "stereo_ez", None),
                "stereo_axial": getattr(bond, "stereo_axial", None),
                "stereo_endo_exo": getattr(bond, "stereo_endo_exo", None),
                "stereo_helical": getattr(bond, "stereo_helical", None),
            }
            for bond in sorted(getattr(graph, "bonds", {}).values(), key=lambda b: b.id)
        ],
    }


def _worker_path() -> Path:
    return Path(__file__).with_name("_rdkit_worker.py")


def _run_worker(request: dict[str, Any], timeout_s: float) -> dict[str, Any]:
    """Ejecuta worker RDKit y devuelve JSON robusto, incluso ante crash."""
    worker = _worker_path()
    cmd = [sys.executable, str(worker)]
    try:
        proc = subprocess.run(
            cmd,
            input=json.dumps(request),
            text=True,
            capture_output=True,
            timeout=max(1.0, float(timeout_s)),
            check=False,
        )
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": "timeout"}
    except Exception as exc:
        return {"ok": False, "error": str(exc)}

    if proc.returncode != 0:
        if proc.returncode < 0:
            return {
                "ok": False,
                "error": f"worker_exit_signal:{-proc.returncode}",
                "stderr": (proc.stderr or "").strip(),
            }
        return {
            "ok": False,
            "error": f"worker_exit_code:{proc.returncode}",
            "stderr": (proc.stderr or "").strip(),
        }
    try:
        data = json.loads(proc.stdout or "{}")
    except Exception:
        return {
            "ok": False,
            "error": "invalid_worker_json",
            "stdout": (proc.stdout or "").strip(),
            "stderr": (proc.stderr or "").strip(),
        }
    if not isinstance(data, dict):
        return {"ok": False, "error": "invalid_worker_payload"}
    return data
