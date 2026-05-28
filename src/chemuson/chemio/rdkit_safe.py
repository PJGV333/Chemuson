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


def text_to_molblock(
    value: str,
    fmt: str = "smiles",
    timeout_s: float = 8.0,
) -> dict[str, Any]:
    """Convierte SMILES/molblock a molblock mediante worker RDKit aislado."""
    request = {
        "mode": "to_molblock",
        "format": str(fmt),
        "value": str(value or ""),
    }
    return _run_worker(request, timeout_s=timeout_s)


def smiles_to_molgraph_isolated(
    smiles: str,
    timeout_s: float = 8.0,
):
    """Convierte SMILES a `MolGraph` usando RDKit aislado + parser local."""
    response = text_to_molblock(smiles, fmt="smiles", timeout_s=timeout_s)
    if not response.get("ok"):
        return None, str(response.get("error", "worker_error"))
    molblock = str(response.get("molblock", "") or "")
    if not molblock.strip():
        return None, "empty_molblock"
    try:
        # Import local para evitar ciclos en carga de módulo.
        from chemuson.chemio.rdkit_io import molfile_to_molgraph

        return molfile_to_molgraph(molblock), None
    except Exception as exc:
        return None, str(exc)


def molgraph_to_smiles_isolated(
    graph,
    timeout_s: float = 5.0,
) -> tuple[str | None, str | None]:
    """Convierte un `MolGraph` a SMILES usando RDKit en un proceso aislado."""
    request = _graph_request_payload(
        graph=graph,
        chain_atom_ids=[],
        mode="graph_to_smiles",
    )
    response = _run_worker(request, timeout_s=timeout_s)
    if not response.get("ok"):
        return None, str(response.get("error", "worker_error"))
    smiles = str(response.get("smiles", "") or "").strip()
    if not smiles:
        return None, "empty_smiles"
    return smiles, None


def molecular_descriptors_isolated(
    graph,
    timeout_s: float = 5.0,
) -> tuple[dict[str, Any] | None, str | None]:
    """Calcula descriptores RDKit en un worker aislado."""
    request = _graph_request_payload(
        graph=graph,
        chain_atom_ids=[],
        mode="graph_descriptors",
    )
    response = _run_worker(request, timeout_s=timeout_s)
    if not response.get("ok"):
        return None, str(response.get("error", "worker_error"))
    descriptors = response.get("descriptors", {})
    if not isinstance(descriptors, dict):
        return None, "invalid_descriptors"
    return descriptors, None


def conformer_3d_isolated(
    graph,
    timeout_s: float = 8.0,
) -> tuple[dict[int, tuple[float, float, float]] | None, dict[str, Any], str | None]:
    """Genera conformación 3D RDKit en un worker aislado."""
    request = _graph_request_payload(
        graph=graph,
        chain_atom_ids=[],
        mode="graph_conformer3d",
    )
    response = _run_worker(request, timeout_s=timeout_s)
    if not response.get("ok"):
        return None, {}, str(response.get("error", "worker_error"))
    raw_positions = response.get("positions", {})
    if not isinstance(raw_positions, dict):
        return None, {}, "invalid_positions"
    positions: dict[int, tuple[float, float, float]] = {}
    try:
        for atom_id, coords in raw_positions.items():
            if not isinstance(coords, (list, tuple)) or len(coords) != 3:
                continue
            positions[int(atom_id)] = (float(coords[0]), float(coords[1]), float(coords[2]))
    except Exception as exc:
        return None, {}, str(exc)
    if not positions:
        return None, {}, "empty_positions"
    metadata = response.get("metadata", {})
    if not isinstance(metadata, dict):
        metadata = {}
    frames = response.get("frames", [])
    if isinstance(frames, list):
        metadata["frames"] = frames
    return positions, metadata, None


def optimize_3d_isolated(
    graph,
    *,
    coordinates: dict[int, tuple[float, float, float]] | None = None,
    forcefield: str = "MMFF94",
    max_iters: int = 200,
    steps_per_update: int = 25,
    seed: int = 0xC0FFEE,
    timeout_s: float = 20.0,
) -> tuple[dict[int, tuple[float, float, float]] | None, dict[str, Any], str | None]:
    """Optimiza conformación 3D RDKit en un proceso aislado."""
    request = _graph_request_payload(
        graph=graph,
        chain_atom_ids=[],
        mode="graph_optimize3d",
    )
    request.update(
        {
            "forcefield": str(forcefield or "MMFF94"),
            "max_iters": int(max_iters or 200),
            "steps_per_update": int(steps_per_update or 25),
            "seed": int(seed),
            "coordinates": _coordinates_payload(coordinates or {}),
        }
    )
    response = _run_worker(request, timeout_s=timeout_s)
    if not response.get("ok"):
        return None, {}, str(response.get("error", "worker_error"))
    raw_positions = response.get("positions", {})
    if not isinstance(raw_positions, dict):
        return None, {}, "invalid_positions"
    positions: dict[int, tuple[float, float, float]] = {}
    try:
        for atom_id, coords in raw_positions.items():
            if not isinstance(coords, (list, tuple)) or len(coords) != 3:
                continue
            positions[int(atom_id)] = (float(coords[0]), float(coords[1]), float(coords[2]))
    except Exception as exc:
        return None, {}, str(exc)
    if not positions:
        return None, {}, "empty_positions"
    metadata = response.get("metadata", {})
    if not isinstance(metadata, dict):
        metadata = {}
    frames = response.get("frames", [])
    if isinstance(frames, list):
        metadata["frames"] = frames
    return positions, metadata, None


def clean2d_isolated(
    graph,
    timeout_s: float = 8.0,
) -> tuple[dict[int, tuple[float, float]] | None, str | None]:
    """Genera coordenadas 2D limpias usando RDKit en un worker aislado.

    Construye el Mol desde el grafo, llama Compute2DCoords y devuelve
    posiciones 2D mapeadas por atom_id del grafo original.

    Args:
        graph: MolGraph de Chemuson.
        timeout_s: Timeout en segundos.

    Returns:
        Tupla (posiciones, error) donde posiciones mapea atom_id -> (x, y).
    """
    request = _graph_request_payload(
        graph=graph,
        chain_atom_ids=[],
        mode="graph_clean2d",
    )
    response = _run_worker(request, timeout_s=timeout_s)
    if not response.get("ok"):
        return None, str(response.get("error", "worker_error"))
    raw_positions = response.get("positions", {})
    if not isinstance(raw_positions, dict):
        return None, "invalid_positions"
    positions: dict[int, tuple[float, float]] = {}
    try:
        for atom_id, coords in raw_positions.items():
            if not isinstance(coords, (list, tuple)) or len(coords) < 2:
                continue
            positions[int(atom_id)] = (float(coords[0]), float(coords[1]))
    except Exception as exc:
        return None, str(exc)
    if not positions:
        return None, "empty_positions"
    return positions, None


def _coordinates_payload(coordinates: dict[int, tuple[float, float, float]]) -> dict[str, list[float]]:
    payload: dict[str, list[float]] = {}
    for atom_id, coords in coordinates.items():
        if not isinstance(coords, (list, tuple)) or len(coords) != 3:
            continue
        payload[str(int(atom_id))] = [float(coords[0]), float(coords[1]), float(coords[2])]
    return payload


def is_rdkit_worker_available(timeout_s: float = 5.0) -> bool:
    """Smoke-check liviano para saber si el worker RDKit está utilizable."""
    result = run_rdkit_stereo_extract("CC", fmt="smiles", timeout_s=timeout_s)
    return bool(result.get("ok"))


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
