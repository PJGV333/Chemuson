from __future__ import annotations

"""Opt-in JSON debug snapshots for Clean 2D runs.

The helpers in this module are reporting-only. They serialize data already
present in a Clean 2D result and never generate, rank, or select candidates.
"""

from collections.abc import Iterable, Mapping
import json
import math
import os
from pathlib import Path
from typing import Any, TYPE_CHECKING

from chemuson.clean2d.quality_reporting import clean2d_candidate_quality_diagnostic
from chemuson.core.model import MolGraph, bond_affects_valence

if TYPE_CHECKING:
    from chemuson.clean2d.engine import Clean2DResult


CLEAN2D_DEBUG_SNAPSHOT_SCHEMA = "chemuson.clean2d.debug-snapshot"
CLEAN2D_DEBUG_SNAPSHOT_VERSION = 1
CLEAN2D_DEBUG_SNAPSHOT_ENV = "CHEMUSON_CLEAN2D_DEBUG_SNAPSHOT"


def clean2d_debug_snapshot_enabled(debug_snapshot: bool | None = None) -> bool:
    """Return whether Clean 2D debug snapshot capture is explicitly enabled."""

    if debug_snapshot is not None:
        return bool(debug_snapshot)
    return str(os.environ.get(CLEAN2D_DEBUG_SNAPSHOT_ENV, "")).strip().lower() in {
        "1",
        "true",
        "yes",
        "on",
    }


def build_clean2d_debug_snapshot(
    graph: MolGraph,
    result: "Clean2DResult",
    *,
    target_whole_structure: bool,
    initial_selection: Mapping[str, Iterable[int]] | None = None,
    metadata: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Build a stable JSON-serializable snapshot from an existing result."""

    atom_ids = sorted(int(atom_id) for atom_id in result.atom_ids)
    final_coords = None
    if result.selected is not None and not result.selected.rejected:
        final_coords = _coords_json(result.selected.coords, atom_ids)

    snapshot = {
        "schema": CLEAN2D_DEBUG_SNAPSHOT_SCHEMA,
        "version": CLEAN2D_DEBUG_SNAPSHOT_VERSION,
        "mode": getattr(result.mode, "value", str(result.mode)),
        "target": {
            "whole_structure": bool(target_whole_structure),
            "atom_ids": None if target_whole_structure else atom_ids,
        },
        "initial_selection": _selection_json(initial_selection),
        "topology": _topology_json(graph, atom_ids),
        "initial_coordinates": _coords_json(result.before_coords, atom_ids),
        "final_coordinates": final_coords,
        "candidates": _candidates_json(result),
        "final_state": str(result.result_state),
        "final_reason": str(result.stable_reason or "") or None,
        "metadata": _json_safe(dict(metadata or {})),
    }
    validate_clean2d_debug_snapshot(snapshot)
    return snapshot


def write_clean2d_debug_snapshot(path: str | Path, snapshot: Mapping[str, Any]) -> Path:
    """Write a snapshot to an explicit caller-controlled JSON path."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    data = validate_clean2d_debug_snapshot(snapshot)
    output.write_text(
        json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return output


def read_clean2d_debug_snapshot(path: str | Path) -> dict[str, Any]:
    """Read and validate a Clean 2D debug snapshot JSON file."""

    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return validate_clean2d_debug_snapshot(data)


def validate_clean2d_debug_snapshot(snapshot: Mapping[str, Any]) -> dict[str, Any]:
    """Return a JSON-safe snapshot or raise if required stable fields are missing."""

    data = _json_safe(dict(snapshot))
    if data.get("schema") != CLEAN2D_DEBUG_SNAPSHOT_SCHEMA:
        raise ValueError("invalid Clean 2D debug snapshot schema")
    if data.get("version") != CLEAN2D_DEBUG_SNAPSHOT_VERSION:
        raise ValueError("invalid Clean 2D debug snapshot version")
    for key in (
        "mode",
        "target",
        "initial_selection",
        "topology",
        "initial_coordinates",
        "candidates",
        "final_state",
    ):
        if key not in data:
            raise ValueError(f"missing Clean 2D debug snapshot field: {key}")
    json.dumps(data, allow_nan=False)
    return data


def run_clean2d_with_debug_snapshot(
    graph: MolGraph,
    atom_ids: Iterable[int] | None = None,
    *,
    path: str | Path | None = None,
    initial_selection: Mapping[str, Iterable[int]] | None = None,
    metadata: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> "Clean2DResult":
    """Test helper that scopes snapshot capture to an explicit call/path."""

    from chemuson.clean2d.engine import run_clean2d_engine

    return run_clean2d_engine(
        graph,
        atom_ids,
        debug_snapshot=True,
        debug_snapshot_path=path,
        debug_initial_selection=initial_selection,
        debug_snapshot_metadata=metadata,
        **kwargs,
    )


def _selection_json(selection: Mapping[str, Iterable[int]] | None) -> dict[str, list[int]]:
    selection = selection or {}
    return {
        "atom_ids": sorted(int(atom_id) for atom_id in selection.get("atom_ids", ()) or ()),
        "bond_ids": sorted(int(bond_id) for bond_id in selection.get("bond_ids", ()) or ()),
    }


def _topology_json(graph: MolGraph, atom_ids: list[int]) -> dict[str, Any]:
    selected = set(atom_ids)
    atoms = [
        {
            "id": int(atom_id),
            "element": str(graph.atoms[atom_id].element),
        }
        for atom_id in atom_ids
        if atom_id in graph.atoms
    ]
    bonds = []
    for bond in sorted(graph.bonds.values(), key=lambda item: int(item.id)):
        if bond.a1_id not in selected or bond.a2_id not in selected:
            continue
        if not bond_affects_valence(bond):
            continue
        bonds.append(
            {
                "id": int(bond.id),
                "atom_ids": [int(bond.a1_id), int(bond.a2_id)],
                "order": int(getattr(bond, "order", 1) or 1),
                "is_aromatic": bool(getattr(bond, "is_aromatic", False)),
            }
        )
    return {"atoms": atoms, "bonds": bonds}


def _coords_json(coords: Mapping[int, tuple[float, float]], atom_ids: Iterable[int]) -> list[dict[str, Any]]:
    out = []
    for atom_id in sorted(int(item) for item in atom_ids):
        if atom_id not in coords:
            continue
        x, y = coords[atom_id]
        out.append({"atom_id": atom_id, "x": _finite_float(x), "y": _finite_float(y)})
    return out


def _candidates_json(result: "Clean2DResult") -> list[dict[str, Any]]:
    selected = result.selected
    rows = []
    for index, candidate in enumerate((*result.candidates, *result.rejected)):
        row = {
            "index": index,
            "source": str(candidate.source),
            "rejected": bool(candidate.rejected),
            "selected": candidate is selected,
            "quality_diagnostic": clean2d_candidate_quality_diagnostic(candidate),
        }
        rows.append(_json_safe(row))
    return rows


def _finite_float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _json_safe(value: Any) -> Any:
    if value is None or isinstance(value, (str, bool)):
        return value
    if isinstance(value, int) and not isinstance(value, bool):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        items = value
        if isinstance(value, (set, frozenset)):
            items = sorted(value)
        return [_json_safe(item) for item in items]
    enum_value = getattr(value, "value", None)
    if isinstance(enum_value, (str, int, float, bool)) or enum_value is None:
        return _json_safe(enum_value)
    return None
