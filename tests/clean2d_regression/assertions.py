from __future__ import annotations

import json
import math
from dataclasses import asdict, is_dataclass
from typing import Any

from chemuson.clean2d import (
    assert_clean2d_invariants,
    capture_clean2d_snapshot,
    classify_clean2d_layout_quality,
    run_clean2d_with_debug_snapshot,
    validate_clean2d_debug_snapshot,
)
from chemuson.core.model import MolGraph

from .cases import Clean2DRegressionCase


def assert_case_metadata(case: Clean2DRegressionCase) -> None:
    assert case.name.strip(), "regression case name is required"
    assert case.family.strip(), f"{case.name}: family is required"
    assert case.mode, f"{case.name}: Clean 2D mode is required"
    assert case.target_label in {"whole", "selection"}
    assert case.expected_states, f"{case.name}: expected contract states are required"


def assert_valid_graph(graph: MolGraph, case: Clean2DRegressionCase) -> None:
    assert graph.atoms, f"{case.name}: graph must contain atoms"
    assert graph.bonds, f"{case.name}: graph must contain bonds"
    for bond in graph.bonds.values():
        assert bond.a1_id in graph.atoms, f"{case.name}: bond {bond.id} has missing a1"
        assert bond.a2_id in graph.atoms, f"{case.name}: bond {bond.id} has missing a2"
    if case.target is not None:
        missing = set(case.target) - set(graph.atoms)
        assert not missing, f"{case.name}: target atoms missing: {sorted(missing)}"


def execute_case(case: Clean2DRegressionCase) -> dict[str, Any]:
    assert_case_metadata(case)
    graph = case.builder()
    assert_valid_graph(graph, case)

    atom_ids = None if case.target is None else set(case.target)
    before_snapshot = capture_clean2d_snapshot(graph, atom_ids)
    before_metrics = _collect_metrics(graph, atom_ids)

    result = run_clean2d_with_debug_snapshot(
        graph,
        atom_ids,
        mode=case.mode,
        target_bond_length=40.0,
        metadata={
            "case": case.name,
            "family": case.family,
            "known_delicate": case.known_delicate,
            "known_failure": case.known_failure,
        },
    )

    assert result.result_state in case.expected_states, (
        f"{case.name}: unexpected state {result.result_state}; "
        f"expected one of {sorted(case.expected_states)}"
    )
    if result.selected is not None and not result.selected.rejected:
        assert_clean2d_invariants(before_snapshot, graph, result.before_coords, result.selected.coords, atom_ids=atom_ids)

    snapshot = result.debug_snapshot
    assert snapshot is not None, f"{case.name}: debug snapshot was not produced"
    validate_clean2d_debug_snapshot(snapshot)
    json.dumps(snapshot, allow_nan=False)

    after_metrics = None
    if result.selected is not None and not result.selected.rejected:
        after_metrics = _collect_metrics(graph, atom_ids, coords=result.selected.coords)

    corpus_snapshot = {
        "case": {
            "name": case.name,
            "family": case.family,
            "mode": case.mode.value,
            "target": case.target_label,
            "target_atom_ids": None if case.target is None else list(case.target),
            "expected_states": sorted(case.expected_states),
            "known_delicate": case.known_delicate,
            "known_failure": case.known_failure,
            "notes": case.notes,
        },
        "result": {
            "state": result.result_state,
            "reason": result.stable_reason or None,
            "candidate_sources": list(result.candidate_sources),
        },
        "identity": {
            "atom_ids": sorted(before_snapshot.atom_ids),
            "bond_ids": sorted(before_snapshot.bond_ids),
        },
        "metrics": {
            "before": before_metrics,
            "after": after_metrics,
        },
        "engine_debug_snapshot": snapshot,
    }
    json.dumps(corpus_snapshot, allow_nan=False)
    return corpus_snapshot


def _collect_metrics(
    graph: MolGraph,
    atom_ids: set[int] | None,
    *,
    coords: dict[int, tuple[float, float]] | None = None,
) -> dict[str, Any]:
    report = classify_clean2d_layout_quality(graph, atom_ids, coords=coords, target_bond_length=40.0)
    return _json_safe(
        {
            "quality_class": report.quality_class,
            "reason": report.reason,
            "crossings": report.crossings,
            "min_nonbonded_distance": report.min_nonbonded_distance,
            "min_ring_degeneracy": report.min_ring_degeneracy,
            "length_rms_error": report.length_rms_error,
            "length_max_error": report.length_max_error,
            "angle_rms_deviation": report.angle_rms_deviation,
            "angle_max_deviation": report.angle_max_deviation,
            "visual_score": report.visual_score,
        }
    )


def _json_safe(value: Any) -> Any:
    if is_dataclass(value):
        return _json_safe(asdict(value))
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_json_safe(item) for item in value]
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    return value
