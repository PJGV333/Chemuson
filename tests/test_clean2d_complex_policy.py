from __future__ import annotations

import math

from chemuson.clean2d import (
    classify_clean2d_complexity,
    clean2d_candidate_quality_diagnostic,
    run_clean2d_engine,
)
from chemuson.core.model import MolGraph


def _add_ring(graph: MolGraph, cx: float, cy: float, *, n: int = 6, radius: float = 40.0, aromatic: bool = True) -> list[int]:
    atoms = []
    for idx in range(n):
        angle = math.radians(idx * 360.0 / n)
        atoms.append(graph.add_atom("C", cx + radius * math.cos(angle), cy + radius * math.sin(angle)).id)
    for idx in range(n):
        graph.add_bond(atoms[idx], atoms[(idx + 1) % n], order=1, is_aromatic=aromatic)
    return atoms


def _hierarchical_graph() -> MolGraph:
    graph = MolGraph()
    rings = [_add_ring(graph, idx * 120.0, 0.0) for idx in range(4)]
    graph.add_bond(rings[0][0], rings[1][3], order=1)
    graph.add_bond(rings[1][1], rings[2][4], order=1)
    graph.add_bond(rings[2][2], rings[3][5], order=1)
    a = graph.add_atom("C", 80.0, 50.0, stereo_cip="R").id
    b = graph.add_atom("C", 160.0, 50.0, stereo_cip="S").id
    graph.add_bond(rings[0][1], a, order=1)
    graph.add_bond(a, b, order=1)
    graph.add_bond(b, rings[3][4], order=1)
    return graph


def _macrocycle_with_cavity_graph() -> MolGraph:
    graph = MolGraph()
    atoms = _add_ring(graph, 0.0, 0.0, n=12, radius=80.0, aromatic=False)
    bridge_a = graph.add_atom("C", 20.0, 20.0).id
    bridge_b = graph.add_atom("C", 40.0, 20.0).id
    graph.add_bond(atoms[2], bridge_a, order=1)
    graph.add_bond(bridge_a, bridge_b, order=1)
    graph.add_bond(bridge_b, atoms[8], order=1)
    return graph


def _cyclophane_like_graph() -> MolGraph:
    graph = MolGraph()
    left = _add_ring(graph, -80.0, 0.0)
    right = _add_ring(graph, 80.0, 0.0)
    top_a = graph.add_atom("C", -20.0, 70.0).id
    top_b = graph.add_atom("C", 20.0, 70.0).id
    bottom_a = graph.add_atom("C", -20.0, -70.0).id
    bottom_b = graph.add_atom("C", 20.0, -70.0).id
    graph.add_bond(left[1], top_a, order=1)
    graph.add_bond(top_a, top_b, order=1)
    graph.add_bond(top_b, right[2], order=1)
    graph.add_bond(left[4], bottom_a, order=1)
    graph.add_bond(bottom_a, bottom_b, order=1)
    graph.add_bond(bottom_b, right[5], order=1)
    return graph


def _simple_graph() -> MolGraph:
    graph = MolGraph()
    a = graph.add_atom("C", 0.0, 0.0).id
    b = graph.add_atom("C", 40.0, 0.0).id
    c = graph.add_atom("C", 80.0, 0.0).id
    graph.add_bond(a, b, order=1)
    graph.add_bond(b, c, order=1)
    return graph


def test_hierarchical_blocks_policy_blocks_global_redraw_and_local_repair() -> None:
    profile = classify_clean2d_complexity(_hierarchical_graph())

    assert profile.has_hierarchical_blocks
    assert profile.preserve_only_required
    assert profile.global_redraw_allowed is False
    assert profile.local_repair_allowed is False
    assert profile.internal_route_allowed is True
    assert profile.reason


def test_macrocycle_internal_cavity_and_bridge_policy_blocks_global_redraw() -> None:
    profile = classify_clean2d_complexity(_macrocycle_with_cavity_graph())

    assert profile.macrocycle_count >= 1
    assert profile.intramolecular_bridge_count >= 1
    assert profile.internal_cavity_count >= 1
    assert profile.global_redraw_allowed is False
    assert profile.local_repair_allowed is False


def test_cyclophane_like_policy_blocks_global_redraw_or_records_high_risk_bridge() -> None:
    profile = classify_clean2d_complexity(_cyclophane_like_graph())

    assert profile.cyclophane_count >= 1 or profile.intramolecular_bridge_count >= 1 or profile.has_hierarchical_blocks
    assert profile.global_redraw_allowed is False


def test_simple_structure_is_not_forced_to_preserve_only() -> None:
    profile = classify_clean2d_complexity(_simple_graph())

    result = run_clean2d_engine(_simple_graph(), mode="quick", target_bond_length=40.0)

    assert profile.preserve_only_required is False
    assert profile.global_redraw_allowed is True
    assert profile.local_repair_allowed is False
    assert result.selected is not None
    assert result.selected.metadata.get("complex_preserve_only") is not True


def test_complex_policy_evidence_is_json_safe_in_diagnostics_and_snapshots() -> None:
    result = run_clean2d_engine(_hierarchical_graph(), mode="quick", target_bond_length=40.0, debug_snapshot=True)

    assert result.selected is not None
    diagnostic = clean2d_candidate_quality_diagnostic(result.selected)
    policy = diagnostic["internal"]["metadata"]["complex_policy"]
    assert policy["global_redraw_allowed"] is False
    assert policy["local_repair_allowed"] is False
    assert result.debug_snapshot is not None
    selected = next(row for row in result.debug_snapshot["candidates"] if row["selected"])
    snapshot_policy = selected["quality_diagnostic"]["internal"]["metadata"]["complex_policy"]
    assert snapshot_policy == policy
