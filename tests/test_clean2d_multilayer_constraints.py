from __future__ import annotations

import math


from chemuson.clean2d import (
    assert_clean2d_invariants,
    capture_clean2d_snapshot,
    classify_clean2d_layout_quality,
    generate_clean2d_candidates,
    has_hierarchical_block_layout_signals,
    has_intramolecular_block_layout_problem,
    run_clean2d_engine,
)
import chemuson.clean2d.engine as clean2d_engine
from chemuson.core.layers import build_multilayer_chemical_graph
from chemuson.core.model import BondStyle, MolGraph


def _distance(a: tuple[float, float], b: tuple[float, float]) -> float:
    return math.hypot(b[0] - a[0], b[1] - a[1])


def test_clean2d_moves_motif_blocks_for_coordination_without_valence_change() -> None:
    graph = MolGraph()
    metal = graph.add_atom("Fe", 0.0, 0.0, is_coordination_center=True).id
    n = graph.add_atom("N", 200.0, 0.0).id
    c = graph.add_atom("C", 240.0, 0.0).id
    graph.add_bond(n, c, order=1)
    graph.add_bond(metal, n, style=BondStyle.COORDINATION, donor_atom_id=n, length_px=48.0)

    before_valence = graph.bond_order_sum(metal)
    result = run_clean2d_engine(graph, mode="publication", target_bond_length=40.0)

    assert result.ok
    assert result.selected is not None
    assert result.selected.source == "motif_constraints"
    assert graph.bond_order_sum(metal) == before_valence == 0.0
    after = result.selected.coords
    assert _distance(after[metal], after[n]) < 80.0
    assert math.isclose(_distance(after[n], after[c]), 40.0, rel_tol=0.25)


def _add_hexagon(graph: MolGraph, cx: float, cy: float, radius: float = 24.0) -> list[int]:
    atoms = []
    for idx in range(6):
        angle = math.radians(60.0 * idx)
        atoms.append(graph.add_atom("C", cx + math.cos(angle) * radius, cy + math.sin(angle) * radius).id)
    for idx, atom_id in enumerate(atoms):
        graph.add_bond(atom_id, atoms[(idx + 1) % len(atoms)], is_aromatic=True)
    return atoms


def test_clean2d_generates_internal_block_transform_across_linker() -> None:
    graph = MolGraph()
    left_ring = _add_hexagon(graph, 0.0, 0.0)
    right_ring = _add_hexagon(graph, 180.0, 0.0)
    linker_a = graph.add_atom("C", 54.0, 0.0).id
    linker_b = graph.add_atom("C", 126.0, 0.0).id
    graph.add_bond(left_ring[0], linker_a, order=1)
    graph.add_bond(linker_a, linker_b, order=1)
    graph.add_bond(linker_b, right_ring[3], order=1)
    graph.add_bond(left_ring[2], right_ring[1], style=BondStyle.INTERACTION, length_px=70.0)

    before = {
        atom_id: (atom.x, atom.y)
        for atom_id, atom in graph.atoms.items()
    }
    before_distance = _distance(before[left_ring[2]], before[right_ring[1]])
    candidates = generate_clean2d_candidates(graph, mode="publication", target_bond_length=40.0)
    block_candidate = next(
        candidate for candidate in candidates if candidate.source in {"block_layout", "block_constraints"}
    )

    assert block_candidate.metadata["block_operation_count"] >= 1
    assert block_candidate.metadata["interaction_constraint_error"] < block_candidate.metadata["interaction_constraint_error_before"]
    assert _distance(block_candidate.coords[left_ring[2]], block_candidate.coords[right_ring[1]]) < before_distance


def test_clean2d_uses_blockgraph_for_three_aromatic_rings_without_interactions() -> None:
    graph = MolGraph()
    rings = [_add_hexagon(graph, float(idx) * 70.0, 0.0) for idx in range(3)]
    for idx in range(2):
        linker = graph.add_atom("C", float(idx) * 70.0 + 35.0, 0.0).id
        graph.add_bond(rings[idx][0], linker, order=1)
        graph.add_bond(linker, rings[idx + 1][3], order=1)

    layers = build_multilayer_chemical_graph(graph)
    quality = classify_clean2d_layout_quality(graph, target_bond_length=40.0)
    assert has_intramolecular_block_layout_problem(layers.block_graph, quality)
    assert has_hierarchical_block_layout_signals(layers.block_graph)

    candidates = generate_clean2d_candidates(graph, mode="publication", target_bond_length=40.0)
    block_candidate = next(candidate for candidate in candidates if candidate.source == "block_layout")

    assert block_candidate.metadata["interaction_constraint_count"] == 0
    assert block_candidate.metadata["internal_block_constraint_count"] > 0
    assert block_candidate.metadata["block_operation_count"] >= 1
    assert block_candidate.metadata["block_constraint_error"] < block_candidate.metadata["block_constraint_error_before"]


def test_block_problem_does_not_route_to_local_graph_cleaner_first(monkeypatch) -> None:
    graph = MolGraph()
    atoms = [graph.add_atom("C", float(idx), 0.0).id for idx in range(12)]
    for idx, atom_id in enumerate(atoms):
        graph.add_bond(atom_id, atoms[(idx + 1) % len(atoms)], order=1)

    def fail_local_graph(*args, **kwargs):
        raise AssertionError("local_graph_cleaner_should_not_run_first")

    monkeypatch.setattr(clean2d_engine, "local_graph_clean2d", fail_local_graph)
    result = run_clean2d_engine(graph, mode="quick", target_bond_length=40.0)

    assert result.selected is None or result.selected.source != "local_graph"


def test_tetrandrine_like_hierarchical_blocks_do_not_select_local_graph() -> None:
    graph = MolGraph()
    rings = [
        _add_hexagon(graph, 0.0, 0.0),
        _add_hexagon(graph, 64.0, 0.0),
        _add_hexagon(graph, 64.0, 52.0),
        _add_hexagon(graph, 0.0, 52.0),
    ]
    stereo_a = graph.add_atom("C", 32.0, 20.0, stereo_cip="R").id
    stereo_b = graph.add_atom("C", 32.0, 32.0, stereo_cip="S").id
    graph.add_bond(rings[0][0], stereo_a, order=1)
    graph.add_bond(stereo_a, rings[1][3], order=1)
    graph.add_bond(rings[2][0], stereo_b, order=1)
    graph.add_bond(stereo_b, rings[3][3], order=1)

    bridge_a = graph.add_atom("C", 32.0, -18.0).id
    bridge_b = graph.add_atom("C", 96.0, 26.0).id
    bridge_c = graph.add_atom("C", 32.0, 70.0).id
    bridge_d = graph.add_atom("C", -32.0, 26.0).id
    graph.add_bond(rings[0][1], bridge_a, order=1)
    graph.add_bond(bridge_a, rings[1][5], order=1)
    graph.add_bond(rings[1][1], bridge_b, order=1)
    graph.add_bond(bridge_b, rings[2][5], order=1)
    graph.add_bond(rings[2][1], bridge_c, order=1)
    graph.add_bond(bridge_c, rings[3][5], order=1)
    graph.add_bond(rings[3][1], bridge_d, order=1)
    graph.add_bond(bridge_d, rings[0][5], order=1)

    before = {atom_id: (atom.x, atom.y) for atom_id, atom in graph.atoms.items()}
    bonds = [bond for bond in graph.bonds.values() if clean2d_engine.bond_affects_valence(bond)]
    before_crossings = clean2d_engine._count_crossings(before, bonds)
    before_quality = classify_clean2d_layout_quality(graph, coords=before, target_bond_length=40.0)
    snapshot = capture_clean2d_snapshot(graph)
    layers = build_multilayer_chemical_graph(graph)

    assert has_hierarchical_block_layout_signals(layers.block_graph)
    assert not any(bond.style in {BondStyle.INTERACTION, BondStyle.COORDINATION} for bond in graph.bonds.values())

    result = run_clean2d_engine(graph, mode="publication", target_bond_length=40.0)

    assert result.selected is not None
    assert result.selected.source in {"block_layout", "block_constraints"}
    after = result.selected.coords
    assert clean2d_engine._count_crossings(after, bonds) <= before_crossings
    quality = classify_clean2d_layout_quality(graph, coords=after, target_bond_length=40.0)
    assert quality.min_ring_degeneracy + 1e-6 >= before_quality.min_ring_degeneracy
    assert_clean2d_invariants(snapshot, graph, before, after)
