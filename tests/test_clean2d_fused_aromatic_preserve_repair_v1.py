from __future__ import annotations

import math

from chemuson.clean2d import (
    assert_clean2d_invariants,
    capture_clean2d_snapshot,
    count_new_bond_crossings,
    evaluate_clean2d_layout,
    ring_degeneracy_score,
    run_clean2d_engine,
)
from chemuson.core.model import MolGraph


def _collapsed_naphthalene() -> tuple[MolGraph, set[int], tuple[set[int], set[int]]]:
    graph = MolGraph()
    for index in range(10):
        graph.add_atom("C", float((index % 5) * 5), float((index // 5) * 3), atom_id=index + 1)
    for a1, a2 in ((1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1), (4, 7), (7, 8), (8, 9), (9, 10), (10, 5)):
        graph.add_bond(a1, a2, order=1, is_aromatic=True)
    return graph, set(range(1, 11)), ({1, 2, 3, 4, 5, 6}, {4, 5, 7, 8, 9, 10})


def test_simple_fused_aromatic_template_repairs_collapsed_naphthalene() -> None:
    graph, atoms, rings = _collapsed_naphthalene()
    before = {atom_id: (atom.x, atom.y) for atom_id, atom in graph.atoms.items()}
    snapshot = capture_clean2d_snapshot(graph)

    result = run_clean2d_engine(graph, atoms, mode="publication", target_bond_length=40.0)

    assert result.ok and result.selected is not None
    assert result.selected.source == "fused_aromatic_template"
    assert all(ring_degeneracy_score(result.selected.coords, ring) > 0.18 for ring in rings)
    report = evaluate_clean2d_layout(atoms, list(graph.bonds.values()), before, result.selected.coords, 40.0)
    assert report.min_bond_length_ratio >= 0.70
    assert report.max_bond_length_ratio <= 1.45
    assert count_new_bond_crossings(before, result.selected.coords, list(graph.bonds.values())) == 0
    assert report.min_nonbonded_after > 10.0
    assert_clean2d_invariants(graph, graph, before, result.selected.coords, atom_ids=atoms)


def test_fused_template_is_deterministic_and_excludes_substituents() -> None:
    graph, atoms, _rings = _collapsed_naphthalene()
    first = run_clean2d_engine(graph, atoms, mode="publication", target_bond_length=40.0)
    second = run_clean2d_engine(graph, atoms, mode="publication", target_bond_length=40.0)
    assert first.selected is not None and second.selected is not None
    assert first.selected.coords == second.selected.coords

    graph.add_atom("C", 60.0, 0.0)
    graph.add_bond(1, 11, order=1)
    excluded = run_clean2d_engine(graph, mode="publication", target_bond_length=40.0)
    assert all(candidate.source != "fused_aromatic_template" for candidate in (*excluded.candidates, *excluded.rejected))
