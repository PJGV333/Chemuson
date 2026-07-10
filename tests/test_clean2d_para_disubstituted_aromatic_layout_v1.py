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
from chemuson.core.model import BondStyle, MolGraph


def _para_benzene() -> tuple[MolGraph, set[int], tuple[int, int], tuple[int, int]]:
    graph = MolGraph()
    ring = [graph.add_atom("C", float(index * 7), float((index % 2) * 3)).id for index in range(6)]
    for index in range(6):
        graph.add_bond(ring[index], ring[(index + 1) % 6], order=1, is_aromatic=True)
    left = graph.add_atom("Cl", -5.0, 2.0).id
    right = graph.add_atom("O", 40.0, -2.0).id
    graph.add_bond(ring[0], left, order=1)
    graph.add_bond(ring[3], right, order=1)
    return graph, set(ring), (ring[0], ring[3]), (left, right)


def _sources(result) -> set[str]:
    return {candidate.source for candidate in (*result.candidates, *result.rejected)}


def test_para_template_regularizes_opposite_terminal_substituents() -> None:
    graph, ring, anchors, substituents = _para_benzene()
    atoms = set(graph.atoms)
    before = {atom_id: (atom.x, atom.y) for atom_id, atom in graph.atoms.items()}
    snapshot = capture_clean2d_snapshot(graph)
    result = run_clean2d_engine(graph, mode="publication", target_bond_length=40.0)

    assert result.ok and result.selected is not None
    assert result.selected.source == "para_disubstituted_aromatic_template"
    coords = result.selected.coords
    center = (sum(coords[atom_id][0] for atom_id in ring) / 6, sum(coords[atom_id][1] for atom_id in ring) / 6)
    vectors = []
    for anchor, substituent in zip(anchors, substituents):
        outward = (coords[anchor][0] - center[0], coords[anchor][1] - center[1])
        branch = (coords[substituent][0] - coords[anchor][0], coords[substituent][1] - coords[anchor][1])
        assert outward[0] * branch[0] + outward[1] * branch[1] > 0
        vectors.append(branch)
    assert vectors[0][0] * vectors[1][0] + vectors[0][1] * vectors[1][1] < 0
    report = evaluate_clean2d_layout(atoms, list(graph.bonds.values()), before, coords, 40.0)
    assert ring_degeneracy_score(coords, ring) > 0.18
    assert report.min_bond_length_ratio >= 0.70 and report.max_bond_length_ratio <= 1.45
    assert report.min_nonbonded_after > 10.0
    assert count_new_bond_crossings(before, coords, list(graph.bonds.values())) == 0
    assert_clean2d_invariants(snapshot, graph, before, coords, atom_ids=atoms)


def test_para_template_is_deterministic_and_excludes_non_para_or_stereo() -> None:
    graph, _ring, _anchors, _substituents = _para_benzene()
    first = run_clean2d_engine(graph, mode="quick", target_bond_length=40.0)
    second = run_clean2d_engine(graph, mode="quick", target_bond_length=40.0)
    assert first.selected is not None and second.selected is not None
    assert first.selected.coords == second.selected.coords

    graph, ring, _anchors, _substituents = _para_benzene()
    graph.bonds[7].a1_id = sorted(ring)[1]
    assert "para_disubstituted_aromatic_template" not in _sources(run_clean2d_engine(graph, mode="quick", target_bond_length=40.0))
    graph, _ring, _anchors, _substituents = _para_benzene()
    graph.bonds[1].style = BondStyle.WEDGE
    assert "para_disubstituted_aromatic_template" not in _sources(run_clean2d_engine(graph, mode="quick", target_bond_length=40.0))
    graph, _ring, _anchors, _substituents = _para_benzene()
    assert "para_disubstituted_aromatic_template" not in _sources(run_clean2d_engine(graph, mode="propose", target_bond_length=40.0))
