from __future__ import annotations

import math
import os
import sys

import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.clean2d import count_new_bond_crossings, ring_degeneracy_score, run_clean2d_engine
from chemuson.core.model import MolGraph


def _scrambled_benzene(substituent: str | None = None) -> tuple[MolGraph, list[int], int | None]:
    graph = MolGraph()
    coords = [(0, 0), (12, 5), (28, -4), (8, 2), (-8, -3), (-24, 6)]
    ring = [
        graph.add_atom("C", float(x), float(y)).id
        for x, y in coords
    ]
    for idx in range(6):
        graph.add_bond(ring[idx], ring[(idx + 1) % 6], order=1, is_aromatic=True)
    subst_id = None
    if substituent is not None:
        subst = graph.add_atom(substituent, 1.0, 1.0)
        subst_id = subst.id
        graph.add_bond(ring[0], subst.id, order=1)
    return graph, ring, subst_id


def _regular_benzene_with_substituent(element: str = "C") -> tuple[MolGraph, list[int], int]:
    graph = MolGraph()
    ring = []
    radius = 40.0
    for idx in range(6):
        angle = math.radians(60.0 * idx + 30.0)
        ring.append(graph.add_atom("C", radius * math.cos(angle), radius * math.sin(angle)).id)
    for idx in range(6):
        graph.add_bond(ring[idx], ring[(idx + 1) % 6], order=1, is_aromatic=True)
    anchor = ring[0]
    ax, ay = graph.atoms[anchor].x, graph.atoms[anchor].y
    subst = graph.add_atom(element, ax * 2.0, ay * 2.0)
    graph.add_bond(anchor, subst.id, order=1)
    return graph, ring, subst.id


def _naphthalene_scrambled() -> tuple[MolGraph, list[int]]:
    graph = MolGraph()
    for idx in range(10):
        graph.add_atom("C", float((idx % 5) * 5), float((idx // 5) * 3), atom_id=idx + 1)
    edges = [
        (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1),
        (4, 7), (7, 8), (8, 9), (9, 10), (10, 5),
    ]
    for a1, a2 in edges:
        graph.add_bond(a1, a2, order=1, is_aromatic=True)
    return graph, list(range(1, 11))


def _external_vector_score(
    coords: dict[int, tuple[float, float]],
    ring: list[int],
    anchor: int,
    substituent: int,
) -> float:
    cx = sum(coords[atom_id][0] for atom_id in ring) / len(ring)
    cy = sum(coords[atom_id][1] for atom_id in ring) / len(ring)
    ax, ay = coords[anchor]
    sx, sy = coords[substituent]
    outward = (ax - cx, ay - cy)
    branch = (sx - ax, sy - ay)
    return outward[0] * branch[0] + outward[1] * branch[1]


def _min_nonbonded_distance(
    coords: dict[int, tuple[float, float]],
    graph: MolGraph,
) -> float:
    bonded = {
        (min(b.a1_id, b.a2_id), max(b.a1_id, b.a2_id))
        for b in graph.bonds.values()
    }
    ids = sorted(coords)
    best = float("inf")
    for idx, a_id in enumerate(ids):
        for b_id in ids[idx + 1:]:
            if (min(a_id, b_id), max(a_id, b_id)) in bonded:
                continue
            best = min(best, math.hypot(coords[a_id][0] - coords[b_id][0], coords[a_id][1] - coords[b_id][1]))
    return best


def test_benzene_substituent_preserves_aromaticity_and_points_outward() -> None:
    graph, ring, methyl_id = _scrambled_benzene("C")

    result = run_clean2d_engine(graph, mode="publication", target_bond_length=40.0)

    assert result.ok
    assert result.selected is not None
    assert all(bond.is_aromatic for bond in graph.bonds.values() if {bond.a1_id, bond.a2_id} <= set(ring))
    assert ring_degeneracy_score(result.selected.coords, set(ring)) > 0.25
    assert _external_vector_score(result.selected.coords, ring, ring[0], methyl_id) > 0.0


@pytest.mark.parametrize("element", ["N", "O"])
def test_aniline_and_phenol_like_substituents_do_not_collide(element: str) -> None:
    graph, ring, subst_id = _scrambled_benzene(element)

    result = run_clean2d_engine(graph, mode="publication", target_bond_length=40.0)

    assert result.ok
    assert result.selected is not None
    assert subst_id in result.selected.coords
    assert _min_nonbonded_distance(result.selected.coords, graph) > 12.0


def test_anisole_like_two_atom_substituent_is_laid_out_from_ring() -> None:
    graph, ring, oxygen_id = _scrambled_benzene("O")
    oxygen = graph.atoms[oxygen_id]
    methyl = graph.add_atom("C", oxygen.x + 1.0, oxygen.y + 1.0)
    graph.add_bond(oxygen_id, methyl.id, order=1)

    result = run_clean2d_engine(graph, mode="publication", target_bond_length=40.0)

    assert result.ok
    assert result.selected is not None
    assert _external_vector_score(result.selected.coords, ring, ring[0], oxygen_id) > 0.0
    assert _min_nonbonded_distance(result.selected.coords, graph) > 12.0


def test_styrene_preserves_external_double_bond_without_new_crossing() -> None:
    graph, ring, vinyl_id = _regular_benzene_with_substituent("C")
    terminal = graph.add_atom("C", graph.atoms[vinyl_id].x + 35.0, graph.atoms[vinyl_id].y)
    graph.add_bond(vinyl_id, terminal.id, order=2)
    before = {atom_id: (atom.x, atom.y) for atom_id, atom in graph.atoms.items()}

    result = run_clean2d_engine(graph, mode="publication", target_bond_length=40.0)

    assert result.ok
    assert result.selected is not None
    assert count_new_bond_crossings(before, result.selected.coords, list(graph.bonds.values())) == 0
    double_len = math.hypot(
        result.selected.coords[terminal.id][0] - result.selected.coords[vinyl_id][0],
        result.selected.coords[terminal.id][1] - result.selected.coords[vinyl_id][1],
    )
    assert double_len == pytest.approx(40.0 * 0.97, rel=0.12)


def test_naphthalene_fused_system_does_not_collapse() -> None:
    graph, atoms = _naphthalene_scrambled()

    result = run_clean2d_engine(graph, atoms, mode="publication", target_bond_length=40.0)

    assert result.ok
    assert result.selected is not None
    first_ring = {1, 2, 3, 4, 5, 6}
    second_ring = {4, 5, 7, 8, 9, 10}
    assert ring_degeneracy_score(result.selected.coords, first_ring) > 0.18
    assert ring_degeneracy_score(result.selected.coords, second_ring) > 0.18


@pytest.mark.parametrize(
    "elements",
    [
        ["N", "C", "C", "C", "C", "C"],  # pyridine-like
        ["N", "C", "N", "C", "C"],        # imidazole-like
    ],
)
def test_simple_heteroaromatics_keep_geometry_and_aromatic_bonds(elements: list[str]) -> None:
    graph = MolGraph()
    ids = []
    for idx, element in enumerate(elements):
        ids.append(graph.add_atom(element, float(idx * 3), float((idx % 2) * 2)).id)
    for idx in range(len(ids)):
        graph.add_bond(ids[idx], ids[(idx + 1) % len(ids)], order=1, is_aromatic=True)

    result = run_clean2d_engine(graph, ids, mode="publication", target_bond_length=40.0)

    assert result.ok
    assert result.selected is not None
    assert all(bond.is_aromatic for bond in graph.bonds.values())
    assert ring_degeneracy_score(result.selected.coords, set(ids)) > 0.18


def test_cyclohexane_aliphatic_ring_is_not_treated_as_aromatic() -> None:
    graph = MolGraph()
    ids = []
    for idx in range(6):
        ids.append(graph.add_atom("C", float(idx * 4), float((idx % 2) * 2)).id)
    for idx in range(6):
        graph.add_bond(ids[idx], ids[(idx + 1) % 6], order=1, is_aromatic=False)

    result = run_clean2d_engine(graph, ids, mode="publication", target_bond_length=40.0)

    assert result.ok
    assert result.selected is not None
    assert not any(bond.is_aromatic for bond in graph.bonds.values())
    assert ring_degeneracy_score(result.selected.coords, set(ids)) > 0.18
