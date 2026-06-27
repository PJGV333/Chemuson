from __future__ import annotations

import math

import pytest


from chemuson.clean2d import (
    assert_clean2d_invariants,
    clean2d_geometry_hash,
    count_new_bond_crossings,
    capture_clean2d_snapshot,
    ring_degeneracy_score,
    run_clean2d_engine,
)
from chemuson.core.model import BondStyle, MolGraph


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


def _apply_coords(graph: MolGraph, coords: dict[int, tuple[float, float]]) -> None:
    for atom_id, (x, y) in coords.items():
        graph.atoms[atom_id].x = float(x)
        graph.atoms[atom_id].y = float(y)


def _mixed_chain_phenyl_cyclopropane() -> tuple[MolGraph, set[int], set[int]]:
    graph = MolGraph()
    phenyl: list[int] = []
    radius = 40.0
    for idx in range(6):
        angle = math.radians(60.0 * idx + 30.0)
        phenyl.append(graph.add_atom("C", radius * math.cos(angle), radius * math.sin(angle)).id)
    for idx in range(6):
        graph.add_bond(phenyl[idx], phenyl[(idx + 1) % 6], order=1, is_aromatic=True)

    chain_coords = [
        (90.0, 20.0),
        (130.0, 20.0),
        (150.0, 54.6),
        (190.0, 54.6),
    ]
    chain = [graph.add_atom("C", x, y).id for x, y in chain_coords]
    graph.add_bond(phenyl[0], chain[0], order=1)
    for left, right in zip(chain, chain[1:]):
        graph.add_bond(left, right, order=1)
    branch = graph.add_atom("C", 130.0, -20.0)
    graph.add_bond(chain[1], branch.id, order=1)

    cyclopropane = [
        graph.add_atom("C", 230.0, 54.6).id,
        graph.add_atom("C", 270.0, 54.6).id,
        graph.add_atom("C", 250.0, 89.2).id,
    ]
    graph.add_bond(chain[-1], cyclopropane[0], order=1)
    graph.add_bond(cyclopropane[0], cyclopropane[1], order=1)
    graph.add_bond(cyclopropane[1], cyclopropane[2], order=1)
    graph.add_bond(cyclopropane[2], cyclopropane[0], order=1)
    return graph, set(phenyl), set(cyclopropane)


def _distort_mixed_graph(graph: MolGraph) -> dict[int, tuple[float, float]]:
    before = {atom_id: (atom.x, atom.y) for atom_id, atom in graph.atoms.items()}
    shifts = {
        7: (50.0, -20.0),
        8: (80.0, 30.0),
        9: (100.0, -40.0),
        10: (130.0, 60.0),
        11: (70.0, -50.0),
        12: (150.0, 20.0),
        13: (170.0, 60.0),
        14: (140.0, 80.0),
    }
    for atom_id, coord in shifts.items():
        graph.atoms[atom_id].x += coord[0]
        graph.atoms[atom_id].y += coord[1]
    return before


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


def test_quick_recleans_mixed_structure_even_if_clean_hash_was_seen_before() -> None:
    graph, phenyl, cyclopropane = _mixed_chain_phenyl_cyclopropane()
    first = run_clean2d_engine(graph, mode="quick", target_bond_length=40.0)
    assert first.ok
    assert first.selected is not None
    clean_coords = dict(first.selected.coords)
    clean_hash = clean2d_geometry_hash(graph, clean_coords)
    _apply_coords(graph, clean_coords)

    clean_snapshot = capture_clean2d_snapshot(graph)
    clean_before = {atom_id: (atom.x, atom.y) for atom_id, atom in graph.atoms.items()}
    _distort_mixed_graph(graph)
    distorted = {atom_id: (atom.x, atom.y) for atom_id, atom in graph.atoms.items()}

    second = run_clean2d_engine(
        graph,
        mode="quick",
        target_bond_length=40.0,
        avoid_hashes={clean_hash},
    )

    assert second.ok
    assert second.selected is not None
    assert second.selected.source != "current"
    assert not any(item.rejection_reason == "geometria_repetida" for item in second.rejected)

    lengths = []
    for bond in graph.bonds.values():
        a = second.selected.coords[bond.a1_id]
        b = second.selected.coords[bond.a2_id]
        lengths.append(math.hypot(b[0] - a[0], b[1] - a[1]))
    assert min(lengths) >= 40.0 * 0.75
    assert max(lengths) <= 40.0 * 1.30
    assert count_new_bond_crossings(distorted, second.selected.coords, list(graph.bonds.values())) == 0
    assert ring_degeneracy_score(second.selected.coords, phenyl) > 0.18
    assert ring_degeneracy_score(second.selected.coords, cyclopropane) > 0.05
    assert_clean2d_invariants(clean_snapshot, graph, clean_before, second.selected.coords)


def test_propose_mode_preserves_wedge_hash_metadata_without_bad_geometry() -> None:
    graph = MolGraph()
    center = graph.add_atom("C", 0.0, 0.0)
    left = graph.add_atom("C", -40.0, 0.0)
    right = graph.add_atom("C", 40.0, 0.0)
    up = graph.add_atom("Cl", 0.0, -40.0)
    down = graph.add_atom("Br", 0.0, 40.0)
    graph.add_bond(center.id, left.id, order=1)
    graph.add_bond(center.id, right.id, order=1)
    wedge = graph.add_bond(center.id, up.id, order=1, style=BondStyle.WEDGE)
    hashed = graph.add_bond(center.id, down.id, order=1, style=BondStyle.HASHED)
    before = {atom_id: (atom.x, atom.y) for atom_id, atom in graph.atoms.items()}
    snapshot = capture_clean2d_snapshot(graph)

    result = run_clean2d_engine(graph, mode="propose", target_bond_length=40.0)

    assert result.ok
    assert result.selected is not None
    coords = result.selected.coords
    assert all(math.isfinite(value) for xy in coords.values() for value in xy)
    assert graph.get_bond(wedge.id).style == BondStyle.WEDGE
    assert graph.get_bond(hashed.id).style == BondStyle.HASHED
    lengths = [
        math.hypot(coords[bond.a2_id][0] - coords[bond.a1_id][0], coords[bond.a2_id][1] - coords[bond.a1_id][1])
        for bond in graph.bonds.values()
    ]
    assert min(lengths) >= 40.0 * 0.7
    assert max(lengths) <= 40.0 * 1.4
    assert count_new_bond_crossings(before, coords, list(graph.bonds.values())) == 0
    assert_clean2d_invariants(snapshot, graph, before, coords)
