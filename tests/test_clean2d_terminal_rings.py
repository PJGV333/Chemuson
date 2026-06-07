from __future__ import annotations

import math
import os
import sys
from unittest.mock import patch

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.clean2d import (
    Clean2DCandidate,
    assert_clean2d_invariants,
    capture_clean2d_snapshot,
    classify_clean2d_layout_quality,
    clean2d_geometry_hash,
    count_new_bond_crossings,
    rank_clean2d_candidates,
    ring_degeneracy_score,
    run_clean2d_engine,
)
from chemuson.core.model import MolGraph


def _coords(graph: MolGraph) -> dict[int, tuple[float, float]]:
    return {atom_id: (atom.x, atom.y) for atom_id, atom in graph.atoms.items()}


def _angle_between(
    center: tuple[float, float],
    left: tuple[float, float],
    right: tuple[float, float],
) -> float:
    a1 = math.atan2(left[1] - center[1], left[0] - center[0])
    a2 = math.atan2(right[1] - center[1], right[0] - center[0])
    return abs((math.degrees(a2 - a1) + 180.0) % 360.0 - 180.0)


def _ring_center(coords: dict[int, tuple[float, float]], ring: set[int]) -> tuple[float, float]:
    return (
        sum(coords[atom_id][0] for atom_id in ring) / len(ring),
        sum(coords[atom_id][1] for atom_id in ring) / len(ring),
    )


def _adjacency(graph: MolGraph) -> dict[int, set[int]]:
    adjacency: dict[int, set[int]] = {atom_id: set() for atom_id in graph.atoms}
    for bond in graph.bonds.values():
        adjacency.setdefault(bond.a1_id, set()).add(bond.a2_id)
        adjacency.setdefault(bond.a2_id, set()).add(bond.a1_id)
    return adjacency


def _terminal_ring_attachment(graph: MolGraph, ring: set[int]) -> tuple[int, int]:
    adjacency = _adjacency(graph)
    attachments: list[tuple[int, int]] = []
    for atom_id in sorted(ring):
        external = sorted(neigh for neigh in adjacency.get(atom_id, set()) if neigh not in ring)
        attachments.extend((atom_id, neigh) for neigh in external)
    assert len(attachments) == 1
    return attachments[0]


def _ring_bond_lengths(
    graph: MolGraph,
    coords: dict[int, tuple[float, float]],
    ring: set[int],
) -> list[float]:
    adjacency = _adjacency(graph)
    return [
        math.dist(coords[atom_id], coords[neigh])
        for atom_id in sorted(ring)
        for neigh in sorted(adjacency.get(atom_id, set()) & ring)
        if atom_id < neigh
    ]


def _assert_terminal_ring_points_outward(
    coords: dict[int, tuple[float, float]],
    *,
    ring: set[int],
    anchor: int,
    external: int,
    min_center_angle: float = 120.0,
    min_edge_angle: float = 90.0,
) -> None:
    center = _ring_center(coords, ring)
    center_angle = _angle_between(coords[anchor], coords[external], center)
    assert center_angle > min_center_angle
    ring_neighbors = sorted(
        atom_id
        for atom_id in ring
        if atom_id != anchor
        and math.isclose(math.dist(coords[anchor], coords[atom_id]), 40.0, rel_tol=0.12)
    )
    assert len(ring_neighbors) >= 2
    edge_angles = [
        _angle_between(coords[anchor], coords[external], coords[neigh])
        for neigh in ring_neighbors[:2]
    ]
    assert min(edge_angles) > min_edge_angle


def _raw_mixed_with_terminal_cyclopropane() -> tuple[MolGraph, set[int], set[int]]:
    graph = MolGraph()
    phenyl: list[int] = []
    radius = 40.0
    for idx in range(6):
        angle = math.radians(60.0 * idx + 30.0)
        phenyl.append(graph.add_atom("C", radius * math.cos(angle), radius * math.sin(angle)).id)
    for idx in range(6):
        graph.add_bond(phenyl[idx], phenyl[(idx + 1) % 6], order=1, is_aromatic=True)

    chain = [
        graph.add_atom("C", 90.0, 20.0).id,
        graph.add_atom("C", 130.0, 20.0).id,
        graph.add_atom("C", 150.0, 54.6).id,
        graph.add_atom("C", 190.0, 54.6).id,
    ]
    graph.add_bond(phenyl[0], chain[0], order=1)
    graph.add_bond(chain[0], chain[1], order=1)
    graph.add_bond(chain[1], chain[2], order=2)
    graph.add_bond(chain[2], chain[3], order=1)
    branch = graph.add_atom("C", 130.0, -20.0)
    graph.add_bond(chain[1], branch.id, order=1)

    cyclopropane = {
        graph.add_atom("C", 230.0, 54.6).id,
        graph.add_atom("C", 270.0, 54.6).id,
        graph.add_atom("C", 250.0, 89.2).id,
    }
    cyc = sorted(cyclopropane)
    graph.add_bond(chain[3], cyc[0], order=1)
    graph.add_bond(cyc[0], cyc[1], order=1)
    graph.add_bond(cyc[1], cyc[2], order=1)
    graph.add_bond(cyc[2], cyc[0], order=1)
    return graph, set(phenyl), set(cyc)


def test_terminal_cyclopropane_canonicalizes_outward() -> None:
    graph, _phenyl, cyclopropane = _raw_mixed_with_terminal_cyclopropane()
    before = _coords(graph)
    snapshot = capture_clean2d_snapshot(graph)

    result = run_clean2d_engine(graph, mode="quick", target_bond_length=40.0)

    assert result.ok
    assert result.selected is not None
    coords = result.selected.coords
    quality = classify_clean2d_layout_quality(graph, coords=coords, target_bond_length=40.0)
    assert quality.quality_class == "good"
    assert quality.bad_exocyclic_ring_count == 0
    anchor, external = _terminal_ring_attachment(graph, cyclopropane)
    _assert_terminal_ring_points_outward(coords, ring=cyclopropane, anchor=anchor, external=external)
    ring_lengths = _ring_bond_lengths(graph, coords, cyclopropane)
    assert len(ring_lengths) == 3
    assert min(ring_lengths) >= 40.0 * 0.88
    assert max(ring_lengths) <= 40.0 * 1.12
    assert count_new_bond_crossings(before, coords, list(graph.bonds.values())) == 0
    assert ring_degeneracy_score(coords, cyclopropane) > 0.18
    assert_clean2d_invariants(snapshot, graph, before, coords)


def _fused_bicyclic_graph() -> MolGraph:
    graph = MolGraph()
    coords = {
        1: (0.0, 0.0),
        2: (40.0, 0.0),
        3: (60.0, 34.6),
        4: (40.0, 69.2),
        5: (0.0, 69.2),
        6: (-20.0, 34.6),
        7: (100.0, 0.0),
        8: (120.0, 34.6),
        9: (100.0, 69.2),
        10: (60.0, 69.2),
    }
    for atom_id, (x, y) in coords.items():
        graph.add_atom("C", x, y, atom_id=atom_id)
    for a_id, b_id in (
        (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1),
        (3, 7), (7, 8), (8, 9), (9, 10), (10, 4),
    ):
        graph.add_bond(a_id, b_id, order=1)
    return graph


def test_single_anchor_logic_does_not_run_for_fused_rings() -> None:
    graph = _fused_bicyclic_graph()

    with patch("chemuson.clean2d.engine._place_ring_from_single_anchor", side_effect=AssertionError):
        result = run_clean2d_engine(graph, mode="quick", target_bond_length=40.0)

    assert result.ok
    assert result.selected is not None


def _terminal_cyclohexyl_graph() -> tuple[MolGraph, set[int], int, int]:
    graph = MolGraph()
    c1 = graph.add_atom("C", 0.0, 0.0, atom_id=1)
    c2 = graph.add_atom("C", 40.0, 0.0, atom_id=2)
    graph.add_bond(c1.id, c2.id, order=1)
    ring: list[int] = []
    for idx in range(6):
        angle = math.radians(60.0 * idx)
        ring.append(graph.add_atom("C", 95.0 + 35.0 * math.cos(angle), 8.0 + 35.0 * math.sin(angle)).id)
    graph.add_bond(c2.id, ring[0], order=1)
    for idx in range(6):
        graph.add_bond(ring[idx], ring[(idx + 1) % 6], order=1)
    return graph, set(ring), ring[0], c2.id


def test_terminal_six_member_ring_points_away_from_chain() -> None:
    graph, ring, anchor, external = _terminal_cyclohexyl_graph()

    result = run_clean2d_engine(graph, mode="quick", target_bond_length=40.0)

    assert result.ok
    assert result.selected is not None
    coords = result.selected.coords
    quality = classify_clean2d_layout_quality(graph, coords=coords, target_bond_length=40.0)
    assert quality.quality_class == "good"
    _assert_terminal_ring_points_outward(
        coords,
        ring=ring,
        anchor=anchor,
        external=external,
        min_center_angle=120.0,
        min_edge_angle=90.0,
    )


def _cyclopropane_orientation_candidate(
    graph: MolGraph,
    ring: set[int],
    *,
    outward: bool,
    base_coords: dict[int, tuple[float, float]] | None = None,
) -> dict[int, tuple[float, float]]:
    coords = dict(base_coords or _coords(graph))
    anchor_id, external_id = _terminal_ring_attachment(graph, ring)
    anchor = coords[anchor_id]
    external = coords[external_id]
    vx = external[0] - anchor[0]
    vy = external[1] - anchor[1]
    norm = math.hypot(vx, vy)
    ux, uy = vx / norm, vy / norm
    radius = 40.0 / (2.0 * math.sin(math.pi / 3.0))
    direction = -1.0 if outward else 1.0
    center = (anchor[0] + direction * ux * radius, anchor[1] + direction * uy * radius)
    anchor_angle = math.atan2(anchor[1] - center[1], anchor[0] - center[0])
    out = dict(coords)
    ordered_ring = [anchor_id] + [atom_id for atom_id in sorted(ring) if atom_id != anchor_id]
    for idx, atom_id in enumerate(ordered_ring):
        angle = anchor_angle + 2.0 * math.pi * idx / 3.0
        out[atom_id] = (center[0] + radius * math.cos(angle), center[1] + radius * math.sin(angle))
    out[anchor_id] = anchor
    return out


def test_ranking_prefers_terminal_ring_outward_orientation() -> None:
    graph, _phenyl, cyclopropane = _raw_mixed_with_terminal_cyclopropane()
    before = _coords(graph)
    canonical = run_clean2d_engine(graph, mode="quick", target_bond_length=40.0)
    assert canonical.selected is not None
    base_coords = canonical.selected.coords
    outward = _cyclopropane_orientation_candidate(graph, cyclopropane, outward=True, base_coords=base_coords)
    inward = _cyclopropane_orientation_candidate(graph, cyclopropane, outward=False, base_coords=base_coords)
    inward_hash = clean2d_geometry_hash(graph, inward)
    anchor, external = _terminal_ring_attachment(graph, cyclopropane)

    result = rank_clean2d_candidates(
        graph,
        [
            Clean2DCandidate(source="rdkit_isolated", coords=inward, message="inward"),
            Clean2DCandidate(source="safe_fallback", coords=outward, message="outward"),
        ],
        before,
        set(graph.atoms),
        mode="quick",
        target_bond_length=40.0,
    )

    assert result.selected is not None
    assert result.selected.geometry_hash != inward_hash
    assert result.selected.message == "Estructura 2D optimizada"
    _assert_terminal_ring_points_outward(result.selected.coords, ring=cyclopropane, anchor=anchor, external=external)
