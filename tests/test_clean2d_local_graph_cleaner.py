from __future__ import annotations

import math
import os
import sys

import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.clean2d import (
    LocalDefect,
    LocalRepair,
    local_graph_clean2d,
    ring_degeneracy_score,
    run_clean2d_engine,
    validate_local_repair,
)
from chemuson.core.model import MolGraph


TARGET = 40.0


def _add_benzene(graph: MolGraph, center: tuple[float, float]) -> list[int]:
    ring: list[int] = []
    for idx in range(6):
        angle = math.radians(60.0 * idx + 30.0)
        ring.append(graph.add_atom("C", center[0] + TARGET * math.cos(angle), center[1] + TARGET * math.sin(angle)).id)
    for idx in range(6):
        graph.add_bond(ring[idx], ring[(idx + 1) % 6], order=1, is_aromatic=True)
    return ring


def _complex_three_ring_graph() -> tuple[MolGraph, dict[str, object]]:
    graph = MolGraph()
    r1 = _add_benzene(graph, (-100.0, 0.0))
    r2 = _add_benzene(graph, (100.0, 0.0))
    r3 = _add_benzene(graph, (0.0, 150.0))
    c19 = graph.add_atom("C", -42.0, 20.0).id
    c20 = graph.add_atom("O", 0.0, 0.0).id
    c21 = graph.add_atom("C", 42.0, -20.0).id
    graph.add_bond(r1[0], c19, order=1)
    graph.add_bond(c19, c20, order=1)
    graph.add_bond(c20, c21, order=1)
    graph.add_bond(c21, r2[3], order=1)
    c22 = graph.add_atom("C", -80.0, 70.0).id
    c23 = graph.add_atom("O", -45.0, 105.0).id
    graph.add_bond(r1[1], c22, order=1)
    graph.add_bond(c22, c23, order=1)
    graph.add_bond(c23, r3[4], order=1)
    c24 = graph.add_atom("C", 80.0, 70.0).id
    c25 = graph.add_atom("O", 45.0, 105.0).id
    graph.add_bond(r2[2], c24, order=1)
    graph.add_bond(c24, c25, order=1)
    graph.add_bond(c25, r3[5], order=1)
    c26 = graph.add_atom("C", 0.0, 105.0).id
    c27 = graph.add_atom("C", 0.0, 65.0).id
    graph.add_bond(r3[3], c26, order=1)
    graph.add_bond(c26, c27, order=1)
    graph.add_bond(c27, c20, order=1)
    return graph, {"rings": {"r1": r1, "r2": r2, "r3": r3}}


def _coords(graph: MolGraph) -> dict[int, tuple[float, float]]:
    return {atom_id: (atom.x, atom.y) for atom_id, atom in graph.atoms.items()}


def _distance(coords: dict[int, tuple[float, float]], a_id: int, b_id: int) -> float:
    return math.hypot(coords[b_id][0] - coords[a_id][0], coords[b_id][1] - coords[a_id][1])


def _bond_lengths(graph: MolGraph, coords: dict[int, tuple[float, float]]) -> dict[int, float]:
    return {
        bond.id: _distance(coords, bond.a1_id, bond.a2_id)
        for bond in graph.bonds.values()
        if bond.a1_id in coords and bond.a2_id in coords
    }


def _center(coords: dict[int, tuple[float, float]], atom_ids: list[int]) -> tuple[float, float]:
    return (
        sum(coords[atom_id][0] for atom_id in atom_ids) / len(atom_ids),
        sum(coords[atom_id][1] for atom_id in atom_ids) / len(atom_ids),
    )


def test_collision_between_flexible_substituent_and_ring_atom_is_separated_locally() -> None:
    graph, info = _complex_three_ring_graph()
    rings = info["rings"]
    assert isinstance(rings, dict)
    r1 = rings["r1"]
    assert isinstance(r1, list)
    anchor = r1[0]
    near_atom = r1[1]
    flexible = graph.add_atom("C", graph.atoms[near_atom].x + 1.0, graph.atoms[near_atom].y + 1.0).id
    graph.add_bond(anchor, flexible, order=1)
    before = _coords(graph)
    before_distance = _distance(before, flexible, near_atom)

    result = local_graph_clean2d(graph, target_bond_length=TARGET, mode="quick")

    assert result.ok
    assert result.report.defects_by_type.get("collision_atom", 0) > 0
    assert _distance(result.coords, flexible, near_atom) > before_distance + TARGET * 0.15
    assert result.report.mean_displacement < TARGET * 0.15
    moved_far = [
        atom_id
        for atom_id, coord in result.coords.items()
        if math.dist(before[atom_id], coord) > TARGET * 0.30
    ]
    assert len(moved_far) <= 3


def test_no_single_atom_ring_move_in_complex() -> None:
    graph, info = _complex_three_ring_graph()
    rings = info["rings"]
    assert isinstance(rings, dict)
    ring = rings["r1"]
    assert isinstance(ring, list)
    before = _coords(graph)
    after = dict(before)
    moved_atom = ring[0]
    after[moved_atom] = (after[moved_atom][0] + TARGET * 0.25, after[moved_atom][1])
    defect = LocalDefect("aromatic_ring", atom_ids=tuple(ring), score=1.0, involves_ring=True)
    repair = LocalRepair(defect, after, (moved_atom,), "repair_aromatic_ring_local")

    validation = validate_local_repair(
        current=before,
        after=after,
        repair=repair,
        graph=graph,
        atom_ids=set(graph.atoms),
        bonds=list(graph.bonds.values()),
        target=TARGET,
        mode="quick",
    )

    assert not validation.accepted
    assert validation.reason in {"movimiento_atomico_anillo", "movimiento_incoherente_bloque_rigido"}


def test_deformed_aromatic_ring_regularizes_preserving_center_and_substituent_direction() -> None:
    graph, info = _complex_three_ring_graph()
    rings = info["rings"]
    assert isinstance(rings, dict)
    r2 = rings["r2"]
    assert isinstance(r2, list)
    anchor = r2[0]
    before = _coords(graph)
    center_before = _center(before, r2)
    outward = (
        (before[anchor][0] - center_before[0]) / TARGET,
        (before[anchor][1] - center_before[1]) / TARGET,
    )
    substituent = graph.add_atom("O", before[anchor][0] + outward[0] * TARGET, before[anchor][1] + outward[1] * TARGET).id
    graph.add_bond(anchor, substituent, order=1)
    graph.atoms[r2[2]].x = center_before[0] + (graph.atoms[r2[2]].x - center_before[0]) * 0.05
    graph.atoms[r2[2]].y = center_before[1] + (graph.atoms[r2[2]].y - center_before[1]) * 0.05
    distorted = _coords(graph)
    center_distorted = _center(distorted, r2)
    score_before = ring_degeneracy_score(distorted, set(r2))

    result = local_graph_clean2d(graph, target_bond_length=TARGET, mode="publication")

    assert result.ok
    assert result.report.defects_by_type.get("aromatic_ring", 0) > 0
    assert ring_degeneracy_score(result.coords, set(r2)) > score_before + 0.15
    center_after = _center(result.coords, r2)
    assert math.dist(center_distorted, center_after) < TARGET * 0.15
    branch = (
        result.coords[substituent][0] - result.coords[anchor][0],
        result.coords[substituent][1] - result.coords[anchor][1],
    )
    outward_after = (
        result.coords[anchor][0] - center_after[0],
        result.coords[anchor][1] - center_after[1],
    )
    assert branch[0] * outward_after[0] + branch[1] * outward_after[1] > 0.0


def test_terminal_substituent_can_move_without_moving_ring() -> None:
    graph, info = _complex_three_ring_graph()
    rings = info["rings"]
    assert isinstance(rings, dict)
    ring = rings["r2"]
    assert isinstance(ring, list)
    anchor = ring[0]
    before_ring_coords = {atom_id: (graph.atoms[atom_id].x, graph.atoms[atom_id].y) for atom_id in ring}
    methoxy_o = graph.add_atom("O", graph.atoms[anchor].x + TARGET * 2.1, graph.atoms[anchor].y).id
    methoxy_c = graph.add_atom("C", graph.atoms[anchor].x + TARGET * 3.1, graph.atoms[anchor].y).id
    graph.add_bond(anchor, methoxy_o, order=1)
    graph.add_bond(methoxy_o, methoxy_c, order=1)
    before = _coords(graph)
    before_len = _distance(before, anchor, methoxy_o)

    result = local_graph_clean2d(graph, target_bond_length=TARGET, mode="quick")

    assert result.ok
    assert _distance(result.coords, anchor, methoxy_o) < before_len
    assert methoxy_o in result.changed_coords
    assert methoxy_c in result.changed_coords
    for atom_id, coord in before_ring_coords.items():
        assert math.dist(coord, result.coords[atom_id]) < 0.6
    assert result.report.bond_integrity_regressions == 0


def test_simple_molecules_keep_existing_clean2d_flow() -> None:
    chain = MolGraph()
    chain_coords = [(0.0, 0.0), (40.0, 0.0), (60.0, 34.6), (100.0, 34.6), (40.0, -40.0)]
    atoms = [chain.add_atom("C", x, y) for x, y in chain_coords]
    for left, right in ((atoms[0], atoms[1]), (atoms[1], atoms[2]), (atoms[2], atoms[3]), (atoms[1], atoms[4])):
        chain.add_bond(left.id, right.id, order=1)

    chain_result = run_clean2d_engine(chain, mode="quick", target_bond_length=TARGET)
    assert chain_result.ok
    assert chain_result.selected is not None
    assert chain_result.selected.source != "local_graph"

    benzene = MolGraph()
    ring = [benzene.add_atom("C", float(idx * 5), float((idx % 2) * 3)).id for idx in range(6)]
    for idx in range(6):
        benzene.add_bond(ring[idx], ring[(idx + 1) % 6], order=1, is_aromatic=True)
    subst = benzene.add_atom("O", 2.0, 2.0)
    benzene.add_bond(ring[0], subst.id, order=1)

    benzene_result = run_clean2d_engine(benzene, mode="publication", target_bond_length=TARGET)
    assert benzene_result.ok
    assert benzene_result.selected is not None
    assert benzene_result.selected.source != "local_graph"
    assert ring_degeneracy_score(benzene_result.selected.coords, set(ring)) > 0.18

    mixed = MolGraph()
    phenyl = _add_benzene(mixed, (0.0, 0.0))
    c1 = mixed.add_atom("C", 90.0, 20.0).id
    c2 = mixed.add_atom("C", 130.0, 20.0).id
    c3 = mixed.add_atom("C", 150.0, 54.6).id
    mixed.add_bond(phenyl[0], c1, order=1)
    mixed.add_bond(c1, c2, order=1)
    mixed.add_bond(c2, c3, order=1)
    cyclopropane = [
        mixed.add_atom("C", 190.0, 54.6).id,
        mixed.add_atom("C", 230.0, 54.6).id,
        mixed.add_atom("C", 210.0, 89.2).id,
    ]
    mixed.add_bond(c3, cyclopropane[0], order=1)
    mixed.add_bond(cyclopropane[0], cyclopropane[1], order=1)
    mixed.add_bond(cyclopropane[1], cyclopropane[2], order=1)
    mixed.add_bond(cyclopropane[2], cyclopropane[0], order=1)

    mixed_result = run_clean2d_engine(mixed, mode="publication", target_bond_length=TARGET)
    assert mixed_result.ok
    assert mixed_result.selected is not None
    assert mixed_result.selected.source != "local_graph"
