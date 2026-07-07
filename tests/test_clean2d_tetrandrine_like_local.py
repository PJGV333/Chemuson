from __future__ import annotations

import math
from unittest.mock import patch

import pytest


from chemuson.clean2d import (
    classify_clean2d_complexity,
    local_graph_clean2d,
    ring_degeneracy_score,
    run_clean2d_engine,
    stereo_layout_signature,
)
from chemuson.core.model import BondStyle, MolGraph


TARGET = 40.0


def _add_benzene(graph: MolGraph, center: tuple[float, float]) -> list[int]:
    ring: list[int] = []
    for idx in range(6):
        angle = math.radians(60.0 * idx + 30.0)
        ring.append(graph.add_atom("C", center[0] + TARGET * math.cos(angle), center[1] + TARGET * math.sin(angle)).id)
    for idx in range(6):
        graph.add_bond(ring[idx], ring[(idx + 1) % 6], order=1, is_aromatic=True)
    return ring


def _unit_from(center: tuple[float, float], point: tuple[float, float]) -> tuple[float, float]:
    dx = point[0] - center[0]
    dy = point[1] - center[1]
    length = math.hypot(dx, dy)
    return dx / length, dy / length


def _tetrandrine_like_graph() -> tuple[MolGraph, dict[str, object]]:
    graph = MolGraph()
    centers = {
        "r1": (-120.0, 0.0),
        "r2": (120.0, 0.0),
        "r3": (-120.0, 140.0),
        "r4": (120.0, 140.0),
    }
    rings = {name: _add_benzene(graph, center) for name, center in centers.items()}

    p25 = graph.add_atom("C", -43.0, 20.0).id
    p26 = graph.add_atom("C", 0.0, 0.0, stereo_cip="R").id
    p27 = graph.add_atom("C", 43.0, -20.0).id
    graph.add_bond(rings["r1"][0], p25, order=1)
    graph.add_bond(p25, p26, order=1)
    graph.add_bond(p26, p27, order=1)
    graph.add_bond(p27, rings["r2"][3], order=1)

    p28 = graph.add_atom("O", -154.6, 60.0).id
    p29 = graph.add_atom("C", -154.6, 100.0).id
    graph.add_bond(rings["r1"][2], p28, order=1)
    graph.add_bond(p28, p29, order=1)
    graph.add_bond(p29, rings["r3"][4], order=1)

    p30 = graph.add_atom("O", 154.6, 60.0).id
    p31 = graph.add_atom("C", 154.6, 100.0).id
    graph.add_bond(rings["r2"][1], p30, order=1)
    graph.add_bond(p30, p31, order=1)
    graph.add_bond(p31, rings["r4"][5], order=1)

    p32 = graph.add_atom("C", -40.0, 160.0).id
    p33 = graph.add_atom("C", 0.0, 160.0).id
    p34 = graph.add_atom("C", 40.0, 160.0).id
    graph.add_bond(rings["r3"][0], p32, order=1)
    graph.add_bond(p32, p33, order=1)
    graph.add_bond(p33, p34, order=1)
    graph.add_bond(p34, rings["r4"][3], order=1)

    cl = graph.add_atom("Cl", 0.0, -38.0).id
    br = graph.add_atom("Br", 0.0, 38.0).id
    graph.add_bond(p26, cl, order=1, style=BondStyle.WEDGE)
    graph.add_bond(p26, br, order=1, style=BondStyle.HASHED)

    anchor = rings["r4"][0]
    ux, uy = _unit_from(centers["r4"], (graph.atoms[anchor].x, graph.atoms[anchor].y))
    tail1 = graph.add_atom("C", graph.atoms[anchor].x + ux * TARGET, graph.atoms[anchor].y + uy * TARGET).id
    tail2 = graph.add_atom("O", graph.atoms[tail1].x + ux * TARGET, graph.atoms[tail1].y + uy * TARGET).id
    graph.add_bond(anchor, tail1, order=1)
    graph.add_bond(tail1, tail2, order=1)

    return graph, {"rings": rings, "centers": centers, "tail_bond": (anchor, tail1), "tail_atoms": (tail1, tail2)}


def _clean_complex_graph() -> MolGraph:
    graph = MolGraph()
    for center in [(-220.0, 0.0), (-110.0, 0.0), (0.0, 0.0), (110.0, 0.0), (220.0, 0.0)]:
        _add_benzene(graph, center)
    return graph


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


def _bbox_diag(coords: dict[int, tuple[float, float]]) -> float:
    xs = [x for x, _y in coords.values()]
    ys = [y for _x, y in coords.values()]
    return math.hypot(max(xs) - min(xs), max(ys) - min(ys))


def test_complex_local_deformed_bond_repairs_without_global_redraw() -> None:
    graph, info = _tetrandrine_like_graph()
    before = _coords(graph)
    signature_before = stereo_layout_signature(graph, before)
    rings = info["rings"]
    assert isinstance(rings, dict)
    ring_scores_before = {
        name: ring_degeneracy_score(before, set(ring))
        for name, ring in rings.items()
    }
    anchor, tail1 = info["tail_bond"]
    tail1_id, tail2_id = info["tail_atoms"]
    assert isinstance(anchor, int) and isinstance(tail1, int)
    assert isinstance(tail1_id, int) and isinstance(tail2_id, int)

    ux = (graph.atoms[tail1_id].x - graph.atoms[anchor].x) / TARGET
    uy = (graph.atoms[tail1_id].y - graph.atoms[anchor].y) / TARGET
    for atom_id in (tail1_id, tail2_id):
        graph.atoms[atom_id].x += ux * 35.0
        graph.atoms[atom_id].y += uy * 35.0
    distorted = _coords(graph)
    before_len = _distance(distorted, anchor, tail1)
    assert before_len > TARGET * 1.8

    result = local_graph_clean2d(graph, target_bond_length=TARGET, mode="quick")

    assert result.ok
    assert result.report.number_of_repairs_accepted > 0
    after_len = _distance(result.coords, anchor, tail1)
    assert after_len < before_len
    assert result.report.bond_integrity_regressions == 0
    assert result.report.mean_displacement < TARGET * 0.12
    assert result.report.max_displacement <= TARGET * 0.75
    assert stereo_layout_signature(graph, result.coords) == signature_before
    assert _bbox_diag(result.coords) == pytest.approx(_bbox_diag(distorted), rel=0.08)
    for name, ring in rings.items():
        assert ring_degeneracy_score(result.coords, set(ring)) >= ring_scores_before[name] - 0.03


def test_bond_integrity_after_local_repair() -> None:
    graph, info = _tetrandrine_like_graph()
    anchor, tail1 = info["tail_bond"]
    tail1_id, tail2_id = info["tail_atoms"]
    assert isinstance(anchor, int) and isinstance(tail1, int)
    assert isinstance(tail1_id, int) and isinstance(tail2_id, int)
    before_clean = _coords(graph)
    ux = (graph.atoms[tail1_id].x - graph.atoms[anchor].x) / TARGET
    uy = (graph.atoms[tail1_id].y - graph.atoms[anchor].y) / TARGET
    for atom_id in (tail1_id, tail2_id):
        graph.atoms[atom_id].x += ux * 35.0
        graph.atoms[atom_id].y += uy * 35.0
    distorted = _coords(graph)
    before_lengths = _bond_lengths(graph, distorted)
    target_bond = next(
        bond.id
        for bond in graph.bonds.values()
        if {bond.a1_id, bond.a2_id} == {anchor, tail1}
    )

    result = local_graph_clean2d(graph, target_bond_length=TARGET, mode="quick")

    assert result.ok
    assert result.report.bond_integrity_regressions == 0
    after_lengths = _bond_lengths(graph, result.coords)
    for bond in graph.bonds.values():
        if bond.id == target_bond:
            continue
        before_len = before_lengths[bond.id]
        after_len = after_lengths[bond.id]
        clean_len = _distance(before_clean, bond.a1_id, bond.a2_id)
        if 0.70 * TARGET <= before_len <= 1.35 * TARGET:
            assert 0.70 * TARGET <= after_len <= 1.35 * TARGET
        assert after_len <= max(before_len * 1.20, clean_len * 1.20)


def test_tetrandrine_like_no_visual_disconnection() -> None:
    graph, info = _tetrandrine_like_graph()
    before = _coords(graph)
    signature_before = stereo_layout_signature(graph, before)
    rings = info["rings"]
    assert isinstance(rings, dict)

    result = local_graph_clean2d(graph, target_bond_length=TARGET, mode="quick")

    assert result.report.bond_integrity_regressions == 0
    assert stereo_layout_signature(graph, result.coords) == signature_before
    for bond in graph.bonds.values():
        before_len = _distance(before, bond.a1_id, bond.a2_id)
        after_len = _distance(result.coords, bond.a1_id, bond.a2_id)
        if before_len <= TARGET * 1.35:
            assert after_len <= TARGET * 1.35
        assert after_len <= max(before_len * 1.20, before_len + 1.0)
    for ring in rings.values():
        assert isinstance(ring, list)
        assert ring_degeneracy_score(result.coords, set(ring)) >= ring_degeneracy_score(before, set(ring)) - 0.03


def test_complex_reasonable_structure_reports_no_safe_local_defects() -> None:
    graph = _clean_complex_graph()

    result = local_graph_clean2d(graph, target_bond_length=TARGET, mode="quick")

    assert not result.ok
    assert result.changed_coords == {}
    assert result.report.number_of_repairs_accepted == 0
    assert result.report.reason in {"no_defects", "no_safe_local_repair"}
    assert "no se detectaron defectos locales seguros" in result.report.message


def test_complex_engine_does_not_call_global_redraw_candidates() -> None:
    graph = _clean_complex_graph()
    profile = classify_clean2d_complexity(graph)

    assert profile.global_redraw_allowed is False
    assert profile.local_repair_allowed is False

    with patch("chemuson.clean2d.engine._candidate_from_rdkit_isolated") as rdkit_candidate:
        rdkit_candidate.side_effect = AssertionError("global redraw must not run for complex quick clean")
        result = run_clean2d_engine(graph, mode="quick", target_bond_length=TARGET)

    rdkit_candidate.assert_not_called()
    assert result.selected is not None
    assert result.selected.source == "current"
    assert result.selected.metadata["complex_policy"]["global_redraw_allowed"] is False
    assert "no se redibujó" in result.message
