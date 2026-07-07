from __future__ import annotations

import json
import math

import pytest

from chemuson.clean2d import (
    CLEAN2D_DEBUG_SNAPSHOT_ENV,
    CLEAN2D_DEBUG_SNAPSHOT_SCHEMA,
    CLEAN2D_DEBUG_SNAPSHOT_VERSION,
    read_clean2d_debug_snapshot,
    run_clean2d_engine,
    run_clean2d_with_debug_snapshot,
    validate_clean2d_debug_snapshot,
)
from chemuson.core.model import BondStyle, MolGraph


def _two_atom_graph(distance: float = 40.0) -> MolGraph:
    graph = MolGraph()
    graph.add_atom("C", 0.0, 0.0, atom_id=1)
    graph.add_atom("O", distance, 0.0, atom_id=2)
    graph.add_bond(1, 2, order=2, bond_id=1)
    return graph


def _add_hexagon(graph: MolGraph, cx: float, cy: float, radius: float = 24.0) -> list[int]:
    atoms = []
    for idx in range(6):
        angle = math.radians(idx * 60.0 + 30.0)
        atoms.append(graph.add_atom("C", cx + math.cos(angle) * radius, cy + math.sin(angle) * radius).id)
    for idx in range(6):
        graph.add_bond(atoms[idx], atoms[(idx + 1) % 6], order=1, is_aromatic=True)
    return atoms


def _naphthalene_like_graph() -> MolGraph:
    graph = MolGraph()
    for idx in range(10):
        graph.add_atom("C", float((idx % 5) * 5), float((idx // 5) * 3), atom_id=idx + 1)
    for a1, a2 in [
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 6),
        (6, 1),
        (4, 7),
        (7, 8),
        (8, 9),
        (9, 10),
        (10, 5),
    ]:
        graph.add_bond(a1, a2, order=1, is_aromatic=True)
    return graph


def _tetrandrine_like_graph() -> MolGraph:
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
    for left, right in [
        (rings[0][1], rings[1][5]),
        (rings[1][1], rings[2][5]),
        (rings[2][1], rings[3][5]),
        (rings[3][1], rings[0][5]),
    ]:
        bridge = graph.add_atom("C", 32.0, 26.0).id
        graph.add_bond(left, bridge, order=1)
        graph.add_bond(bridge, right, order=1)
    return graph


def _distorted_mixed_graph() -> MolGraph:
    graph = MolGraph()
    phenyl = _add_hexagon(graph, 0.0, 0.0, radius=40.0)
    chain = [graph.add_atom("C", 180.0 + idx * 90.0, float((idx % 2) * 8)).id for idx in range(4)]
    graph.add_bond(phenyl[0], chain[0], order=1)
    graph.add_bond(chain[0], chain[1], order=1, style=BondStyle.BOLD)
    graph.add_bond(chain[1], chain[2], order=2, stereo_ez="E")
    graph.add_bond(chain[2], chain[3], order=1)
    cyclopropane = [graph.add_atom("C", 620.0 + idx * 70.0, float(idx * 2)).id for idx in range(3)]
    graph.add_bond(chain[3], cyclopropane[0], order=1)
    graph.add_bond(cyclopropane[0], cyclopropane[1], order=1)
    graph.add_bond(cyclopropane[1], cyclopropane[2], order=1)
    graph.add_bond(cyclopropane[2], cyclopropane[0], order=1)
    return graph


def test_debug_snapshots_are_disabled_by_default_and_do_not_change_result(monkeypatch) -> None:
    monkeypatch.delenv(CLEAN2D_DEBUG_SNAPSHOT_ENV, raising=False)
    graph = _two_atom_graph(distance=75.0)

    baseline = run_clean2d_engine(graph, mode="quick", target_bond_length=40.0)
    debug = run_clean2d_engine(graph, mode="quick", target_bond_length=40.0, debug_snapshot=True)

    assert baseline.debug_snapshot is None
    assert debug.debug_snapshot is not None
    assert baseline.result_state == debug.result_state
    assert baseline.stable_reason == debug.stable_reason
    assert baseline.candidate_sources == debug.candidate_sources
    assert (baseline.selected.source if baseline.selected else None) == (debug.selected.source if debug.selected else None)
    assert (baseline.selected.coords if baseline.selected else None) == (debug.selected.coords if debug.selected else None)


def test_environment_variable_enables_debug_snapshot(monkeypatch) -> None:
    monkeypatch.setenv(CLEAN2D_DEBUG_SNAPSHOT_ENV, "1")
    result = run_clean2d_engine(_two_atom_graph(), mode="quick", target_bond_length=40.0)

    assert result.debug_snapshot is not None
    assert result.debug_snapshot["schema"] == CLEAN2D_DEBUG_SNAPSHOT_SCHEMA
    assert result.debug_snapshot["version"] == CLEAN2D_DEBUG_SNAPSHOT_VERSION


def test_explicit_debug_parameter_and_test_helper_write_read_snapshot(tmp_path, monkeypatch) -> None:
    monkeypatch.delenv(CLEAN2D_DEBUG_SNAPSHOT_ENV, raising=False)
    path = tmp_path / "snapshot.json"
    result = run_clean2d_with_debug_snapshot(
        _two_atom_graph(),
        path=path,
        mode="quick",
        target_bond_length=40.0,
        initial_selection={"atom_ids": [2], "bond_ids": [1]},
        metadata={"case": "unit"},
    )

    assert result.debug_snapshot is not None
    assert path.exists()
    loaded = read_clean2d_debug_snapshot(path)
    assert loaded == result.debug_snapshot
    assert loaded["initial_selection"] == {"atom_ids": [2], "bond_ids": [1]}
    assert loaded["metadata"] == {"case": "unit"}


def test_snapshot_records_topology_coordinates_target_candidates_and_decision(monkeypatch) -> None:
    monkeypatch.delenv(CLEAN2D_DEBUG_SNAPSHOT_ENV, raising=False)
    result = run_clean2d_engine(
        _two_atom_graph(distance=40.0),
        [1, 2],
        mode="quick",
        target_bond_length=40.0,
        debug_snapshot=True,
        debug_initial_selection={"atom_ids": [1], "bond_ids": [1]},
    )

    snapshot = result.debug_snapshot
    assert snapshot is not None
    validate_clean2d_debug_snapshot(snapshot)
    json.dumps(snapshot, allow_nan=False)
    assert snapshot["mode"] == "quick"
    assert snapshot["target"] == {"whole_structure": False, "atom_ids": [1, 2]}
    assert snapshot["topology"]["atoms"] == [{"id": 1, "element": "C"}, {"id": 2, "element": "O"}]
    assert snapshot["topology"]["bonds"] == [{"id": 1, "atom_ids": [1, 2], "order": 2, "is_aromatic": False}]
    assert snapshot["initial_coordinates"] == [
        {"atom_id": 1, "x": 0.0, "y": 0.0},
        {"atom_id": 2, "x": 40.0, "y": 0.0},
    ]
    assert snapshot["final_coordinates"] is not None
    assert [row["source"] for row in snapshot["candidates"]] == list(result.candidate_sources)
    assert all("quality_diagnostic" in row for row in snapshot["candidates"])
    assert snapshot["final_state"] == result.result_state
    assert snapshot["final_reason"] == (result.stable_reason or None)


def test_snapshot_represents_missing_final_coordinates_safely(monkeypatch) -> None:
    monkeypatch.delenv(CLEAN2D_DEBUG_SNAPSHOT_ENV, raising=False)
    result = run_clean2d_engine(
        _two_atom_graph(distance=200.0),
        mode="propose",
        target_bond_length=40.0,
        debug_snapshot=True,
    )

    snapshot = result.debug_snapshot
    assert snapshot is not None
    assert result.selected is None
    assert snapshot["final_coordinates"] is None
    assert snapshot["final_state"] == "failed-controlled"
    assert snapshot["final_reason"]
    json.dumps(snapshot, allow_nan=False)


@pytest.mark.parametrize(
    ("name", "builder"),
    [
        ("naphthalene_fused_system", _naphthalene_like_graph),
        ("tetrandrine_like_hierarchical_blocks", _tetrandrine_like_graph),
        ("smart_clean_distorted_mixed_structure", _distorted_mixed_graph),
        ("complex_engine_global_redraw_guard_shape", _tetrandrine_like_graph),
    ],
)
def test_known_complex_shapes_can_emit_structural_snapshots_without_algorithm_assertions(
    name: str,
    builder,
    monkeypatch,
) -> None:
    monkeypatch.delenv(CLEAN2D_DEBUG_SNAPSHOT_ENV, raising=False)
    graph = builder()

    result = run_clean2d_engine(
        graph,
        mode="quick",
        target_bond_length=40.0,
        debug_snapshot=True,
        debug_snapshot_metadata={"case": name},
    )

    snapshot = result.debug_snapshot
    assert snapshot is not None
    assert snapshot["metadata"] == {"case": name}
    assert snapshot["topology"]["atoms"]
    assert snapshot["initial_coordinates"]
    assert snapshot["candidates"] or snapshot["final_state"] == "failed-controlled"
    validate_clean2d_debug_snapshot(snapshot)
