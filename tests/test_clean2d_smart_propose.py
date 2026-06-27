from __future__ import annotations

from dataclasses import asdict
import math
from pathlib import Path
from unittest.mock import MagicMock

import pytest
from PyQt6.QtWidgets import QApplication


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
    summarize_clean2d_candidates,
)
from chemuson.chemio.persistence import PersistenceManager
from chemuson.core.model import BondStyle, MolGraph
from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.controllers.clean2d_controller import Clean2DController


def _ensure_app() -> None:
    QApplication.instance() or QApplication([])


def _populate_clean_mixed_structure(graph: MolGraph) -> tuple[set[int], set[int]]:
    phenyl: list[int] = []
    radius = 40.0
    for idx in range(6):
        angle = math.radians(60.0 * idx + 30.0)
        phenyl.append(
            graph.add_atom(
                "C",
                radius * math.cos(angle),
                radius * math.sin(angle),
                isotope=13 if idx == 1 else None,
            ).id
        )
    for idx in range(6):
        graph.add_bond(phenyl[idx], phenyl[(idx + 1) % 6], order=1, is_aromatic=True)

    chain_coords = [
        (69.2820323027551, 40.0),
        (103.92304845413264, 20.0),
        (138.5640646055102, 40.0),
        (173.20508075688775, 20.0),
    ]
    chain = [
        graph.add_atom("C", x, y, mapping=99 if idx == 1 else None).id
        for idx, (x, y) in enumerate(chain_coords)
    ]
    graph.add_bond(phenyl[0], chain[0], order=1)
    graph.add_bond(chain[0], chain[1], order=1, style=BondStyle.BOLD)
    graph.add_bond(chain[1], chain[2], order=2, stereo_ez="E")
    graph.add_bond(chain[2], chain[3], order=1)
    branch = graph.add_atom("C", 103.92304845413264, -20.0, charge=-1)
    graph.add_bond(chain[1], branch.id, order=1)

    cyclopropane = [
        graph.add_atom("C", 207.84609690826528, 40.0).id,
        graph.add_atom("C", 247.84609690826528, 40.0).id,
        graph.add_atom("C", 227.84609690826528, 74.64101615137754).id,
    ]
    graph.add_bond(chain[3], cyclopropane[0], order=1)
    graph.add_bond(cyclopropane[0], cyclopropane[1], order=1)
    graph.add_bond(cyclopropane[1], cyclopropane[2], order=1)
    graph.add_bond(cyclopropane[2], cyclopropane[0], order=1)
    return set(phenyl), set(cyclopropane)


def _populate_raw_mixed_structure(graph: MolGraph) -> tuple[set[int], set[int]]:
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

    cyclopropane = [
        graph.add_atom("C", 230.0, 54.6).id,
        graph.add_atom("C", 270.0, 54.6).id,
        graph.add_atom("C", 250.0, 89.2).id,
    ]
    graph.add_bond(chain[3], cyclopropane[0], order=1)
    graph.add_bond(cyclopropane[0], cyclopropane[1], order=1)
    graph.add_bond(cyclopropane[1], cyclopropane[2], order=1)
    graph.add_bond(cyclopropane[2], cyclopropane[0], order=1)
    return set(phenyl), set(cyclopropane)


def _canvas_with_clean_mixed_structure() -> tuple[ChemusonCanvas, set[int], set[int]]:
    _ensure_app()
    canvas = ChemusonCanvas()
    phenyl, cyclopropane = _populate_clean_mixed_structure(canvas.model)
    canvas.state.bond_length = 40.0
    canvas._sync_scene_with_model()
    canvas.state.selected_atoms = set(canvas.model.atoms)
    return canvas, phenyl, cyclopropane


def _canvas_with_raw_mixed_structure() -> tuple[ChemusonCanvas, set[int], set[int]]:
    _ensure_app()
    canvas = ChemusonCanvas()
    phenyl, cyclopropane = _populate_raw_mixed_structure(canvas.model)
    canvas.state.bond_length = 40.0
    canvas._sync_scene_with_model()
    canvas.state.selected_atoms = set(canvas.model.atoms)
    return canvas, phenyl, cyclopropane


def _canvas_from_fixture(filename: str) -> ChemusonCanvas:
    _ensure_app()
    canvas = ChemusonCanvas()
    fixture = Path(__file__).parent / "fixtures" / filename
    PersistenceManager.load_from_file(str(fixture), canvas)
    canvas.state.bond_length = 40.0
    canvas.state.selected_atoms = set(canvas.model.atoms)
    return canvas


def _window_for_canvas(canvas: ChemusonCanvas) -> MagicMock:
    window = MagicMock()
    window.canvas = canvas
    status = MagicMock()
    window.statusBar.return_value = status
    return window


def _coords(graph: MolGraph) -> dict[int, tuple[float, float]]:
    return {atom_id: (atom.x, atom.y) for atom_id, atom in graph.atoms.items()}


def _last_status(window: MagicMock) -> str:
    return str(window.statusBar.return_value.showMessage.call_args.args[0])


def _lengths(graph: MolGraph, coords: dict[int, tuple[float, float]]) -> list[float]:
    return [
        math.hypot(
            coords[bond.a2_id][0] - coords[bond.a1_id][0],
            coords[bond.a2_id][1] - coords[bond.a1_id][1],
        )
        for bond in graph.bonds.values()
    ]


def _mean_displacement(
    before: dict[int, tuple[float, float]],
    after: dict[int, tuple[float, float]],
) -> float:
    common = set(before) & set(after)
    if not common:
        return 0.0
    return sum(math.dist(before[atom_id], after[atom_id]) for atom_id in common) / len(common)


def _candidate_debug_rows(result) -> list[dict]:
    keys = (
        "source",
        "rejected",
        "reason",
        "score",
        "quality_class",
        "quality_reason",
        "visual_score",
        "length_rms_error",
        "angle_rms_deviation",
        "exocyclic_ring_angle_min",
        "bad_exocyclic_ring_count",
        "geometry_hash",
        "novelty",
    )
    return [{key: row.get(key) for key in keys} for row in summarize_clean2d_candidates(result)]


def _chemical_signature(graph: MolGraph) -> tuple[tuple[tuple[int, dict], ...], tuple[tuple[int, dict], ...]]:
    atoms = []
    for atom_id, atom in sorted(graph.atoms.items()):
        data = asdict(atom)
        data.pop("x")
        data.pop("y")
        atoms.append((atom_id, data))
    bonds = [(bond_id, asdict(bond)) for bond_id, bond in sorted(graph.bonds.items())]
    return tuple(atoms), tuple(bonds)


def _distort_mixed_structure(graph: MolGraph) -> None:
    shifts = {
        7: (45.0, -55.0),
        8: (95.0, 35.0),
        9: (130.0, -45.0),
        10: (165.0, 80.0),
        11: (70.0, -65.0),
        12: (115.0, 75.0),
        13: (150.0, -15.0),
        14: (120.0, 95.0),
    }
    for atom_id, (dx, dy) in shifts.items():
        atom = graph.get_atom(atom_id)
        atom.x += dx
        atom.y += dy


def _nudge_to_suboptimal_geometry(graph: MolGraph) -> None:
    graph.atoms[8].x += 8.0
    graph.atoms[8].y += 5.0
    graph.atoms[9].x += 12.0


def test_smart_clean_first_press_canonicalizes_raw_mixed_structure() -> None:
    canvas, _phenyl, _cyclopropane = _canvas_with_raw_mixed_structure()
    try:
        window = _window_for_canvas(canvas)
        controller = Clean2DController()
        before = _coords(canvas.model)
        before_quality = classify_clean2d_layout_quality(canvas.graph, target_bond_length=40.0)

        quick = run_clean2d_engine(canvas.graph, set(canvas.model.atoms), mode="quick", target_bond_length=40.0)
        assert quick.selected is not None
        assert quick.selected.source != "current"
        assert before_quality.quality_class == "needs_rebuild"

        controller.run_clean_2d(window, 1.0, 200, "(test)", mode="quick")

        after = _coords(canvas.model)
        after_quality = classify_clean2d_layout_quality(canvas.graph, target_bond_length=40.0)
        assert after_quality.quality_class == "good"
        assert after_quality.visual_score < before_quality.visual_score
        assert max(_lengths(canvas.model, after)) < max(_lengths(canvas.model, before))
        assert "optimizada" in _last_status(window) or "limpiada" in _last_status(window)
        assert "alternativa propuesta" not in _last_status(window)
        assert canvas.undo_stack.count() == 1
    finally:
        canvas.close()


def test_smart_clean_second_press_proposes_after_canonicalization() -> None:
    canvas, _phenyl, _cyclopropane = _canvas_with_raw_mixed_structure()
    try:
        window = _window_for_canvas(canvas)
        controller = Clean2DController()

        controller.run_clean_2d(window, 1.0, 200, "(test)", mode="quick")
        canonical_hash = clean2d_geometry_hash(canvas.graph, _coords(canvas.model))
        assert classify_clean2d_layout_quality(canvas.graph, target_bond_length=40.0).quality_class == "good"

        controller.run_clean_2d(window, 1.0, 200, "(test)", mode="quick")

        assert clean2d_geometry_hash(canvas.graph, _coords(canvas.model)) != canonical_hash
        assert "alternativa propuesta" in _last_status(window)
        assert canvas.undo_stack.count() == 2
    finally:
        canvas.close()


def test_ctrl_k_on_real_crude_chemuson_file_canonicalizes_first() -> None:
    canvas = _canvas_from_fixture("test_clean2d_crude_chemuson.cmsn")
    try:
        window = _window_for_canvas(canvas)
        controller = Clean2DController()
        atom_ids = set(canvas.model.atoms)
        before = _coords(canvas.model)
        before_snapshot = capture_clean2d_snapshot(canvas.model)
        before_signature = _chemical_signature(canvas.model)

        engine_result = run_clean2d_engine(
            canvas.graph,
            atom_ids,
            mode="quick",
            target_bond_length=40.0,
        )
        debug_rows = _candidate_debug_rows(engine_result)
        assert engine_result.selected is not None, debug_rows
        assert engine_result.selected.source != "current", debug_rows

        controller.run_clean_2d(window, 1.0, 200, "(test)", mode="quick")

        after = _coords(canvas.model)
        after_quality = classify_clean2d_layout_quality(
            canvas.graph,
            atom_ids,
            target_bond_length=40.0,
        )
        status = _last_status(window)
        assert "optimizada" in status or "limpiada" in status
        assert "alternativa propuesta" not in status
        assert after_quality.quality_class == "good"
        assert _mean_displacement(before, after) >= 40.0 * 0.05
        assert_clean2d_invariants(before_snapshot, canvas.model, before, after)
        assert _chemical_signature(canvas.model) == before_signature
        assert canvas.undo_stack.count() == 1

        canvas.undo_stack.undo()
        assert _coords(canvas.model) == pytest.approx(before)
        canvas.undo_stack.redo()
        assert _coords(canvas.model) == pytest.approx(after)
    finally:
        canvas.close()


def test_smart_clean_repairs_distorted_structure_before_proposing() -> None:
    canvas, phenyl, cyclopropane = _canvas_with_clean_mixed_structure()
    try:
        _distort_mixed_structure(canvas.model)
        canvas._sync_scene_with_model()
        window = _window_for_canvas(canvas)
        controller = Clean2DController()
        distorted = _coords(canvas.model)
        before_lengths = _lengths(canvas.model, distorted)

        controller.run_clean_2d(window, 1.0, 200, "(test)", mode="quick")

        after = _coords(canvas.model)
        after_lengths = _lengths(canvas.model, after)
        assert max(before_lengths) > 40.0 * 2.0
        assert min(after_lengths) >= 40.0 * 0.70
        assert max(after_lengths) <= 40.0 * 1.45
        assert ring_degeneracy_score(after, phenyl) > 0.18
        assert ring_degeneracy_score(after, cyclopropane) > 0.05
        assert "optimizada" in _last_status(window) or "limpiada" in _last_status(window)
        assert "alternativa propuesta" not in _last_status(window)
    finally:
        canvas.close()


def test_successive_smart_cleans_avoid_immediate_geometry_repeat() -> None:
    canvas, _phenyl, _cyclopropane = _canvas_with_clean_mixed_structure()
    try:
        window = _window_for_canvas(canvas)
        controller = Clean2DController()
        initial_hash = clean2d_geometry_hash(canvas.graph, _coords(canvas.model))

        controller.run_clean_2d(window, 1.0, 200, "(test)", mode="quick")
        first_hash = clean2d_geometry_hash(canvas.graph, _coords(canvas.model))
        controller.run_clean_2d(window, 1.0, 200, "(test)", mode="quick")
        second_hash = clean2d_geometry_hash(canvas.graph, _coords(canvas.model))

        assert first_hash != initial_hash
        assert second_hash != first_hash
        assert canvas.undo_stack.count() == 2
    finally:
        canvas.close()


def test_smart_clean_preserves_chemical_invariants_and_metadata() -> None:
    canvas, _phenyl, _cyclopropane = _canvas_with_raw_mixed_structure()
    try:
        window = _window_for_canvas(canvas)
        controller = Clean2DController()
        before_coords = _coords(canvas.model)
        before_snapshot = capture_clean2d_snapshot(canvas.model)
        before_signature = _chemical_signature(canvas.model)

        controller.run_clean_2d(window, 1.0, 200, "(test)", mode="quick")

        after_coords = _coords(canvas.model)
        assert_clean2d_invariants(before_snapshot, canvas.model, before_coords, after_coords)
        assert _chemical_signature(canvas.model) == before_signature
    finally:
        canvas.close()


def test_smart_clean_alternative_has_safe_geometry() -> None:
    canvas, phenyl, cyclopropane = _canvas_with_clean_mixed_structure()
    try:
        window = _window_for_canvas(canvas)
        controller = Clean2DController()
        before = _coords(canvas.model)

        controller.run_clean_2d(window, 1.0, 200, "(test)", mode="quick")

        after = _coords(canvas.model)
        lengths = _lengths(canvas.model, after)
        assert all(math.isfinite(value) for xy in after.values() for value in xy)
        assert sum(lengths) / len(lengths) == pytest.approx(40.0, rel=0.10)
        assert min(lengths) >= 40.0 * 0.70
        assert max(lengths) <= 40.0 * 1.45
        assert count_new_bond_crossings(before, after, list(canvas.model.bonds.values())) == 0
        assert ring_degeneracy_score(after, phenyl) > 0.18
        assert ring_degeneracy_score(after, cyclopropane) > 0.05
    finally:
        canvas.close()


def test_smart_clean_application_is_undoable_and_redoable() -> None:
    canvas, _phenyl, _cyclopropane = _canvas_with_raw_mixed_structure()
    try:
        window = _window_for_canvas(canvas)
        controller = Clean2DController()
        before = _coords(canvas.model)

        controller.run_clean_2d(window, 1.0, 200, "(test)", mode="quick")
        optimized = _coords(canvas.model)
        controller.run_clean_2d(window, 1.0, 200, "(test)", mode="quick")
        proposed = _coords(canvas.model)

        assert optimized != before
        assert proposed != optimized
        canvas.undo_stack.undo()
        assert _coords(canvas.model) == pytest.approx(optimized)
        canvas.undo_stack.undo()
        assert _coords(canvas.model) == pytest.approx(before)
        canvas.undo_stack.redo()
        assert _coords(canvas.model) == pytest.approx(optimized)
        canvas.undo_stack.redo()
        assert _coords(canvas.model) == pytest.approx(proposed)
    finally:
        canvas.close()


def test_propose_rejects_suboptimal_base_instead_of_dragging_bad_geometry() -> None:
    graph = MolGraph()
    _populate_clean_mixed_structure(graph)
    _nudge_to_suboptimal_geometry(graph)
    quality = classify_clean2d_layout_quality(graph, target_bond_length=40.0)
    assert quality.quality_class == "needs_polish"

    propose = run_clean2d_engine(graph, mode="propose", target_bond_length=40.0)
    quick = run_clean2d_engine(graph, mode="quick", target_bond_length=40.0)

    assert not propose.ok
    assert propose.message == "La estructura debe optimizarse antes de proponer conformeros 2D"
    assert quick.ok
    assert quick.selected is not None
    canonical_quality = classify_clean2d_layout_quality(
        graph,
        coords=quick.selected.coords,
        target_bond_length=40.0,
    )
    assert canonical_quality.visual_score < quality.visual_score


def test_candidate_exit_gate_rejects_motion_without_canonicalization() -> None:
    graph = MolGraph()
    _populate_clean_mixed_structure(graph)
    _nudge_to_suboptimal_geometry(graph)
    before = _coords(graph)
    assert classify_clean2d_layout_quality(graph, target_bond_length=40.0).quality_class == "needs_polish"
    moved_same_shape = {
        atom_id: (x + 12.0, y - 8.0)
        for atom_id, (x, y) in before.items()
    }
    moving_candidate = Clean2DCandidate(
        source="rdkit_isolated",
        coords=moved_same_shape,
        message="bad motion",
    )

    result = rank_clean2d_candidates(
        graph,
        [moving_candidate],
        before,
        set(graph.atoms),
        mode="quick",
        target_bond_length=40.0,
    )
    summary = summarize_clean2d_candidates(result)

    assert result.selected is None
    assert any(row["source"] == "rdkit_isolated" and row["reason"] == "candidato_no_mejora_suficiente" for row in summary)
    assert all("quality_class" in row and "angle_rms_deviation" in row for row in summary)


def test_best_candidate_chosen_by_final_quality_not_source_priority() -> None:
    graph = MolGraph()
    _populate_clean_mixed_structure(graph)
    good_graph = MolGraph()
    _populate_clean_mixed_structure(good_graph)
    _nudge_to_suboptimal_geometry(graph)
    before = _coords(graph)
    good_coords = _coords(good_graph)
    assert classify_clean2d_layout_quality(graph, target_bond_length=40.0).quality_class == "needs_polish"

    high_priority_bad = Clean2DCandidate(
        source="rdkit_isolated",
        coords={atom_id: (x + 6.0, y + 6.0) for atom_id, (x, y) in before.items()},
        message="high priority but still bad",
    )
    low_priority_good = Clean2DCandidate(
        source="safe_fallback",
        coords=good_coords,
        message="low priority good",
    )

    result = rank_clean2d_candidates(
        graph,
        [high_priority_bad, low_priority_good],
        before,
        set(graph.atoms),
        mode="quick",
        target_bond_length=40.0,
    )

    assert result.selected is not None
    assert result.selected.source == "safe_fallback"
    assert result.selected.metadata["quality_class"] == "good"


def test_smart_clean_reports_clear_noop_when_no_safe_alternative_exists() -> None:
    _ensure_app()
    canvas = ChemusonCanvas()
    try:
        canvas.model.add_atom("He", 0.0, 0.0)
        canvas.state.bond_length = 40.0
        canvas._sync_scene_with_model()
        window = _window_for_canvas(canvas)
        controller = Clean2DController()

        controller.run_clean_2d(window, 1.0, 200, "(test)", mode="quick")

        assert canvas.undo_stack.count() == 0
        assert _last_status(window) == "Estructura 2D ya estaba limpia; no se encontró alternativa segura"
    finally:
        canvas.close()
