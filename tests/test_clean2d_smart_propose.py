from __future__ import annotations

from dataclasses import asdict
import math
import os
import sys
from unittest.mock import MagicMock

import pytest
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.clean2d import (
    assert_clean2d_invariants,
    capture_clean2d_snapshot,
    clean2d_geometry_hash,
    count_new_bond_crossings,
    ring_degeneracy_score,
    run_clean2d_engine,
)
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


def _canvas_with_clean_mixed_structure() -> tuple[ChemusonCanvas, set[int], set[int]]:
    _ensure_app()
    canvas = ChemusonCanvas()
    phenyl, cyclopropane = _populate_clean_mixed_structure(canvas.model)
    canvas.state.bond_length = 40.0
    canvas._sync_scene_with_model()
    canvas.state.selected_atoms = set(canvas.model.atoms)
    return canvas, phenyl, cyclopropane


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


def test_smart_clean_proposes_alternative_for_initially_clean_mixed_structure() -> None:
    canvas, _phenyl, _cyclopropane = _canvas_with_clean_mixed_structure()
    try:
        window = _window_for_canvas(canvas)
        controller = Clean2DController()
        before = _coords(canvas.model)
        before_hash = clean2d_geometry_hash(canvas.graph, before)

        quick = run_clean2d_engine(canvas.graph, set(canvas.model.atoms), mode="quick", target_bond_length=40.0)
        assert quick.selected is not None
        assert quick.selected.source == "current"

        controller.run_clean_2d(window, 1.0, 200, "(test)", mode="quick")

        after = _coords(canvas.model)
        assert clean2d_geometry_hash(canvas.graph, after) != before_hash
        assert "alternativa propuesta" in _last_status(window)
        assert canvas.undo_stack.count() == 1
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
        assert "limpiada" in _last_status(window)
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
    canvas, _phenyl, _cyclopropane = _canvas_with_clean_mixed_structure()
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
    canvas, _phenyl, _cyclopropane = _canvas_with_clean_mixed_structure()
    try:
        window = _window_for_canvas(canvas)
        controller = Clean2DController()
        before = _coords(canvas.model)

        controller.run_clean_2d(window, 1.0, 200, "(test)", mode="quick")
        after = _coords(canvas.model)

        assert after != before
        canvas.undo_stack.undo()
        assert _coords(canvas.model) == pytest.approx(before)
        canvas.undo_stack.redo()
        assert _coords(canvas.model) == pytest.approx(after)
    finally:
        canvas.close()


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
