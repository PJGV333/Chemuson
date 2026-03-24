"""Regresiones para arrastre de selección con undo/redo."""

from __future__ import annotations

import os
import sys

import pytest
from PyQt6.QtCore import QPoint, Qt
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.gui.canvas import ChemusonCanvas


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def _drag_selection(canvas: ChemusonCanvas, start_scene, dx: int, dy: int) -> None:
    canvas.resize(900, 700)
    canvas.show()
    QApplication.processEvents()

    start_view = canvas.mapFromScene(start_scene)
    end_view = QPoint(start_view.x() + dx, start_view.y() + dy)

    QTest.mousePress(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        start_view,
    )
    QTest.mouseMove(canvas.viewport(), end_view)
    QTest.mouseRelease(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        end_view,
    )
    QApplication.processEvents()


def test_mouse_drag_pushes_single_undoable_move() -> None:
    canvas = ChemusonCanvas()
    atom_a = canvas.model.add_atom("C", 100.0, 100.0)
    atom_b = canvas.model.add_atom("C", 140.0, 100.0)
    canvas.model.add_bond(atom_a.id, atom_b.id)
    canvas._rebuild_items_from_model()

    canvas.atom_items[atom_a.id].setSelected(True)
    canvas.atom_items[atom_b.id].setSelected(True)
    canvas._sync_selection_from_scene()

    _drag_selection(canvas, canvas.atom_items[atom_a.id].scenePos(), 28, 16)

    moved_a = canvas.model.get_atom(atom_a.id)
    moved_b = canvas.model.get_atom(atom_b.id)
    assert canvas.undo_stack.count() == 1
    assert canvas.undo_stack.index() == 1
    assert (moved_a.x, moved_a.y) == pytest.approx((128.0, 116.0))
    assert (moved_b.x, moved_b.y) == pytest.approx((168.0, 116.0))

    canvas.undo_stack.undo()
    QApplication.processEvents()

    restored_a = canvas.model.get_atom(atom_a.id)
    restored_b = canvas.model.get_atom(atom_b.id)
    assert (restored_a.x, restored_a.y) == pytest.approx((100.0, 100.0))
    assert (restored_b.x, restored_b.y) == pytest.approx((140.0, 100.0))

    canvas.undo_stack.redo()
    QApplication.processEvents()

    redone_a = canvas.model.get_atom(atom_a.id)
    redone_b = canvas.model.get_atom(atom_b.id)
    assert (redone_a.x, redone_a.y) == pytest.approx((128.0, 116.0))
    assert (redone_b.x, redone_b.y) == pytest.approx((168.0, 116.0))


def test_drag_clears_live_state_before_undo_stack_index_changes() -> None:
    canvas = ChemusonCanvas()
    atom_a = canvas.model.add_atom("C", 110.0, 110.0)
    atom_b = canvas.model.add_atom("C", 150.0, 110.0)
    canvas.model.add_bond(atom_a.id, atom_b.id)
    canvas._rebuild_items_from_model()

    canvas.atom_items[atom_a.id].setSelected(True)
    canvas.atom_items[atom_b.id].setSelected(True)
    canvas._sync_selection_from_scene()

    states_during_push: list[bool] = []
    canvas.undo_stack.indexChanged.connect(
        lambda _index: states_during_push.append(bool(canvas._dragging_selection))
    )

    _drag_selection(canvas, canvas.atom_items[atom_a.id].scenePos(), 22, 12)

    assert states_during_push
    assert states_during_push == [False]


def test_live_drag_defers_expensive_refresh_until_release(monkeypatch) -> None:
    canvas = ChemusonCanvas()
    atom_a = canvas.model.add_atom("N", 100.0, 100.0, is_explicit=True)
    atom_b = canvas.model.add_atom("C", 140.0, 100.0, is_explicit=True)
    canvas.model.add_bond(atom_a.id, atom_b.id)
    canvas.state.show_implicit_hydrogens = True
    canvas._rebuild_items_from_model()

    canvas.atom_items[atom_a.id].setSelected(True)
    canvas.atom_items[atom_b.id].setSelected(True)
    canvas._sync_selection_from_scene()

    counts = {"implicit_h": 0, "numbering": 0, "bbox": 0}

    orig_implicit = canvas._refresh_implicit_h_overlays
    orig_numbering = canvas.recompute_numbering
    orig_bbox = canvas._selected_items_bbox

    def _count_implicit(*args, **kwargs):
        counts["implicit_h"] += 1
        return orig_implicit(*args, **kwargs)

    def _count_numbering(*args, **kwargs):
        counts["numbering"] += 1
        return orig_numbering(*args, **kwargs)

    def _count_bbox(*args, **kwargs):
        counts["bbox"] += 1
        return orig_bbox(*args, **kwargs)

    monkeypatch.setattr(canvas, "_refresh_implicit_h_overlays", _count_implicit)
    monkeypatch.setattr(canvas, "recompute_numbering", _count_numbering)
    monkeypatch.setattr(canvas, "_selected_items_bbox", _count_bbox)

    canvas.resize(900, 700)
    canvas.show()
    QApplication.processEvents()

    start_scene = canvas.atom_items[atom_a.id].scenePos()
    start_view = canvas.mapFromScene(start_scene)
    end_view = QPoint(start_view.x() + 24, start_view.y() + 14)

    QTest.mousePress(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        start_view,
    )
    QApplication.processEvents()
    after_press = dict(counts)

    QTest.mouseMove(canvas.viewport(), end_view)
    QApplication.processEvents()

    assert counts["implicit_h"] == after_press["implicit_h"]
    assert counts["numbering"] == after_press["numbering"]
    assert counts["bbox"] == after_press["bbox"]

    QTest.mouseRelease(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        end_view,
    )
    QApplication.processEvents()

    assert counts["implicit_h"] > after_press["implicit_h"]
    assert counts["numbering"] > after_press["numbering"]
    assert counts["bbox"] > after_press["bbox"]
