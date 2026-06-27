"""Regresiones para rotación rígida de fragmentos alrededor de un átomo pivote."""

from __future__ import annotations

import math

import pytest
from PyQt6.QtWidgets import QApplication


from chemuson.gui.canvas import ChemusonCanvas


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def _distance(canvas: ChemusonCanvas, atom_a: int, atom_b: int) -> float:
    a = canvas.model.get_atom(atom_a)
    b = canvas.model.get_atom(atom_b)
    return math.hypot(a.x - b.x, a.y - b.y)


def test_rotate_selected_fragment_around_pivot_preserves_geometry_and_undo() -> None:
    canvas = ChemusonCanvas()

    atom_a = canvas.model.add_atom("C", 80.0, 100.0)
    atom_b = canvas.model.add_atom("C", 120.0, 100.0)
    atom_c = canvas.model.add_atom("C", 160.0, 100.0)
    atom_d = canvas.model.add_atom("C", 200.0, 100.0)

    canvas.model.add_bond(atom_a.id, atom_b.id, order=1)
    canvas.model.add_bond(atom_b.id, atom_c.id, order=1)
    canvas.model.add_bond(atom_c.id, atom_d.id, order=1)
    canvas._rebuild_items_from_model()

    before_bc = _distance(canvas, atom_b.id, atom_c.id)
    before_cd = _distance(canvas, atom_c.id, atom_d.id)

    canvas.atom_items[atom_c.id].setSelected(True)
    canvas.atom_items[atom_d.id].setSelected(True)
    canvas._sync_selection_from_scene()

    ok, _message = canvas.rotate_selected_fragment_degrees(atom_b.id, 90.0)

    assert ok
    assert canvas.undo_stack.count() == 1
    assert (canvas.model.get_atom(atom_a.id).x, canvas.model.get_atom(atom_a.id).y) == pytest.approx(
        (80.0, 100.0)
    )
    assert (canvas.model.get_atom(atom_b.id).x, canvas.model.get_atom(atom_b.id).y) == pytest.approx(
        (120.0, 100.0)
    )
    assert (canvas.model.get_atom(atom_c.id).x, canvas.model.get_atom(atom_c.id).y) == pytest.approx(
        (120.0, 140.0)
    )
    assert (canvas.model.get_atom(atom_d.id).x, canvas.model.get_atom(atom_d.id).y) == pytest.approx(
        (120.0, 180.0)
    )
    assert _distance(canvas, atom_b.id, atom_c.id) == pytest.approx(before_bc)
    assert _distance(canvas, atom_c.id, atom_d.id) == pytest.approx(before_cd)

    canvas.undo_stack.undo()

    assert (canvas.model.get_atom(atom_c.id).x, canvas.model.get_atom(atom_c.id).y) == pytest.approx(
        (160.0, 100.0)
    )
    assert (canvas.model.get_atom(atom_d.id).x, canvas.model.get_atom(atom_d.id).y) == pytest.approx(
        (200.0, 100.0)
    )

    canvas.undo_stack.redo()

    assert (canvas.model.get_atom(atom_c.id).x, canvas.model.get_atom(atom_c.id).y) == pytest.approx(
        (120.0, 140.0)
    )
    assert (canvas.model.get_atom(atom_d.id).x, canvas.model.get_atom(atom_d.id).y) == pytest.approx(
        (120.0, 180.0)
    )


def test_rotate_selected_fragment_allows_pivot_to_be_part_of_selection() -> None:
    canvas = ChemusonCanvas()

    atom_a = canvas.model.add_atom("C", 100.0, 100.0)
    atom_b = canvas.model.add_atom("C", 140.0, 100.0)
    atom_c = canvas.model.add_atom("C", 180.0, 100.0)

    canvas.model.add_bond(atom_a.id, atom_b.id, order=1)
    canvas.model.add_bond(atom_b.id, atom_c.id, order=1)
    canvas._rebuild_items_from_model()

    canvas.atom_items[atom_a.id].setSelected(True)
    canvas.atom_items[atom_b.id].setSelected(True)
    canvas.atom_items[atom_c.id].setSelected(True)
    canvas._sync_selection_from_scene()

    ok, _message = canvas.rotate_selected_fragment_degrees(atom_a.id, 180.0)

    assert ok
    assert (canvas.model.get_atom(atom_a.id).x, canvas.model.get_atom(atom_a.id).y) == pytest.approx(
        (100.0, 100.0)
    )
    assert (canvas.model.get_atom(atom_b.id).x, canvas.model.get_atom(atom_b.id).y) == pytest.approx(
        (60.0, 100.0)
    )
    assert (canvas.model.get_atom(atom_c.id).x, canvas.model.get_atom(atom_c.id).y) == pytest.approx(
        (20.0, 100.0)
    )


def test_rotate_selected_fragment_rejects_external_attachment_outside_pivot() -> None:
    canvas = ChemusonCanvas()

    atom_a = canvas.model.add_atom("C", 80.0, 100.0)
    atom_b = canvas.model.add_atom("C", 120.0, 100.0)
    atom_c = canvas.model.add_atom("C", 160.0, 100.0)
    atom_d = canvas.model.add_atom("C", 200.0, 100.0)

    canvas.model.add_bond(atom_a.id, atom_b.id, order=1)
    canvas.model.add_bond(atom_b.id, atom_c.id, order=1)
    canvas.model.add_bond(atom_c.id, atom_d.id, order=1)
    canvas._rebuild_items_from_model()

    before = {
        atom.id: (atom.x, atom.y)
        for atom in (atom_a, atom_b, atom_c, atom_d)
    }

    canvas.atom_items[atom_c.id].setSelected(True)
    canvas._sync_selection_from_scene()

    ok, message = canvas.rotate_selected_fragment_degrees(atom_b.id, 90.0)

    assert not ok
    assert "pivote" in message
    after = {
        atom.id: (canvas.model.get_atom(atom.id).x, canvas.model.get_atom(atom.id).y)
        for atom in (atom_a, atom_b, atom_c, atom_d)
    }
    assert after == pytest.approx(before)
    assert canvas.undo_stack.count() == 0


def test_rotate_selected_fragment_rejects_other_disconnected_component() -> None:
    canvas = ChemusonCanvas()

    pivot = canvas.model.add_atom("C", 100.0, 100.0)
    attached = canvas.model.add_atom("C", 140.0, 100.0)
    other_a = canvas.model.add_atom("C", 240.0, 100.0)
    other_b = canvas.model.add_atom("C", 280.0, 100.0)

    canvas.model.add_bond(pivot.id, attached.id, order=1)
    canvas.model.add_bond(other_a.id, other_b.id, order=1)
    canvas._rebuild_items_from_model()

    canvas.atom_items[other_a.id].setSelected(True)
    canvas.atom_items[other_b.id].setSelected(True)
    canvas._sync_selection_from_scene()

    ok, message = canvas.rotate_selected_fragment_degrees(pivot.id, 90.0)

    assert not ok
    assert "pivote" in message
    assert canvas.undo_stack.count() == 0
