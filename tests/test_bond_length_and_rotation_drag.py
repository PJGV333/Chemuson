"""Regresiones para longitud visual de enlace y rotación por arrastre."""

from __future__ import annotations

import os
import sys

import pytest
from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.commands import ChangeBondLengthCommand


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def test_change_bond_length_updates_visible_bond_geometry() -> None:
    canvas = ChemusonCanvas()

    atom_a = canvas.model.add_atom("C", 100.0, 100.0)
    atom_b = canvas.model.add_atom("C", 140.0, 100.0)
    bond = canvas.model.add_bond(atom_a.id, atom_b.id, order=1)
    canvas._rebuild_items_from_model()

    item = canvas.bond_items[bond.id]
    assert item.path().boundingRect().width() == pytest.approx(40.0, abs=0.2)

    canvas.undo_stack.push(ChangeBondLengthCommand(canvas.model, canvas, bond.id, 10.0))

    updated_item = canvas.bond_items[bond.id]
    assert canvas.model.get_bond(bond.id).length_px == pytest.approx(10.0)
    assert updated_item.path().boundingRect().width() == pytest.approx(10.0, abs=0.2)
    assert (canvas.model.get_atom(atom_a.id).x, canvas.model.get_atom(atom_a.id).y) == pytest.approx(
        (100.0, 100.0)
    )
    assert (canvas.model.get_atom(atom_b.id).x, canvas.model.get_atom(atom_b.id).y) == pytest.approx(
        (140.0, 100.0)
    )

    canvas.undo_stack.undo()

    restored_item = canvas.bond_items[bond.id]
    assert canvas.model.get_bond(bond.id).length_px is None
    assert restored_item.path().boundingRect().width() == pytest.approx(40.0, abs=0.2)


def test_change_bond_length_repositions_terminal_explicit_label_without_overlap() -> None:
    canvas = ChemusonCanvas()

    atom_a = canvas.model.add_atom("N", 100.0, 100.0, is_explicit=True)
    atom_h = canvas.model.add_atom("H", 140.0, 100.0, is_explicit=True)
    bond = canvas.model.add_bond(atom_a.id, atom_h.id, order=1)
    canvas._rebuild_items_from_model()

    label_before = canvas.atom_items[atom_h.id].label.sceneBoundingRect().center()
    assert label_before.x() == pytest.approx(140.0, abs=1.5)

    canvas.undo_stack.push(ChangeBondLengthCommand(canvas.model, canvas, bond.id, 10.0))

    label_after = canvas.atom_items[atom_h.id].label.sceneBoundingRect().center()
    assert label_after.x() < label_before.x() - 2.0

    n_rect = canvas.atom_items[atom_a.id].label.sceneBoundingRect().adjusted(0.5, 0.5, -0.5, -0.5)
    h_rect = canvas.atom_items[atom_h.id].label.sceneBoundingRect().adjusted(0.5, 0.5, -0.5, -0.5)
    assert not n_rect.intersects(h_rect)

    canvas.undo_stack.undo()

    label_restored = canvas.atom_items[atom_h.id].label.sceneBoundingRect().center()
    assert label_restored.x() == pytest.approx(140.0, abs=1.5)


def test_rotation_drag_pushes_undo_and_redo_for_selected_fragment() -> None:
    canvas = ChemusonCanvas()

    atom_a = canvas.model.add_atom("C", 100.0, 100.0)
    atom_b = canvas.model.add_atom("C", 140.0, 100.0)
    bond = canvas.model.add_bond(atom_a.id, atom_b.id, order=1)
    canvas._rebuild_items_from_model()

    canvas.bond_items[bond.id].setSelected(True)
    canvas._sync_selection_from_scene()

    center = canvas._selected_items_bbox().center()
    canvas._begin_rotation_drag(QPointF(center.x() + 20.0, center.y()))
    canvas._update_rotation_drag(QPointF(center.x(), center.y() + 20.0))
    canvas._finalize_rotation_drag()

    assert canvas.undo_stack.count() == 1
    assert (canvas.model.get_atom(atom_a.id).x, canvas.model.get_atom(atom_a.id).y) == pytest.approx(
        (120.0, 80.0),
        abs=0.3,
    )
    assert (canvas.model.get_atom(atom_b.id).x, canvas.model.get_atom(atom_b.id).y) == pytest.approx(
        (120.0, 120.0),
        abs=0.3,
    )

    canvas.undo_stack.undo()

    assert (canvas.model.get_atom(atom_a.id).x, canvas.model.get_atom(atom_a.id).y) == pytest.approx(
        (100.0, 100.0),
        abs=0.2,
    )
    assert (canvas.model.get_atom(atom_b.id).x, canvas.model.get_atom(atom_b.id).y) == pytest.approx(
        (140.0, 100.0),
        abs=0.2,
    )

    canvas.undo_stack.redo()

    assert (canvas.model.get_atom(atom_a.id).x, canvas.model.get_atom(atom_a.id).y) == pytest.approx(
        (120.0, 80.0),
        abs=0.3,
    )
    assert (canvas.model.get_atom(atom_b.id).x, canvas.model.get_atom(atom_b.id).y) == pytest.approx(
        (120.0, 120.0),
        abs=0.3,
    )
