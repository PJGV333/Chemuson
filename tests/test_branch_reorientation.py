"""Regresiones para reorientación de ramas alrededor de enlaces pivote."""

from __future__ import annotations

import math

import pytest
from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QApplication


from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.geom import angle_deg, endpoint_from_angle_len


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def _distance(canvas: ChemusonCanvas, atom_a: int, atom_b: int) -> float:
    a = canvas.model.get_atom(atom_a)
    b = canvas.model.get_atom(atom_b)
    return math.hypot(a.x - b.x, a.y - b.y)


def test_invert_selected_branch_uses_smaller_side_and_preserves_lengths() -> None:
    canvas = ChemusonCanvas()

    atom_a = canvas.model.add_atom("C", 80.0, 100.0)
    atom_b = canvas.model.add_atom("C", 120.0, 100.0)
    atom_c = canvas.model.add_atom("C", 160.0, 100.0)
    atom_d = canvas.model.add_atom("C", 200.0, 60.0)
    atom_e = canvas.model.add_atom("C", 240.0, 60.0)

    canvas.model.add_bond(atom_a.id, atom_b.id, order=1)
    canvas.model.add_bond(atom_b.id, atom_c.id, order=1)
    bond_cd = canvas.model.add_bond(atom_c.id, atom_d.id, order=1)
    canvas.model.add_bond(atom_d.id, atom_e.id, order=1)
    canvas._rebuild_items_from_model()

    before_cd = _distance(canvas, atom_c.id, atom_d.id)
    before_de = _distance(canvas, atom_d.id, atom_e.id)

    canvas.bond_items[bond_cd.id].setSelected(True)
    canvas._sync_selection_from_scene()

    ok, _message = canvas.invert_selected_branch()

    assert ok
    assert canvas.undo_stack.count() == 1
    assert (canvas.model.get_atom(atom_a.id).x, canvas.model.get_atom(atom_a.id).y) == pytest.approx(
        (80.0, 100.0)
    )
    assert (canvas.model.get_atom(atom_b.id).x, canvas.model.get_atom(atom_b.id).y) == pytest.approx(
        (120.0, 100.0)
    )
    assert (canvas.model.get_atom(atom_c.id).x, canvas.model.get_atom(atom_c.id).y) == pytest.approx(
        (160.0, 100.0)
    )
    assert (canvas.model.get_atom(atom_d.id).x, canvas.model.get_atom(atom_d.id).y) == pytest.approx(
        (120.0, 140.0)
    )
    assert (canvas.model.get_atom(atom_e.id).x, canvas.model.get_atom(atom_e.id).y) == pytest.approx(
        (80.0, 140.0)
    )
    assert _distance(canvas, atom_c.id, atom_d.id) == pytest.approx(before_cd)
    assert _distance(canvas, atom_d.id, atom_e.id) == pytest.approx(before_de)

    canvas.undo_stack.undo()

    assert (canvas.model.get_atom(atom_d.id).x, canvas.model.get_atom(atom_d.id).y) == pytest.approx(
        (200.0, 60.0)
    )
    assert (canvas.model.get_atom(atom_e.id).x, canvas.model.get_atom(atom_e.id).y) == pytest.approx(
        (240.0, 60.0)
    )


def test_branch_rotation_snaps_to_next_chemical_orientation() -> None:
    canvas = ChemusonCanvas()

    center = canvas.model.add_atom("C", 200.0, 200.0)
    carbonyl_o = canvas.model.add_atom("O", 240.0, 200.0)
    moving = canvas.model.add_atom("N", 180.0, 165.359)
    tail = canvas.model.add_atom("C", 160.0, 130.718)

    canvas.model.add_bond(center.id, carbonyl_o.id, order=2)
    pivot = canvas.model.add_bond(center.id, moving.id, order=1)
    canvas.model.add_bond(moving.id, tail.id, order=1)
    canvas._rebuild_items_from_model()

    ok, _message = canvas.rotate_branch_side_degrees(pivot.id, center.id, 60.0)

    assert ok
    snapped_angle = angle_deg(
        QPointF(canvas.model.get_atom(center.id).x, canvas.model.get_atom(center.id).y),
        QPointF(canvas.model.get_atom(moving.id).x, canvas.model.get_atom(moving.id).y),
    )
    assert snapped_angle == pytest.approx(240.0, abs=0.4)
    assert _distance(canvas, center.id, moving.id) == pytest.approx(40.0, abs=0.2)
    assert _distance(canvas, moving.id, tail.id) == pytest.approx(40.0, abs=0.3)


def test_auto_arrange_branch_picks_less_congested_side() -> None:
    canvas = ChemusonCanvas()

    center = canvas.model.add_atom("C", 200.0, 200.0)
    carbonyl_o = canvas.model.add_atom("O", 240.0, 200.0)
    moving_pos = endpoint_from_angle_len(QPointF(200.0, 200.0), 120.0, 40.0)
    moving = canvas.model.add_atom("N", moving_pos.x(), moving_pos.y())
    tail_pos = endpoint_from_angle_len(moving_pos, 120.0, 40.0)
    tail = canvas.model.add_atom("C", tail_pos.x(), tail_pos.y())
    obstacle = canvas.model.add_atom("Cl", tail_pos.x() + 2.0, tail_pos.y() + 1.0)

    canvas.model.add_bond(center.id, carbonyl_o.id, order=2)
    pivot = canvas.model.add_bond(center.id, moving.id, order=1)
    canvas.model.add_bond(moving.id, tail.id, order=1)
    canvas._rebuild_items_from_model()

    before_angle = angle_deg(
        QPointF(canvas.model.get_atom(center.id).x, canvas.model.get_atom(center.id).y),
        QPointF(canvas.model.get_atom(moving.id).x, canvas.model.get_atom(moving.id).y),
    )
    assert before_angle == pytest.approx(120.0, abs=0.2)

    ok, _message = canvas.auto_arrange_branch_side(pivot.id, center.id)

    assert ok
    after_angle = angle_deg(
        QPointF(canvas.model.get_atom(center.id).x, canvas.model.get_atom(center.id).y),
        QPointF(canvas.model.get_atom(moving.id).x, canvas.model.get_atom(moving.id).y),
    )
    assert after_angle == pytest.approx(240.0, abs=0.4)
    obstacle_atom = canvas.model.get_atom(obstacle.id)
    tail_atom = canvas.model.get_atom(tail.id)
    assert math.hypot(tail_atom.x - obstacle_atom.x, tail_atom.y - obstacle_atom.y) > 20.0


def test_selected_branch_rotation_rejects_cyclic_bond() -> None:
    canvas = ChemusonCanvas()

    atom_a = canvas.model.add_atom("C", 100.0, 100.0)
    atom_b = canvas.model.add_atom("C", 140.0, 100.0)
    atom_c = canvas.model.add_atom("C", 120.0, 65.359)
    bond_ab = canvas.model.add_bond(atom_a.id, atom_b.id, order=1)
    canvas.model.add_bond(atom_b.id, atom_c.id, order=1)
    canvas.model.add_bond(atom_c.id, atom_a.id, order=1)
    canvas._rebuild_items_from_model()

    before = {
        atom.id: (atom.x, atom.y)
        for atom in (atom_a, atom_b, atom_c)
    }

    canvas.bond_items[bond_ab.id].setSelected(True)
    canvas._sync_selection_from_scene()

    ok, message = canvas.rotate_selected_branch_degrees(60.0)

    assert not ok
    assert "acíclico" in message
    after = {
        atom.id: (canvas.model.get_atom(atom.id).x, canvas.model.get_atom(atom.id).y)
        for atom in (atom_a, atom_b, atom_c)
    }
    assert after == pytest.approx(before)
