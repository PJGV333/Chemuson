"""Regresiones para radio de esferas de coordinación."""

from __future__ import annotations


import pytest
from PyQt6.QtWidgets import QApplication


from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.commands import ChangeCoordinationSphereStyleCommand


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def test_coordination_sphere_radius_command_updates_resin_visual_size():
    """El radio explícito debe reflejarse visualmente incluso con etiquetas largas."""
    canvas = ChemusonCanvas()
    atom = canvas.model.add_atom(
        "Resin",
        120.0,
        120.0,
        is_explicit=True,
        is_coordination_center=True,
        sphere_radius=16.0,
    )
    canvas._rebuild_items_from_model()

    item = canvas.atom_items[atom.id]
    assert item._coordination_draw_radius() == pytest.approx(16.0)

    canvas.undo_stack.push(
        ChangeCoordinationSphereStyleCommand(
            canvas.model,
            canvas,
            atom.id,
            new_radius=4.0,
        )
    )

    assert canvas.model.get_atom(atom.id).sphere_radius == pytest.approx(4.0)
    assert item._coordination_draw_radius() == pytest.approx(4.0)

    canvas.undo_stack.undo()

    assert canvas.model.get_atom(atom.id).sphere_radius == pytest.approx(16.0)
    assert item._coordination_draw_radius() == pytest.approx(16.0)
