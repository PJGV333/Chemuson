"""Regresiones para selección rectangular ceñida."""

import os
import sys

import pytest
from PyQt6.QtCore import QPoint, QPointF, Qt
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.core.model import BondStereo, BondStyle
from chemuson.gui.canvas import ChemusonCanvas


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def test_rect_selection_uses_item_centers_to_avoid_neighbor_leak():
    """Rectángulo de selección no debe arrastrar enlaces vecinos por solapamiento parcial."""
    canvas = ChemusonCanvas()

    a1 = canvas.model.add_atom("N", 120.0, 120.0).id
    a2 = canvas.model.add_atom("N", 160.0, 120.0).id
    a3 = canvas.model.add_atom("N", 174.0, 120.0).id
    a4 = canvas.model.add_atom("N", 214.0, 120.0).id
    b1 = canvas.model.add_bond(
        a1,
        a2,
        order=1,
        style=BondStyle.PLAIN,
        stereo=BondStereo.NONE,
        is_aromatic=False,
    ).id
    b2 = canvas.model.add_bond(
        a3,
        a4,
        order=1,
        style=BondStyle.PLAIN,
        stereo=BondStereo.NONE,
        is_aromatic=False,
    ).id
    canvas._rebuild_items_from_model()

    # El rectángulo toca una porción del segundo enlace, pero su centro queda fuera.
    canvas._begin_selection_drag(QPointF(108.0, 108.0), free_select=False, additive=False)
    canvas._last_scene_pos = QPointF(171.0, 132.0)
    canvas._finalize_selection_drag()

    assert b1 in canvas.state.selected_bonds
    assert b2 not in canvas.state.selected_bonds


def test_switching_from_selection_tool_clears_selection():
    """Al pasar de selección a otra herramienta, la selección activa debe limpiarse."""
    canvas = ChemusonCanvas()

    atom = canvas.model.add_atom("C", 120.0, 120.0)
    canvas._rebuild_items_from_model()
    canvas.atom_items[atom.id].setSelected(True)
    canvas._sync_selection_from_scene()
    assert atom.id in canvas.state.selected_atoms

    canvas.set_current_tool("tool_select")
    canvas.set_current_tool("tool_bond")

    assert canvas.state.selected_atoms == set()
    assert canvas.state.selected_bonds == set()
    assert not canvas.scene.selectedItems()


def test_switching_from_precise_rotation_tool_clears_selection_on_other_tool():
    """La selección usada por rotación precisa no debe quedarse al pasar a otra herramienta."""
    canvas = ChemusonCanvas()

    atom = canvas.model.add_atom("C", 120.0, 120.0)
    canvas._rebuild_items_from_model()
    canvas.atom_items[atom.id].setSelected(True)
    canvas._sync_selection_from_scene()
    assert atom.id in canvas.state.selected_atoms

    canvas.set_current_tool("tool_rotate_3d_precise")
    canvas.set_current_tool("tool_bond")

    assert canvas.state.selected_atoms == set()
    assert canvas.state.selected_bonds == set()
    assert not canvas.scene.selectedItems()


def test_precise_rotation_prompts_after_delimiting_selection(monkeypatch):
    """Rotación 3D precisa abre su diálogo al finalizar delimitación de selección."""
    canvas = ChemusonCanvas()

    atom = canvas.model.add_atom("C", 120.0, 120.0)
    canvas._rebuild_items_from_model()
    canvas.set_current_tool("tool_rotate_3d_precise")

    prompts: list[dict] = []
    monkeypatch.setattr(
        canvas,
        "_prompt_precise_3d_rotation",
        lambda *args, **kwargs: prompts.append(dict(kwargs)),
    )

    canvas._begin_selection_drag(QPointF(100.0, 100.0), free_select=False, additive=False)
    canvas._last_scene_pos = QPointF(140.0, 140.0)
    canvas.mouseReleaseEvent(object())

    assert atom.id in canvas.state.selected_atoms
    assert len(prompts) == 1


def test_precise_rotation_starts_free_selection_by_default(monkeypatch):
    """Rotación 3D precisa debe iniciar delimitación como selección libre."""
    canvas = ChemusonCanvas()
    canvas.set_current_tool("tool_rotate_3d_precise")

    calls: list[tuple[bool, bool]] = []
    monkeypatch.setattr(canvas, "_is_on_paper", lambda x, y: True)
    monkeypatch.setattr(canvas, "_get_item_at", lambda scene_pos: None)
    monkeypatch.setattr(
        canvas,
        "_begin_selection_drag",
        lambda scene_pos, free_select, additive: calls.append((bool(free_select), bool(additive))),
    )

    class _FakeMouseEvent:
        def pos(self):
            return QPoint(100, 100)

        def button(self):
            return Qt.MouseButton.LeftButton

        def modifiers(self):
            return Qt.KeyboardModifier.NoModifier

    canvas.mousePressEvent(_FakeMouseEvent())

    assert calls == [(True, False)]
