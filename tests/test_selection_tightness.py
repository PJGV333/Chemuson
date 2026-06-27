"""Regresiones para selección rectangular ceñida."""


import pytest
from PyQt6.QtCore import QPoint, QPointF, Qt
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QApplication


from chemuson.core.model import BondStereo, BondStyle
from chemuson.gui import canvas as canvas_module
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


def test_single_click_selects_implicit_carbon_with_hover_hit_radius():
    """El clic simple debe dejar feedback visible al seleccionar un carbono implícito."""
    canvas = ChemusonCanvas()
    atom = canvas.model.add_atom("C", 120.0, 120.0)
    neighbor = canvas.model.add_atom("C", 160.0, 120.0)
    canvas.model.add_bond(atom.id, neighbor.id)
    canvas._rebuild_items_from_model()

    scene_pos = QPointF(atom.x + 14.0, atom.y)
    assert canvas._pick_hover_target(scene_pos) == (atom.id, None)
    assert canvas.atom_items[atom.id].label.isVisible() is False

    canvas.resize(900, 700)
    canvas.show()
    QApplication.processEvents()

    QTest.mouseClick(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        canvas.mapFromScene(scene_pos),
    )
    QApplication.processEvents()

    assert canvas.state.selected_atoms == {atom.id}
    assert canvas.state.selected_bonds == set()
    assert canvas.atom_items[atom.id].pen().style() != Qt.PenStyle.NoPen
    assert canvas._selected_items_bbox() is not None


def test_double_click_on_implicit_carbon_opens_atom_label_dialog(monkeypatch):
    """El doble clic sobre carbono implícito debe abrir el editor de elemento, no el de enlace."""
    canvas = ChemusonCanvas()
    atom = canvas.model.add_atom("C", 120.0, 120.0)
    neighbor = canvas.model.add_atom("C", 160.0, 120.0)
    canvas.model.add_bond(atom.id, neighbor.id)
    canvas._rebuild_items_from_model()

    scene_pos = QPointF(atom.x + 14.0, atom.y)
    assert canvas._pick_hover_target(scene_pos) == (atom.id, None)

    prompts: list[dict] = []

    def _fake_prompt(current_label, atom_id=None, initial=None):
        prompts.append(
            {
                "current_label": current_label,
                "atom_id": atom_id,
                "initial": initial,
            }
        )
        return None, None

    monkeypatch.setattr(canvas, "_prompt_atom_label", _fake_prompt)
    monkeypatch.setattr(
        canvas_module.QInputDialog,
        "getDouble",
        lambda *args, **kwargs: pytest.fail("No debe abrir el diálogo de longitud de enlace"),
    )

    canvas.resize(900, 700)
    canvas.show()
    QApplication.processEvents()

    QTest.mouseDClick(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        canvas.mapFromScene(scene_pos),
    )
    QApplication.processEvents()

    assert len(prompts) == 1
    assert prompts[0]["atom_id"] == atom.id
