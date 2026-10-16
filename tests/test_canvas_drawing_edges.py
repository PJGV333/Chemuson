"""Qt regressions for molecular drawing at and beyond the visible paper edge."""
from __future__ import annotations

import math

import pytest
from PyQt6.QtCore import QPointF, Qt
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QApplication

from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.commands import AddAtomCommand
from chemuson.gui.templates import build_haworth_template


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def _show_canvas(canvas: ChemusonCanvas, app: QApplication) -> None:
    canvas.resize(900, 650)
    canvas.show()
    app.processEvents()
    canvas.center_on_paper()
    app.processEvents()


def _drag_scene(
    canvas: ChemusonCanvas,
    app: QApplication,
    start: QPointF,
    end: QPointF,
    *,
    allow_end_outside: bool = False,
) -> None:
    start_pos = canvas.mapFromScene(start)
    end_pos = canvas.mapFromScene(end)
    assert canvas.viewport().rect().contains(start_pos)
    if not allow_end_outside:
        assert canvas.viewport().rect().contains(end_pos)
    QTest.mousePress(canvas.viewport(), Qt.MouseButton.LeftButton, pos=start_pos)
    QTest.mouseMove(canvas.viewport(), end_pos, 10)
    QTest.mouseRelease(canvas.viewport(), Qt.MouseButton.LeftButton, pos=end_pos)
    app.processEvents()


@pytest.mark.parametrize(
    ("start_x", "end_x", "allow_end_outside"),
    ((20.0, 60.0, False), (-40.0, 0.0, False), (-40.0, -80.0, True)),
    ids=("inside-paper-margin", "visible-viewport-edge", "beyond-viewport"),
)
def test_bond_tool_can_start_near_or_beyond_paper_edge(
    start_x: float,
    end_x: float,
    allow_end_outside: bool,
) -> None:
    app = QApplication.instance() or QApplication([])
    canvas = ChemusonCanvas()
    _show_canvas(canvas, app)
    canvas.set_current_tool("tool_bond")

    _drag_scene(
        canvas,
        app,
        QPointF(start_x, 300.0),
        QPointF(end_x, 300.0),
        allow_end_outside=allow_end_outside,
    )

    assert len(canvas.model.atoms) == 2
    assert len(canvas.model.bonds) == 1
    assert canvas.undo_stack.count() == 1
    assert all(
        canvas.scene.sceneRect().contains(QPointF(atom.x, atom.y))
        for atom in canvas.model.atoms.values()
    )
    edge_snapshot = tuple(
        (bond.a1_id, bond.a2_id, bond.order, bond.style)
        for bond in canvas.model.bonds.values()
    )
    canvas.undo_stack.undo()
    app.processEvents()
    assert not canvas.model.atoms
    assert not canvas.model.bonds
    canvas.undo_stack.redo()
    app.processEvents()
    assert tuple(
        (bond.a1_id, bond.a2_id, bond.order, bond.style)
        for bond in canvas.model.bonds.values()
    ) == edge_snapshot
    canvas.close()
    app.processEvents()


def test_ring_tool_can_create_closed_cycle_at_viewport_edge_and_undo() -> None:
    app = QApplication.instance() or QApplication([])
    canvas = ChemusonCanvas()
    _show_canvas(canvas, app)
    canvas.set_current_tool("tool_ring")
    canvas.state.active_ring_size = 6

    _drag_scene(canvas, app, QPointF(-40.0, 300.0), QPointF(0.0, 300.0))

    assert len(canvas.model.atoms) == 6
    assert len(canvas.model.bonds) == 6
    degrees = {
        atom_id: sum(atom_id in (bond.a1_id, bond.a2_id) for bond in canvas.model.bonds.values())
        for atom_id in canvas.model.atoms
    }
    assert set(degrees.values()) == {2}
    edge_snapshot = tuple(
        (bond.a1_id, bond.a2_id, bond.order, bond.style)
        for bond in canvas.model.bonds.values()
    )
    canvas.undo_stack.undo()
    app.processEvents()
    assert not canvas.model.atoms
    canvas.undo_stack.redo()
    app.processEvents()
    assert len(canvas.model.atoms) == 6
    assert tuple(
        (bond.a1_id, bond.a2_id, bond.order, bond.style)
        for bond in canvas.model.bonds.values()
    ) == edge_snapshot
    canvas.close()
    app.processEvents()


def test_scene_rect_expands_for_outlying_molecule_and_undo_keeps_view() -> None:
    app = QApplication.instance() or QApplication([])
    canvas = ChemusonCanvas()
    _show_canvas(canvas, app)
    old_rect = canvas.scene.sceneRect()
    far_position = QPointF(old_rect.left() - 240.0, old_rect.center().y())

    canvas.undo_stack.push(
        AddAtomCommand(
            canvas.model,
            canvas,
            "C",
            far_position.x(),
            far_position.y(),
            is_explicit=False,
        )
    )
    app.processEvents()
    assert canvas.scene.sceneRect().contains(far_position)

    canvas.centerOn(far_position)
    app.processEvents()
    center_before_undo = canvas.mapToScene(canvas.viewport().rect().center())
    assert math.isclose(center_before_undo.x(), far_position.x(), abs_tol=2.0)
    assert math.isclose(center_before_undo.y(), far_position.y(), abs_tol=2.0)

    canvas.undo_stack.undo()
    app.processEvents()
    center_after_undo = canvas.mapToScene(canvas.viewport().rect().center())
    assert not canvas.model.atoms
    assert math.hypot(
        center_after_undo.x() - center_before_undo.x(),
        center_after_undo.y() - center_before_undo.y(),
    ) < 2.0
    canvas.close()
    app.processEvents()


def test_haworth_template_can_be_placed_at_viewport_edge_and_undone() -> None:
    app = QApplication.instance() or QApplication([])
    canvas = ChemusonCanvas()
    _show_canvas(canvas, app)
    canvas.set_current_tool("tool_select")
    canvas.begin_template_insert_mode(
        build_haworth_template(40.0, anomeric_up=True, bold_front=True),
        "Haworth β",
    )
    click_pos = canvas.mapFromScene(QPointF(-40.0, 300.0))
    assert canvas.viewport().rect().contains(click_pos)

    QTest.mouseClick(canvas.viewport(), Qt.MouseButton.LeftButton, pos=click_pos)
    app.processEvents()
    assert len(canvas.model.atoms) == 11
    assert len(canvas.model.bonds) == 11
    assert canvas._pending_template_graph is None
    edge_snapshot = tuple(
        (bond.a1_id, bond.a2_id, bond.order, bond.style, bond.is_aromatic)
        for bond in canvas.model.bonds.values()
    )

    canvas.undo_stack.undo()
    app.processEvents()
    assert not canvas.model.atoms
    canvas.undo_stack.redo()
    app.processEvents()
    assert len(canvas.model.atoms) == 11
    assert tuple(
        (bond.a1_id, bond.a2_id, bond.order, bond.style, bond.is_aromatic)
        for bond in canvas.model.bonds.values()
    ) == edge_snapshot
    canvas.close()
    app.processEvents()
