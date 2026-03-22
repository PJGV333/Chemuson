"""Regresiones para dibujo de corchetes según el rectángulo arrastrado."""

import os
import sys

import pytest
from PyQt6.QtCore import QPointF, QRectF
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.gui.canvas import ChemusonCanvas


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def _assert_path_stays_within_rect(item, rect: QRectF, tol: float = 0.01) -> None:
    """La traza visible no debe salirse del rectángulo dibujado."""
    bounds = item.path().boundingRect()
    assert bounds.left() >= rect.left() - tol
    assert bounds.top() >= rect.top() - tol
    assert bounds.right() <= rect.right() + tol
    assert bounds.bottom() <= rect.bottom() + tol


@pytest.mark.parametrize("kind", ["[]", "()", "{}"])
def test_bracket_tool_uses_dragged_rect_even_with_structure_inside(kind):
    """La herramienta debe respetar el rectángulo arrastrado, no el bbox del contenido."""
    canvas = ChemusonCanvas()

    a1 = canvas.model.add_atom("C", 120.0, 120.0).id
    a2 = canvas.model.add_atom("C", 160.0, 120.0).id
    canvas.model.add_bond(a1, a2, order=1)
    canvas._rebuild_items_from_model()

    drag_rect = QRectF(100.0, 90.0, 120.0, 80.0)
    canvas.set_current_tool("tool_brackets")
    canvas.state.active_bracket_type = kind
    canvas._bracket_drag_start = drag_rect.topLeft()
    canvas._last_scene_pos = drag_rect.bottomRight()

    canvas.mouseReleaseEvent(object())

    expected_count = 2 if len(kind) == 2 else 1
    assert len(canvas.bracket_items) == expected_count
    for item in canvas.bracket_items:
        base = item.base_rect()
        assert base.left() == pytest.approx(drag_rect.left())
        assert base.top() == pytest.approx(drag_rect.top())
        assert base.width() == pytest.approx(drag_rect.width())
        assert base.height() == pytest.approx(drag_rect.height())
        assert float(item._padding) == pytest.approx(0.0)
        _assert_path_stays_within_rect(item, drag_rect)


def test_bracket_tool_keeps_dragged_rect_when_structure_bbox_is_offset():
    """Aunque la estructura quede solo en una esquina, el bracket debe usar el drag completo."""
    canvas = ChemusonCanvas()

    a1 = canvas.model.add_atom("C", 210.0, 210.0).id
    a2 = canvas.model.add_atom("C", 245.0, 210.0).id
    canvas.model.add_bond(a1, a2, order=1)
    canvas._rebuild_items_from_model()

    drag_rect = QRectF(QPointF(180.0, 170.0), QPointF(320.0, 300.0)).normalized()
    canvas.set_current_tool("tool_brackets")
    canvas.state.active_bracket_type = "[]"
    canvas._bracket_drag_start = drag_rect.topLeft()
    canvas._last_scene_pos = drag_rect.bottomRight()

    canvas.mouseReleaseEvent(object())

    assert len(canvas.bracket_items) == 2
    left, right = canvas.bracket_items
    assert left.base_rect().width() == pytest.approx(drag_rect.width())
    assert left.base_rect().height() == pytest.approx(drag_rect.height())
    assert right.base_rect().width() == pytest.approx(drag_rect.width())
    assert right.base_rect().height() == pytest.approx(drag_rect.height())
    for item in (left, right):
        assert float(item._padding) == pytest.approx(0.0)
        _assert_path_stays_within_rect(item, drag_rect)
