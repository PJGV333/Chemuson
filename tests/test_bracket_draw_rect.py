"""Regresiones para dibujo de corchetes según el rectángulo arrastrado."""


import pytest
from PyQt6.QtCore import QPointF, QRectF
from PyQt6.QtWidgets import QApplication


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


@pytest.mark.parametrize(
    ("kind", "expected_count"),
    [
        ("[]", 2),
        ("[", 1),
        ("]", 1),
        ("()", 2),
        ("{}", 2),
        ("{", 1),
        ("}", 1),
        ("corner", 1),
        ("frame", 1),
        ("rounded_frame", 1),
    ],
)
def test_bracket_tool_uses_dragged_rect_even_with_structure_inside(kind, expected_count):
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

    assert len(canvas.bracket_items) == expected_count
    for item in canvas.bracket_items:
        base = item.base_rect()
        assert base.left() == pytest.approx(drag_rect.left())
        assert base.top() == pytest.approx(drag_rect.top())
        assert base.width() == pytest.approx(drag_rect.width())
        assert base.height() == pytest.approx(drag_rect.height())
        assert float(item._padding) == pytest.approx(0.0)
        _assert_path_stays_within_rect(item, drag_rect)


@pytest.mark.parametrize(
    ("tool_id", "expected_kind"),
    [
        ("tool_brackets_square", "[]"),
        ("tool_brackets_square_left", "["),
        ("tool_brackets_square_right", "]"),
        ("tool_brackets_corner", "corner"),
        ("tool_brackets_curly", "{}"),
        ("tool_brackets_curly_left", "{"),
        ("tool_brackets_curly_right", "}"),
        ("tool_brackets_frame", "frame"),
        ("tool_brackets_frame_rounded", "rounded_frame"),
        ("tool_brackets_round", "()"),
    ],
)
def test_bracket_tool_id_maps_to_expected_decorator(tool_id, expected_kind):
    """Cada opción visible de la paleta debe activar el decorador correcto."""
    canvas = ChemusonCanvas()

    canvas.set_current_tool(tool_id)

    assert canvas.current_tool == "tool_brackets"
    assert canvas.state.active_bracket_type == expected_kind


def test_corner_decorator_matches_single_top_right_angle():
    """Corner debe dibujar un solo ángulo superior derecho, no cuatro esquinas."""
    canvas = ChemusonCanvas()

    drag_rect = QRectF(150.0, 110.0, 70.0, 50.0)
    canvas.set_current_tool("tool_brackets_corner")
    canvas._bracket_drag_start = drag_rect.topLeft()
    canvas._last_scene_pos = drag_rect.bottomRight()

    canvas.mouseReleaseEvent(object())

    assert len(canvas.bracket_items) == 1
    item = canvas.bracket_items[0]
    path = item.path()
    assert path.elementCount() == 3
    start = path.elementAt(0)
    elbow = path.elementAt(1)
    end = path.elementAt(2)
    assert start.x == pytest.approx(drag_rect.left())
    assert start.y == pytest.approx(drag_rect.top())
    assert elbow.x == pytest.approx(drag_rect.right())
    assert elbow.y == pytest.approx(drag_rect.top())
    assert end.x == pytest.approx(drag_rect.right())
    assert end.y == pytest.approx(drag_rect.bottom())


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
