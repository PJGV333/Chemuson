"""Functional contracts for selection overlay geometry and handle policy."""

from __future__ import annotations

from dataclasses import dataclass

import pytest
from PyQt6.QtCore import QPoint, QPointF, QRectF

from chemuson.gui.canvas.selection_overlay import (
    handle_item_distance_sq,
    handle_item_hit_radius,
    offset_scene_point,
    padded_selection_bbox,
    selection_handle_hit_kind,
    selection_handle_scene_positions,
)


@dataclass
class FakeRect:
    rect: QRectF

    def center(self):
        return self.rect.center()


class FakeHandle:
    def __init__(self, center: QPointF, visible: bool = True) -> None:
        self._center = center
        self._visible = visible

    def boundingRect(self):
        return QRectF(-2.0, -2.0, 4.0, 4.0)

    def mapToScene(self, point):
        return QPointF(self._center.x() + point.x(), self._center.y() + point.y())

    def isVisible(self):
        return self._visible


class ErrorHandle(FakeHandle):
    def mapToScene(self, _point):
        raise RuntimeError("deleted")


def test_padding_preserves_current_minimum_and_drag_translation() -> None:
    padded = padded_selection_bbox(QRectF(10.0, 20.0, 30.0, 40.0), 0.5)

    assert padded == QRectF(8.0, 18.0, 34.0, 44.0)
    assert padded.translated(QPointF(3.0, -4.0)) == QRectF(11.0, 14.0, 34.0, 44.0)


def test_offset_scene_point_uses_rounded_view_pixels() -> None:
    def map_from_scene(point):
        return QPointF(point.x() * 2.0, point.y() * 2.0)

    def map_to_scene(point):
        assert isinstance(point, QPoint)
        return QPointF(point.x() / 2.0, point.y() / 2.0)

    result = offset_scene_point(
        QPointF(10.0, 20.0),
        1.4,
        -2.6,
        map_from_scene=map_from_scene,
        map_to_scene=map_to_scene,
    )

    assert result == QPointF(10.5, 18.5)


def test_selection_handle_positions_preserve_three_handle_offsets() -> None:
    padded = QRectF(0.0, 0.0, 20.0, 30.0)
    positions = selection_handle_scene_positions(
        padded,
        offset_in_scene=lambda point, dx, dy: QPointF(point.x() + dx, point.y() + dy),
        rotate_offset=13.0,
        move_offset=0.0,
        handle_radius=6.0,
    )

    assert positions["rotate"] == QPointF(10.0, -13.0)
    assert positions["move"] == QPointF(10.0, 0.0)
    assert positions["scale"] == QPointF(14.0, 24.0)


def test_handle_distance_is_measured_in_screen_space() -> None:
    handle = FakeHandle(QPointF(10.0, 20.0))
    result = handle_item_distance_sq(
        handle,
        QPointF(13.0, 24.0),
        map_from_scene=lambda point: QPointF(point.x() * 2.0, point.y() * 3.0),
    )

    assert result == pytest.approx((6.0 * 6.0) + (12.0 * 12.0))


def test_handle_distance_returns_none_for_deleted_handle() -> None:
    assert handle_item_distance_sq(
        ErrorHandle(QPointF()),
        QPointF(),
        map_from_scene=lambda point: point,
    ) is None


def test_handle_radius_has_current_visual_and_minimum_policy() -> None:
    radius = handle_item_hit_radius(
        FakeHandle(QPointF()),
        map_from_scene=lambda point: QPointF(point.x() * 0.1, point.y() * 0.1),
        selection_handle_radius=6.0,
    )

    assert radius == 18.0


def test_handle_hit_kind_filters_none_invisible_and_picks_closest() -> None:
    handles = [
        ("scale", FakeHandle(QPointF(), visible=False)),
        ("rotate", FakeHandle(QPointF())),
        ("move", FakeHandle(QPointF())),
    ]
    distances = {id(handles[1][1]): 4.0, id(handles[2][1]): 9.0}

    assert selection_handle_hit_kind(
        QPointF(),
        handles,
        distance_sq=lambda handle, _point: distances[id(handle)],
        hit_radius=lambda _handle: 4.0,
    ) == "rotate"


def test_handle_hit_kind_keeps_current_order_on_equal_distance() -> None:
    handles = [("scale", FakeHandle(QPointF())), ("rotate", FakeHandle(QPointF()))]

    assert selection_handle_hit_kind(
        QPointF(),
        handles,
        distance_sq=lambda _handle, _point: 4.0,
        hit_radius=lambda _handle: 4.0,
    ) == "scale"
