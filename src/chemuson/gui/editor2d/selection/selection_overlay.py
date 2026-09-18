"""Deterministic geometry and hit policy for selection overlays.

Scene creation, item visibility and interaction coordination remain owned by
``CanvasSelectionMixin``. These helpers only calculate values or query handles.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable

from PyQt6.QtCore import QPoint, QPointF, QRectF


def padded_selection_bbox(bbox: QRectF, stroke_px: float) -> QRectF:
    """Return the current selection bbox with its visual padding."""
    padded = QRectF(bbox)
    pad = max(2.0, float(stroke_px))
    padded.adjust(-pad, -pad, pad, pad)
    return padded


def offset_scene_point(
    base: QPointF,
    dx_view: float,
    dy_view: float,
    *,
    map_from_scene: Callable[[QPointF], QPointF],
    map_to_scene: Callable[[QPoint], QPointF],
) -> QPointF:
    """Apply a rounded view-pixel offset and convert it back to scene space."""
    view_pt = map_from_scene(base)
    view_x = float(view_pt.x()) + dx_view
    view_y = float(view_pt.y()) + dy_view
    rounded_view_pt = QPoint(int(round(view_x)), int(round(view_y)))
    return map_to_scene(rounded_view_pt)


def selection_handle_scene_positions(
    padded: QRectF,
    *,
    offset_in_scene: Callable[[QPointF, float, float], QPointF],
    rotate_offset: float,
    move_offset: float,
    handle_radius: float,
) -> dict[str, QPointF]:
    """Calculate scene positions for rotate, move and scale handles."""
    top_center = QPointF(padded.center().x(), padded.top())
    return {
        "rotate": offset_in_scene(top_center, 0.0, -rotate_offset),
        "move": offset_in_scene(top_center, 0.0, move_offset),
        "scale": offset_in_scene(
            QPointF(padded.right(), padded.bottom()),
            -handle_radius,
            -handle_radius,
        ),
    }


def handle_item_distance_sq(
    handle: object,
    scene_pos: QPointF,
    *,
    map_from_scene: Callable[[QPointF], QPointF],
) -> float | None:
    """Return pointer-to-handle distance squared in view pixels."""
    view_pos = map_from_scene(scene_pos)
    try:
        center_scene = handle.mapToScene(handle.boundingRect().center())
    except RuntimeError:
        return None
    center_view = map_from_scene(center_scene)
    dx = float(view_pos.x() - center_view.x())
    dy = float(view_pos.y() - center_view.y())
    return dx * dx + dy * dy


def handle_item_hit_radius(
    handle: object,
    *,
    map_from_scene: Callable[[QPointF], QPointF],
    selection_handle_radius: float,
) -> float:
    """Return the effective view-pixel radius for a handle hit."""
    try:
        handle_rect_scene = handle.mapToScene(handle.boundingRect()).boundingRect()
        top_left_view = map_from_scene(handle_rect_scene.topLeft())
        bottom_right_view = map_from_scene(handle_rect_scene.bottomRight())
        visual_radius = max(
            abs(bottom_right_view.x() - top_left_view.x()),
            abs(bottom_right_view.y() - top_left_view.y()),
        ) * 0.75
    except Exception:
        visual_radius = 0.0
    return max(float(visual_radius), float(selection_handle_radius) * 3.0, 18.0)


def selection_handle_hit_kind(
    scene_pos: QPointF,
    handles: Iterable[tuple[str, object | None]],
    *,
    distance_sq: Callable[[object, QPointF], float | None],
    hit_radius: Callable[[object], float],
) -> str | None:
    """Return the closest visible handle kind, preserving iteration tie order."""
    candidates: list[tuple[float, str]] = []
    for kind, handle in handles:
        if handle is None:
            continue
        try:
            if not handle.isVisible():
                continue
        except RuntimeError:
            continue
        current_distance_sq = distance_sq(handle, scene_pos)
        if current_distance_sq is None:
            continue
        radius = hit_radius(handle)
        if current_distance_sq <= (radius * radius):
            candidates.append((current_distance_sq, kind))
    if not candidates:
        return None
    candidates.sort(key=lambda item: item[0])
    return candidates[0][1]
