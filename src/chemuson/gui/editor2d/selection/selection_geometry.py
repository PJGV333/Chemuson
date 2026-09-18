"""Deterministic value geometry used by canvas selection transforms."""

from __future__ import annotations

import math

from PyQt6.QtCore import QPointF


def signed_angle_delta_deg(start_deg: float, end_deg: float) -> float:
    """Return the minimum signed angular delta ``end - start``."""
    return (float(end_deg) - float(start_deg) + 180.0) % 360.0 - 180.0


def rotate_scene_point(point: QPointF, center: QPointF, delta_deg: float) -> QPointF:
    """Rotate a scene point around a center."""
    rad = math.radians(float(delta_deg))
    dx = point.x() - center.x()
    dy = point.y() - center.y()
    cos_t = math.cos(rad)
    sin_t = math.sin(rad)
    return QPointF(
        center.x() + dx * cos_t - dy * sin_t,
        center.y() + dx * sin_t + dy * cos_t,
    )


def optional_float_equal(
    a: float | None,
    b: float | None,
    tol: float = 0.05,
) -> bool:
    """Compare two optional floats using the visual tolerance."""
    if a is None and b is None:
        return True
    if a is None or b is None:
        return False
    return abs(float(a) - float(b)) <= tol


def point_equal(a: QPointF, b: QPointF, tol: float = 1e-4) -> bool:
    """Compare two points componentwise using a tolerance."""
    return abs(a.x() - b.x()) <= tol and abs(a.y() - b.y()) <= tol


def scale_point_from_anchor(
    anchor: QPointF,
    point: QPointF,
    scale: float,
) -> QPointF:
    """Scale a point around an anchor."""
    return QPointF(
        anchor.x() + (point.x() - anchor.x()) * scale,
        anchor.y() + (point.y() - anchor.y()) * scale,
    )


def normalize_label_scale(value: float) -> float | None:
    """Normalize a local label scale for inheritance."""
    scale = max(0.2, float(value))
    return None if abs(scale - 1.0) < 0.02 else scale


def normalize_custom_stroke(
    value: float,
    inherited_stroke_px: float,
) -> float | None:
    """Normalize an effective stroke width to an override or inheritance."""
    stroke = max(0.6, float(value))
    return None if abs(stroke - float(inherited_stroke_px)) < 0.05 else stroke
