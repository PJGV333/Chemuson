"""
Geometry helpers for Chemuson drawing engine.
Pure functions for snapping and picking.
"""
from __future__ import annotations

import math
from typing import Iterable, Optional, Tuple

from PyQt6.QtCore import QPointF


def angle_deg(p0: QPointF, p1: QPointF) -> float:
    """Angle in degrees (0-360) using math coordinates (Y up)."""
    dx = p1.x() - p0.x()
    dy = -(p1.y() - p0.y())
    if dx == 0 and dy == 0:
        return 0.0
    return (math.degrees(math.atan2(dy, dx)) + 360.0) % 360.0


def snap_angle_deg(theta_deg: float, step_deg: float) -> float:
    if step_deg <= 0:
        return theta_deg
    return (round(theta_deg / step_deg) * step_deg) % 360.0


def endpoint_from_angle_len(p0: QPointF, theta_deg: float, length: float) -> QPointF:
    rad = math.radians(theta_deg)
    dx = math.cos(rad) * length
    dy = -math.sin(rad) * length
    return QPointF(p0.x() + dx, p0.y() + dy)


def choose_optimal_direction(angles_deg: Iterable[float]) -> float:
    """Return the midpoint of the largest gap between angles."""
    angles = sorted(a % 360.0 for a in angles_deg)
    if not angles:
        return 0.0
    if len(angles) == 1:
        return (angles[0] + 180.0) % 360.0

    best_gap = -1.0
    best_angle = 0.0
    for i in range(len(angles)):
        a1 = angles[i]
        a2 = angles[(i + 1) % len(angles)]
        gap = (a2 - a1) % 360.0
        if gap > best_gap:
            best_gap = gap
            best_angle = (a1 + gap / 2.0) % 360.0
    return best_angle


def distance_point_to_segment(p: QPointF, a: QPointF, b: QPointF) -> float:
    """Distance from point to segment AB."""
    dx = b.x() - a.x()
    dy = b.y() - a.y()
    if dx == 0 and dy == 0:
        return math.hypot(p.x() - a.x(), p.y() - a.y())
    t = ((p.x() - a.x()) * dx + (p.y() - a.y()) * dy) / (dx * dx + dy * dy)
    t = max(0.0, min(1.0, t))
    proj = QPointF(a.x() + t * dx, a.y() + t * dy)
    return math.hypot(p.x() - proj.x(), p.y() - proj.y())


def closest_atom(
    pos: QPointF, atoms: Iterable[Tuple[int, float, float]], threshold: float
) -> Optional[int]:
    best_id = None
    best_dist = threshold
    for atom_id, x, y in atoms:
        dist = math.hypot(pos.x() - x, pos.y() - y)
        if dist <= best_dist:
            best_id = atom_id
            best_dist = dist
    return best_id


def closest_bond(
    pos: QPointF,
    bonds: Iterable[Tuple[int, QPointF, QPointF]],
    threshold: float,
) -> Optional[int]:
    best_id = None
    best_dist = threshold
    for bond_id, a, b in bonds:
        dist = distance_point_to_segment(pos, a, b)
        if dist <= best_dist:
            best_id = bond_id
            best_dist = dist
    return best_id


def bond_side(a: QPointF, b: QPointF, cursor: QPointF) -> int:
    """Return 1 for left of AB, -1 for right."""
    dx1 = b.x() - a.x()
    dy1 = b.y() - a.y()
    dx2 = cursor.x() - a.x()
    dy2 = cursor.y() - a.y()
    cross = dx1 * dy2 - dy1 * dx2
    return 1 if cross >= 0 else -1
