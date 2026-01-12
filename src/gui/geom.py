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


def normalize_angle_deg(theta_deg: float) -> float:
    return theta_deg % 360.0


def angle_distance_deg(a_deg: float, b_deg: float) -> float:
    diff = (a_deg - b_deg + 180.0) % 360.0 - 180.0
    return abs(diff)


def nearly_colinear_deg(a_deg: float, b_deg: float, tolerance_deg: float = 20.0) -> bool:
    return abs(angle_distance_deg(a_deg, b_deg) - 180.0) <= tolerance_deg


def endpoint_from_angle_len(p0: QPointF, theta_deg: float, length: float) -> QPointF:
    rad = math.radians(theta_deg)
    dx = math.cos(rad) * length
    dy = -math.sin(rad) * length
    return QPointF(p0.x() + dx, p0.y() + dy)


def geometry_for_bond(
    bond_order: int,
    is_aromatic: bool,
    neighbor_angles_deg: Iterable[float],
    colinear_tolerance_deg: float = 20.0,
) -> str:
    if bond_order >= 3:
        base = "sp"
    elif bond_order == 2 or is_aromatic:
        base = "sp2"
    else:
        base = "sp3"

    angles = list(neighbor_angles_deg)
    if len(angles) >= 4:
        return "sp3"
    if len(angles) == 3:
        return "sp2"
    if len(angles) == 2 and nearly_colinear_deg(angles[0], angles[1], colinear_tolerance_deg):
        return "sp"
    return base


def candidate_directions_deg(
    geometry: str,
    existing_angles_deg: Iterable[float],
    base_angle_deg: float,
) -> list[float]:
    if geometry == "sp":
        step = 180.0
        count = 2
    elif geometry == "sp2":
        step = 120.0
        count = 3
    else:
        step = 60.0
        count = 6

    angles = list(existing_angles_deg)
    if angles:
        base = angles[0]
    else:
        base = snap_angle_deg(base_angle_deg, step)

    candidates = [normalize_angle_deg(base + step * i) for i in range(count)]
    return candidates


def filter_occupied_angles_deg(
    candidates_deg: Iterable[float],
    occupied_deg: Iterable[float],
    tolerance_deg: float = 20.0,
) -> list[float]:
    occupied = list(occupied_deg)
    if not occupied:
        return list(candidates_deg)
    filtered = []
    for cand in candidates_deg:
        if all(angle_distance_deg(cand, occ) > tolerance_deg for occ in occupied):
            filtered.append(cand)
    return filtered


def pick_closest_direction_deg(
    candidates_deg: Iterable[float],
    mouse_angle_deg: float,
    preferred_deg: Optional[Iterable[float]] = None,
    preferred_tolerance_deg: float = 15.0,
    preferred_bonus: float = 20.0,
) -> Optional[float]:
    candidates = list(candidates_deg)
    if not candidates:
        return None
    preferred = list(preferred_deg or [])

    def score(angle_deg: float) -> float:
        score_value = angle_distance_deg(angle_deg, mouse_angle_deg)
        if preferred:
            if min(angle_distance_deg(angle_deg, pref) for pref in preferred) <= preferred_tolerance_deg:
                score_value -= preferred_bonus
        return score_value

    return min(candidates, key=score)


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


def segments_intersect(
    a: QPointF,
    b: QPointF,
    c: QPointF,
    d: QPointF,
    exclude_endpoints: bool = True,
    eps: float = 1e-6,
) -> bool:
    def orient(p: QPointF, q: QPointF, r: QPointF) -> float:
        return (q.x() - p.x()) * (r.y() - p.y()) - (q.y() - p.y()) * (r.x() - p.x())

    def on_segment(p: QPointF, q: QPointF, r: QPointF) -> bool:
        return (
            min(p.x(), r.x()) - eps <= q.x() <= max(p.x(), r.x()) + eps
            and min(p.y(), r.y()) - eps <= q.y() <= max(p.y(), r.y()) + eps
        )

    o1 = orient(a, b, c)
    o2 = orient(a, b, d)
    o3 = orient(c, d, a)
    o4 = orient(c, d, b)

    if o1 == 0 and on_segment(a, c, b):
        return not exclude_endpoints or (c != a and c != b)
    if o2 == 0 and on_segment(a, d, b):
        return not exclude_endpoints or (d != a and d != b)
    if o3 == 0 and on_segment(c, a, d):
        return not exclude_endpoints or (a != c and a != d)
    if o4 == 0 and on_segment(c, b, d):
        return not exclude_endpoints or (b != c and b != d)

    return (o1 > 0) != (o2 > 0) and (o3 > 0) != (o4 > 0)


def segment_min_distance(a: QPointF, b: QPointF, c: QPointF, d: QPointF) -> float:
    return min(
        distance_point_to_segment(a, c, d),
        distance_point_to_segment(b, c, d),
        distance_point_to_segment(c, a, b),
        distance_point_to_segment(d, a, b),
    )


def segments_nearly_equal(
    a: QPointF,
    b: QPointF,
    c: QPointF,
    d: QPointF,
    tolerance: float = 5.0,
) -> bool:
    def dist(p: QPointF, q: QPointF) -> float:
        return math.hypot(p.x() - q.x(), p.y() - q.y())

    return (dist(a, c) <= tolerance and dist(b, d) <= tolerance) or (
        dist(a, d) <= tolerance and dist(b, c) <= tolerance
    )
