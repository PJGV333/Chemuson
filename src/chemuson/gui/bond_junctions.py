"""Pure geometry helpers for molecular bond junction profiles."""
from __future__ import annotations

import math

Point = tuple[float, float]
Segment = tuple[Point, Point]


def hashed_wedge_segments(
    start: Point,
    end: Point,
    count: int,
    min_width: float,
    max_width: float,
) -> tuple[Segment, ...]:
    """Return uniformly spaced, monotonically growing hash strokes.

    Each mark is centered in one of ``count`` equal intervals along the bond,
    so neither end places a full transverse stroke directly on the atom vertex.
    Width changes linearly across those same interval centers. The function is
    presentation-only and does not depend on Qt or mutate its inputs.

    Invalid/non-finite geometry and non-positive counts produce no segments.
    Widths must be finite, non-negative, and ordered from narrow to wide.
    """
    try:
        steps = int(count)
        x0, y0 = float(start[0]), float(start[1])
        x1, y1 = float(end[0]), float(end[1])
        narrow = float(min_width)
        wide = float(max_width)
    except (TypeError, ValueError, IndexError, OverflowError):
        return ()

    if steps <= 0:
        return ()
    if not all(math.isfinite(value) for value in (x0, y0, x1, y1, narrow, wide)):
        return ()
    if narrow < 0.0 or wide < narrow:
        return ()

    dx = x1 - x0
    dy = y1 - y0
    length = math.hypot(dx, dy)
    if length <= 1e-9:
        return ()

    nx = -dy / length
    ny = dx / length
    segments: list[Segment] = []
    for index in range(steps):
        fraction = (index + 0.5) / steps
        cx = x0 + dx * fraction
        cy = y0 + dy * fraction
        width = narrow + (wide - narrow) * fraction
        half_width = width * 0.5
        segments.append(
            (
                (cx + nx * half_width, cy + ny * half_width),
                (cx - nx * half_width, cy - ny * half_width),
            )
        )
    return tuple(segments)
