"""
Utilidades geométricas para enlaces en cuña.
"""
from __future__ import annotations

import math
from typing import Tuple

Point = Tuple[float, float]


def compute_wedge_points(
    p0: Point,
    p1: Point,
    width: float,
    *,
    trim_start: float = 0.0,
    trim_end: float = 0.0,
) -> tuple[Point, Point, Point]:
    """Calcula los puntos de una cuña (tip, base1, base2).

    Args:
        p0: Punto del vértice (punta de la cuña).
        p1: Punto del centro de la base.
        width: Ancho total de la base.
        trim_start: Recorte desde la punta (evita solapar etiquetas).
        trim_end: Recorte desde la base.

    Returns:
        Tupla `(tip, base1, base2)` con coordenadas 2D.
    """
    dx = p1[0] - p0[0]
    dy = p1[1] - p0[1]
    length = math.hypot(dx, dy)
    if length <= 1e-9:
        return p0, p0, p0
    ux = dx / length
    uy = dy / length
    tip_x = p0[0] + ux * trim_start
    tip_y = p0[1] + uy * trim_start
    base_x = p1[0] - ux * trim_end
    base_y = p1[1] - uy * trim_end
    nx = -uy
    ny = ux
    half_w = width / 2.0
    base1 = (base_x + nx * half_w, base_y + ny * half_w)
    base2 = (base_x - nx * half_w, base_y - ny * half_w)
    return (tip_x, tip_y), base1, base2
