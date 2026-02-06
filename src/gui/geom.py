"""
Utilidades geométricas para el motor de dibujo de Chemuson.

Incluye funciones puras para ángulos, snapping y selección de elementos.
"""
from __future__ import annotations

import math
from typing import Iterable, Optional, Tuple

from PyQt6.QtCore import QPointF


SP3_BOND_ANGLE_DEG = 109.5


def angle_deg(p0: QPointF, p1: QPointF) -> float:
    """Calcula el ángulo en grados (0-360) usando coordenadas matemáticas."""
    dx = p1.x() - p0.x()
    dy = -(p1.y() - p0.y())
    if dx == 0 and dy == 0:
        return 0.0
    return (math.degrees(math.atan2(dy, dx)) + 360.0) % 360.0


def snap_angle_deg(theta_deg: float, step_deg: float) -> float:
    """Ajusta un ángulo a un múltiplo del paso dado."""
    if step_deg <= 0:
        return theta_deg
    return (round(theta_deg / step_deg) * step_deg) % 360.0


def normalize_angle_deg(theta_deg: float) -> float:
    """Normaliza un ángulo al rango [0, 360)."""
    return theta_deg % 360.0


def angle_distance_deg(a_deg: float, b_deg: float) -> float:
    """Distancia angular mínima entre dos ángulos."""
    diff = (a_deg - b_deg + 180.0) % 360.0 - 180.0
    return abs(diff)


def nearly_colinear_deg(a_deg: float, b_deg: float, tolerance_deg: float = 20.0) -> bool:
    """Indica si dos direcciones son casi colineales (180° +/- tolerancia)."""
    return abs(angle_distance_deg(a_deg, b_deg) - 180.0) <= tolerance_deg


def endpoint_from_angle_len(p0: QPointF, theta_deg: float, length: float) -> QPointF:
    """Calcula el punto final desde un origen, ángulo y longitud."""
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
    """Infere la geometría (sp/sp2/sp3) según órdenes y vecinos.

    Args:
        bond_order: Orden del enlace principal.
        is_aromatic: Si el enlace es aromático.
        neighbor_angles_deg: Ángulos de enlaces vecinos.
        colinear_tolerance_deg: Tolerancia para detectar colinealidad.

    Returns:
        Cadena "sp", "sp2" o "sp3".
    """
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
    """Genera direcciones candidatas según la geometría.

    Args:
        geometry: Tipo de geometría ("sp", "sp2", "sp3").
        existing_angles_deg: Ángulos ya ocupados por enlaces existentes.
        base_angle_deg: Ángulo base sugerido por el cursor.

    Returns:
        Lista de ángulos candidatos (grados).
    """
    if geometry == "sp3":
        angles = list(existing_angles_deg)
        if angles:
            candidates = []
            for base in angles:
                candidates.append(normalize_angle_deg(base + SP3_BOND_ANGLE_DEG))
                candidates.append(normalize_angle_deg(base - SP3_BOND_ANGLE_DEG))
            seen = set()
            deduped = []
            for cand in candidates:
                key = round(cand, 6)
                if key in seen:
                    continue
                seen.add(key)
                deduped.append(cand)
            return deduped
        step = 60.0
        count = 6
    elif geometry == "sp2":
        step = 120.0
        count = 3
    else:
        step = 180.0
        count = 2

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
    """Filtra direcciones cercanas a ángulos ya ocupados."""
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
    """Elige el ángulo candidato más cercano al cursor.

    Args:
        candidates_deg: Ángulos candidatos.
        mouse_angle_deg: Ángulo del cursor.
        preferred_deg: Lista de ángulos preferidos (opcional).
        preferred_tolerance_deg: Tolerancia para preferir un ángulo.
        preferred_bonus: Bonificación (reducción) en la puntuación.

    Returns:
        Ángulo seleccionado o `None` si no hay candidatos.
    """
    candidates = list(candidates_deg)
    if not candidates:
        return None
    preferred = list(preferred_deg or [])

    def score(angle_deg: float) -> float:
        """Calcula la puntuación de un ángulo candidato."""
        score_value = angle_distance_deg(angle_deg, mouse_angle_deg)
        if preferred:
            if min(angle_distance_deg(angle_deg, pref) for pref in preferred) <= preferred_tolerance_deg:
                score_value -= preferred_bonus
        return score_value

    return min(candidates, key=score)


def choose_optimal_direction(angles_deg: Iterable[float]) -> float:
    """Devuelve el punto medio del mayor hueco angular."""
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
    """Distancia de un punto al segmento AB."""
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
    """Devuelve el ID del átomo más cercano dentro de un umbral."""
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
    """Devuelve el ID del enlace más cercano dentro de un umbral."""
    best_id = None
    best_dist = threshold
    for bond_id, a, b in bonds:
        dist = distance_point_to_segment(pos, a, b)
        if dist <= best_dist:
            best_id = bond_id
            best_dist = dist
    return best_id


def bond_side(a: QPointF, b: QPointF, cursor: QPointF) -> int:
    """Indica si el cursor está a la izquierda (1) o derecha (-1) del segmento."""
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
    """Determina si dos segmentos se intersectan.

    Args:
        a: Punto inicial del segmento AB.
        b: Punto final del segmento AB.
        c: Punto inicial del segmento CD.
        d: Punto final del segmento CD.
        exclude_endpoints: Si se excluyen intersecciones en extremos.
        eps: Tolerancia numérica.

    Returns:
        `True` si hay intersección según los criterios.
    """
    def orient(p: QPointF, q: QPointF, r: QPointF) -> float:
        """Orientación (producto cruzado) para test de intersección."""
        return (q.x() - p.x()) * (r.y() - p.y()) - (q.y() - p.y()) * (r.x() - p.x())

    def on_segment(p: QPointF, q: QPointF, r: QPointF) -> bool:
        """Comprueba si q está sobre el segmento pr (con tolerancia)."""
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
    """Distancia mínima entre dos segmentos."""
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
    """Indica si dos segmentos son casi iguales (dentro de tolerancia)."""
    def dist(p: QPointF, q: QPointF) -> float:
        """Distancia euclídea entre dos puntos."""
        return math.hypot(p.x() - q.x(), p.y() - q.y())

    return (dist(a, c) <= tolerance and dist(b, d) <= tolerance) or (
        dist(a, d) <= tolerance and dist(b, c) <= tolerance
    )
