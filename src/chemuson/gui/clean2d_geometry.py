"""Geometría pura para integrar coordenadas Clean2D en la GUI."""

from __future__ import annotations

import math


def coords_center(
    coords: dict[int, tuple[float, float]],
) -> tuple[float, float]:
    """Calcula el centro geométrico de un conjunto de coordenadas."""
    xs = [x for x, _ in coords.values()]
    ys = [y for _, y in coords.values()]
    if not xs:
        return 0.0, 0.0
    return sum(xs) / len(xs), sum(ys) / len(ys)


def average_bond_length(
    coords: dict[int, tuple[float, float]],
    bonds,
) -> float:
    """Promedio de longitudes de enlaces disponibles en ``coords``."""
    lengths: list[float] = []
    for bond in bonds:
        p1 = coords.get(bond.a1_id)
        p2 = coords.get(bond.a2_id)
        if p1 is None or p2 is None:
            continue
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        dist = math.hypot(dx, dy)
        if dist > 1e-6:
            lengths.append(dist)
    if not lengths:
        return 0.0
    return sum(lengths) / len(lengths)


def rescale_coords_to_bond_length(
    coords: dict[int, tuple[float, float]],
    bonds,
    target_length: float,
) -> dict[int, tuple[float, float]]:
    """Reescala coordenadas alrededor de su centro para ajustar longitud media."""
    if not coords:
        return {}
    target = float(target_length)
    if target <= 1e-6:
        return dict(coords)
    current = average_bond_length(coords, bonds)
    if current <= 1e-6:
        return dict(coords)
    scale = target / current
    if not math.isfinite(scale):
        return dict(coords)
    # Evita explosiones numéricas por geometrías degeneradas.
    scale = max(0.05, min(200.0, scale))
    cx, cy = coords_center(coords)
    return {
        atom_id: (cx + (x - cx) * scale, cy + (y - cy) * scale)
        for atom_id, (x, y) in coords.items()
    }


def align_coords_to_reference(
    reference: dict[int, tuple[float, float]],
    coords: dict[int, tuple[float, float]],
) -> dict[int, tuple[float, float]]:
    """Alinea ``coords`` a la pose actual usando rotación rígida o reflexión."""
    if not reference or not coords:
        return dict(coords)

    common_ids = [atom_id for atom_id in reference if atom_id in coords]
    if not common_ids:
        return dict(coords)

    ref_common = {atom_id: reference[atom_id] for atom_id in common_ids}
    coord_common = {atom_id: coords[atom_id] for atom_id in common_ids}
    ref_cx, ref_cy = coords_center(ref_common)
    src_cx, src_cy = coords_center(coord_common)

    if len(common_ids) == 1:
        dx = ref_cx - src_cx
        dy = ref_cy - src_cy
        return {
            atom_id: (x + dx, y + dy)
            for atom_id, (x, y) in coords.items()
        }

    def _candidate(
        mirror_x: bool,
    ) -> tuple[dict[int, tuple[float, float]], float]:
        sum_cos = 0.0
        sum_sin = 0.0
        for atom_id in common_ids:
            qx, qy = coords[atom_id]
            px, py = reference[atom_id]
            qx -= src_cx
            qy -= src_cy
            px -= ref_cx
            py -= ref_cy
            if mirror_x:
                qx = -qx
            sum_cos += qx * px + qy * py
            sum_sin += qx * py - qy * px

        theta = math.atan2(sum_sin, sum_cos)
        cos_t = math.cos(theta)
        sin_t = math.sin(theta)

        transformed: dict[int, tuple[float, float]] = {}
        error = 0.0
        for atom_id, (x, y) in coords.items():
            qx = x - src_cx
            qy = y - src_cy
            if mirror_x:
                qx = -qx
            rx = qx * cos_t - qy * sin_t + ref_cx
            ry = qx * sin_t + qy * cos_t + ref_cy
            transformed[atom_id] = (rx, ry)
            if atom_id in reference:
                px, py = reference[atom_id]
                error += (rx - px) ** 2 + (ry - py) ** 2
        return transformed, error

    direct, direct_error = _candidate(False)
    mirrored, mirrored_error = _candidate(True)
    return mirrored if mirrored_error < direct_error else direct


def project_missing_hydrogen_coords(
    before: dict[int, tuple[float, float]],
    after: dict[int, tuple[float, float]],
    bonds,
    atom_elements: dict[int, str],
) -> dict[int, tuple[float, float]]:
    """Proyecta H faltantes trasladando su vector local desde átomos ancla."""
    projected = dict(after)
    if not before:
        return projected

    adjacency: dict[int, list[int]] = {}
    for bond in bonds:
        adjacency.setdefault(bond.a1_id, []).append(bond.a2_id)
        adjacency.setdefault(bond.a2_id, []).append(bond.a1_id)

    for atom_id, element in atom_elements.items():
        if element != "H" or atom_id in projected:
            continue
        before_pos = before.get(atom_id)
        if before_pos is None:
            continue
        candidate_positions: list[tuple[float, float]] = []
        for anchor_id in adjacency.get(atom_id, []):
            anchor_before = before.get(anchor_id)
            anchor_after = projected.get(anchor_id)
            if anchor_before is None or anchor_after is None:
                continue
            dx = before_pos[0] - anchor_before[0]
            dy = before_pos[1] - anchor_before[1]
            candidate_positions.append(
                (anchor_after[0] + dx, anchor_after[1] + dy)
            )
        if candidate_positions:
            avg_x = sum(pos[0] for pos in candidate_positions) / len(
                candidate_positions
            )
            avg_y = sum(pos[1] for pos in candidate_positions) / len(
                candidate_positions
            )
            projected[atom_id] = (avg_x, avg_y)

    for atom_id, coord in before.items():
        projected.setdefault(atom_id, coord)
    return projected
