from __future__ import annotations

import math
from typing import Any

from chemuson.clean2d.safety import has_cycles


def _bond_a1_helper(bond: Any) -> int:
    return bond.a1_id if hasattr(bond, "a1_id") else bond.get("a1_id", 0)

def _bond_a2_helper(bond: Any) -> int:
    return bond.a2_id if hasattr(bond, "a2_id") else bond.get("a2_id", 0)

def length_only_polish(
    positions: dict[int, tuple[float, float]],
    bonds: list[Any],
    target_len: float,
    *,
    max_iterations: int = 12,
    max_displacement_per_atom: float | None = None,
    damping: float = 0.6,
) -> dict[int, tuple[float, float]]:
    """Ajusta suavemente longitudes de enlace sin cambiar topología global.

    Esta función NO reconstruye componentes, NO aplica fuerzas de ángulo
    ni colisiones, y NO modifica la forma interna de anillos.

    Args:
        positions: Coordenadas actuales por atom_id.
        bonds: Lista de enlaces (objetos con a1_id, a2_id, order).
        target_len: Longitud de enlace objetivo en píxeles.
        max_iterations: Iteraciones máximas de relajación.
        max_displacement_per_atom: Límite de desplazamiento por átomo.
        damping: Factor de amortiguación (0-1).

    Returns:
        Diccionario de coordenadas ajustadas. Solo los átomos que se
        movieron significativamente aparecen (los demás mantienen su
        posición original si se consulta el dict original).
    """
    if not positions or not bonds:
        return {}

    target_len = max(8.0, float(target_len))
    max_disp = (
        float(max_displacement_per_atom)
        if max_displacement_per_atom is not None
        else target_len * 2.0
    )

    atom_ids = set(positions.keys())
    cyclic = has_cycles(atom_ids, bonds)

    out: dict[int, tuple[float, float]] = {}
    for aid, (x, y) in positions.items():
        out[aid] = (float(x), float(y))

    for _iteration in range(max(1, int(max_iterations))):
        max_move = 0.0
        for bond in bonds:
            a1 = _bond_a1_helper(bond)
            a2 = _bond_a2_helper(bond)
            if a1 not in out or a2 not in out:
                continue

            x1, y1 = out[a1]
            x2, y2 = out[a2]
            dx = x2 - x1
            dy = y2 - y1
            dist = math.hypot(dx, dy)
            if dist < 1e-6:
                continue

            desired = _target_length_for_order(bond, target_len)
            delta = (dist - desired) * damping * 0.5

            fx = dx / dist * delta
            fy = dy / dist * delta

            # Para moléculas con ciclos, reducimos el movimiento a la mitad
            scale = 0.5 if cyclic else 1.0

            out[a1] = (out[a1][0] + fx * scale, out[a1][1] + fy * scale)
            out[a2] = (out[a2][0] - fx * scale, out[a2][1] - fy * scale)

            move = abs(fx) + abs(fy)
            if move > max_move:
                max_move = move

        if max_move < 0.5:
            break

    changed: dict[int, tuple[float, float]] = {}
    for aid, (x, y) in out.items():
        ox, oy = positions.get(aid, (x, y))
        d = math.hypot(x - ox, y - oy)
        if d > 0.5:
            clamped_x = ox + max(-max_disp, min(max_disp, x - ox))
            clamped_y = oy + max(-max_disp, min(max_disp, y - oy))
            changed[aid] = (clamped_x, clamped_y)

    return changed


def _target_length_for_order(bond: Any, target_len: float) -> float:
    is_aromatic = bool(bond.is_aromatic) if hasattr(bond, "is_aromatic") else bool(bond.get("is_aromatic", False))
    if is_aromatic:
        return target_len * 0.98
    order = int(bond.order if hasattr(bond, "order") else bond.get("order", 1) or 1)
    if order >= 3:
        return target_len * 0.92
    if order == 2:
        return target_len * 0.96
    return target_len
