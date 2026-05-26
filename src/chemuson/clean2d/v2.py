from __future__ import annotations

"""Motor puro de pulido geométrico 2D con restricciones químicas simples."""

from dataclasses import dataclass
import math
from typing import Iterable

from chemuson.core.model import Bond, MolGraph, bond_affects_valence


@dataclass(frozen=True)
class Clean2DParameters:
    """Parámetros configurables para el motor ``clean2d_v2``."""

    mode: str = "quick"
    iterations: int = 24
    target_bond_length: float | None = None
    length_weight: float = 0.34
    angle_weight: float = 0.018
    collision_weight: float = 0.18
    label_bond_weight: float = 0.10
    damping: float = 0.82
    max_step: float = 10.0

    @classmethod
    def quick(cls, target_bond_length: float | None = None) -> "Clean2DParameters":
        return cls(mode="quick", iterations=18, target_bond_length=target_bond_length)

    @classmethod
    def publication(cls, target_bond_length: float | None = None) -> "Clean2DParameters":
        return cls(
            mode="publication",
            iterations=72,
            target_bond_length=target_bond_length,
            length_weight=0.42,
            angle_weight=0.030,
            collision_weight=0.28,
            label_bond_weight=0.18,
            damping=0.78,
            max_step=7.0,
        )


def optimize_clean2d_positions(
    graph: MolGraph,
    atom_ids: Iterable[int] | None = None,
    params: Clean2DParameters | None = None,
) -> dict[int, tuple[float, float]]:
    """Devuelve coordenadas optimizadas sin mutar el grafo.

    El motor aplica cuatro restricciones ligeras:
    longitud objetivo por orden de enlace, ángulos preferidos por hibridación
    aproximada, separación átomo-átomo y separación átomo-enlace.
    """
    params = params or Clean2DParameters.quick()
    selected = set(atom_ids or graph.atoms.keys())
    selected = {atom_id for atom_id in selected if atom_id in graph.atoms}
    if not selected:
        return {}

    positions = {
        atom_id: (float(graph.atoms[atom_id].x), float(graph.atoms[atom_id].y))
        for atom_id in selected
    }
    bonds = [
        bond
        for bond in graph.bonds.values()
        if bond.a1_id in selected and bond.a2_id in selected and bond_affects_valence(bond)
    ]
    if not bonds:
        return positions

    target_len = float(params.target_bond_length or _average_bond_length(positions, bonds) or 42.0)
    target_len = max(8.0, target_len)
    center_before = _center(positions)

    velocities = {atom_id: (0.0, 0.0) for atom_id in selected}
    for _ in range(max(1, int(params.iterations))):
        forces = {atom_id: (0.0, 0.0) for atom_id in selected}
        _add_length_forces(forces, positions, bonds, target_len, params.length_weight)
        _add_angle_forces(forces, positions, graph, selected, bonds, params.angle_weight)
        _add_atom_collision_forces(forces, positions, target_len, params.collision_weight)
        _add_label_bond_forces(forces, positions, graph, selected, bonds, target_len, params.label_bond_weight)

        next_positions: dict[int, tuple[float, float]] = {}
        for atom_id, (x, y) in positions.items():
            fx, fy = forces[atom_id]
            vx, vy = velocities[atom_id]
            vx = (vx + fx) * params.damping
            vy = (vy + fy) * params.damping
            step = math.hypot(vx, vy)
            if step > params.max_step:
                scale = params.max_step / step
                vx *= scale
                vy *= scale
            velocities[atom_id] = (vx, vy)
            next_positions[atom_id] = (x + vx, y + vy)
        positions = next_positions

    center_after = _center(positions)
    shift_x = center_before[0] - center_after[0]
    shift_y = center_before[1] - center_after[1]
    return {
        atom_id: (x + shift_x, y + shift_y)
        for atom_id, (x, y) in positions.items()
    }


def _add_length_forces(
    forces: dict[int, tuple[float, float]],
    positions: dict[int, tuple[float, float]],
    bonds: list[Bond],
    target_len: float,
    weight: float,
) -> None:
    for bond in bonds:
        x1, y1 = positions[bond.a1_id]
        x2, y2 = positions[bond.a2_id]
        dx = x2 - x1
        dy = y2 - y1
        dist = math.hypot(dx, dy) or 1e-6
        desired = _target_length_for_bond(bond, target_len)
        delta = (dist - desired) * weight
        fx = dx / dist * delta
        fy = dy / dist * delta
        _accumulate(forces, bond.a1_id, fx, fy)
        _accumulate(forces, bond.a2_id, -fx, -fy)


def _add_angle_forces(
    forces: dict[int, tuple[float, float]],
    positions: dict[int, tuple[float, float]],
    graph: MolGraph,
    selected: set[int],
    bonds: list[Bond],
    weight: float,
) -> None:
    neighbors: dict[int, list[int]] = {}
    for bond in bonds:
        neighbors.setdefault(bond.a1_id, []).append(bond.a2_id)
        neighbors.setdefault(bond.a2_id, []).append(bond.a1_id)

    for center_id, neighs in neighbors.items():
        if len(neighs) < 2:
            continue
        cx, cy = positions[center_id]
        target = _target_angle_for_atom(graph, center_id, bonds)
        for idx, left_id in enumerate(neighs):
            for right_id in neighs[idx + 1 :]:
                if left_id not in selected or right_id not in selected:
                    continue
                lx, ly = positions[left_id]
                rx, ry = positions[right_id]
                a1 = math.atan2(ly - cy, lx - cx)
                a2 = math.atan2(ry - cy, rx - cx)
                diff = _angle_delta(a2 - a1)
                abs_diff = abs(diff)
                if abs_diff < 1e-5:
                    continue
                delta = abs_diff - target
                direction = 1.0 if diff >= 0.0 else -1.0
                torque = max(-0.12, min(0.12, delta * weight)) * direction
                _rotate_force(forces, positions, center_id, left_id, -torque)
                _rotate_force(forces, positions, center_id, right_id, torque)


def _add_atom_collision_forces(
    forces: dict[int, tuple[float, float]],
    positions: dict[int, tuple[float, float]],
    target_len: float,
    weight: float,
) -> None:
    ids = sorted(positions)
    threshold = target_len * 0.62
    for idx, a_id in enumerate(ids):
        ax, ay = positions[a_id]
        for b_id in ids[idx + 1 :]:
            bx, by = positions[b_id]
            dx = bx - ax
            dy = by - ay
            dist = math.hypot(dx, dy)
            if dist >= threshold:
                continue
            if dist <= 1e-6:
                dx = 1.0
                dy = 0.0
                dist = 1.0
            push = (threshold - dist) * weight
            fx = dx / dist * push
            fy = dy / dist * push
            _accumulate(forces, a_id, -fx, -fy)
            _accumulate(forces, b_id, fx, fy)


def _add_label_bond_forces(
    forces: dict[int, tuple[float, float]],
    positions: dict[int, tuple[float, float]],
    graph: MolGraph,
    selected: set[int],
    bonds: list[Bond],
    target_len: float,
    weight: float,
) -> None:
    threshold = target_len * 0.34
    for atom_id in selected:
        atom = graph.atoms.get(atom_id)
        if atom is None or atom.element == "C":
            continue
        ax, ay = positions[atom_id]
        for bond in bonds:
            if atom_id in {bond.a1_id, bond.a2_id}:
                continue
            x1, y1 = positions[bond.a1_id]
            x2, y2 = positions[bond.a2_id]
            dist, px, py = _point_segment_distance(ax, ay, x1, y1, x2, y2)
            if dist >= threshold:
                continue
            dx = ax - px
            dy = ay - py
            norm = math.hypot(dx, dy) or 1.0
            push = (threshold - dist) * weight
            _accumulate(forces, atom_id, dx / norm * push, dy / norm * push)


def _target_length_for_bond(bond: Bond, target_len: float) -> float:
    if getattr(bond, "is_aromatic", False):
        return target_len * 0.98
    order = int(getattr(bond, "order", 1) or 1)
    if order >= 3:
        return target_len * 0.92
    if order == 2:
        return target_len * 0.96
    return target_len


def _target_angle_for_atom(graph: MolGraph, atom_id: int, bonds: list[Bond]) -> float:
    has_triple = False
    has_double = False
    for bond in bonds:
        if atom_id not in {bond.a1_id, bond.a2_id}:
            continue
        has_triple = has_triple or int(getattr(bond, "order", 1) or 1) >= 3
        has_double = has_double or int(getattr(bond, "order", 1) or 1) == 2 or bool(getattr(bond, "is_aromatic", False))
    atom = graph.atoms.get(atom_id)
    if has_triple:
        return math.radians(180.0)
    if has_double or (atom is not None and bool(getattr(atom, "is_aromatic", False))):
        return math.radians(120.0)
    return math.radians(109.5)


def _rotate_force(
    forces: dict[int, tuple[float, float]],
    positions: dict[int, tuple[float, float]],
    center_id: int,
    atom_id: int,
    radians_step: float,
) -> None:
    cx, cy = positions[center_id]
    ax, ay = positions[atom_id]
    dx = ax - cx
    dy = ay - cy
    _accumulate(forces, atom_id, -dy * radians_step, dx * radians_step)


def _point_segment_distance(
    px: float,
    py: float,
    x1: float,
    y1: float,
    x2: float,
    y2: float,
) -> tuple[float, float, float]:
    dx = x2 - x1
    dy = y2 - y1
    denom = dx * dx + dy * dy
    if denom <= 1e-12:
        return math.hypot(px - x1, py - y1), x1, y1
    t = ((px - x1) * dx + (py - y1) * dy) / denom
    t = max(0.0, min(1.0, t))
    closest_x = x1 + t * dx
    closest_y = y1 + t * dy
    return math.hypot(px - closest_x, py - closest_y), closest_x, closest_y


def _average_bond_length(positions: dict[int, tuple[float, float]], bonds: list[Bond]) -> float:
    lengths = []
    for bond in bonds:
        x1, y1 = positions[bond.a1_id]
        x2, y2 = positions[bond.a2_id]
        lengths.append(math.hypot(x2 - x1, y2 - y1))
    return sum(lengths) / len(lengths) if lengths else 0.0


def _center(positions: dict[int, tuple[float, float]]) -> tuple[float, float]:
    if not positions:
        return 0.0, 0.0
    return (
        sum(x for x, _y in positions.values()) / len(positions),
        sum(y for _x, y in positions.values()) / len(positions),
    )


def _angle_delta(value: float) -> float:
    while value <= -math.pi:
        value += math.tau
    while value > math.pi:
        value -= math.tau
    return value


def _accumulate(forces: dict[int, tuple[float, float]], atom_id: int, fx: float, fy: float) -> None:
    old_x, old_y = forces[atom_id]
    forces[atom_id] = (old_x + fx, old_y + fy)
