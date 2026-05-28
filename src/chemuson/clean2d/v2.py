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
    max_component_scale: float = 3.0

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
            max_component_scale=4.0,
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

    positions = _normalize_component_lengths(
        positions,
        bonds,
        target_len,
        max_component_scale=max(1.0, float(params.max_component_scale)),
    )
    positions = _rebuild_acyclic_components_if_needed(positions, graph, selected, bonds, target_len)
    _nudge_degenerate_angles(positions, graph, selected, bonds, target_len)
    for _ in range(max(1, int(params.iterations))):
        before_step = dict(positions)
        forces = {atom_id: (0.0, 0.0) for atom_id in positions}
        _add_length_forces(forces, positions, bonds, target_len, params.length_weight)
        _add_angle_forces(forces, positions, graph, selected, bonds, params.angle_weight)
        _apply_forces(positions, forces, damping=params.damping, max_step=params.max_step)
        _resolve_nonbonded_atom_collisions(positions, bonds, target_len, params.collision_weight, params.max_step)
        _resolve_label_bond_collisions(positions, graph, selected, bonds, target_len, params.label_bond_weight, params.max_step)
        if _max_position_delta(before_step, positions) < 0.01:
            break

    center_after = _center(positions)
    shift_x = center_before[0] - center_after[0]
    shift_y = center_before[1] - center_after[1]
    return {
        atom_id: (x + shift_x, y + shift_y)
        for atom_id, (x, y) in positions.items()
    }


def _normalize_component_lengths(
    positions: dict[int, tuple[float, float]],
    bonds: list[Bond],
    target_len: float,
    *,
    max_component_scale: float,
) -> dict[int, tuple[float, float]]:
    """Escala cada componente de forma uniforme, preservando ángulos relativos."""
    adjacency: dict[int, list[int]] = {atom_id: [] for atom_id in positions}
    for bond in bonds:
        adjacency.setdefault(bond.a1_id, []).append(bond.a2_id)
        adjacency.setdefault(bond.a2_id, []).append(bond.a1_id)

    out = dict(positions)
    visited: set[int] = set()
    for start in sorted(positions):
        if start in visited:
            continue
        stack = [start]
        component: set[int] = set()
        while stack:
            atom_id = stack.pop()
            if atom_id in component:
                continue
            component.add(atom_id)
            visited.add(atom_id)
            stack.extend(adjacency.get(atom_id, []))
        comp_bonds = [bond for bond in bonds if bond.a1_id in component and bond.a2_id in component]
        if not comp_bonds:
            continue
        current = _average_bond_length(out, comp_bonds)
        if current <= 1e-6:
            continue
        scale = target_len / current
        scale = max(1.0 / max_component_scale, min(max_component_scale, scale))
        if abs(scale - 1.0) <= 1e-6:
            continue
        cx, cy = _center({atom_id: out[atom_id] for atom_id in component})
        for atom_id in component:
            x, y = out[atom_id]
            out[atom_id] = (cx + (x - cx) * scale, cy + (y - cy) * scale)
    return out


def _resolve_nonbonded_atom_collisions(
    positions: dict[int, tuple[float, float]],
    bonds: list[Bond],
    target_len: float,
    weight: float,
    max_step: float,
) -> None:
    bonded = {tuple(sorted((bond.a1_id, bond.a2_id))) for bond in bonds}
    ids = sorted(positions)
    threshold = target_len * 0.55
    for idx, a_id in enumerate(ids):
        ax, ay = positions[a_id]
        for b_id in ids[idx + 1 :]:
            if tuple(sorted((a_id, b_id))) in bonded:
                continue
            bx, by = positions[b_id]
            dx = bx - ax
            dy = by - ay
            dist = math.hypot(dx, dy)
            if dist >= threshold:
                continue
            if dist <= 1e-6:
                dx, dy, dist = _deterministic_unit_vector(a_id, b_id)
            push = min(max_step, (threshold - dist) * max(0.0, weight))
            fx = dx / dist * push * 0.5
            fy = dy / dist * push * 0.5
            positions[a_id] = (ax - fx, ay - fy)
            positions[b_id] = (bx + fx, by + fy)


def _resolve_label_bond_collisions(
    positions: dict[int, tuple[float, float]],
    graph: MolGraph,
    selected: set[int],
    bonds: list[Bond],
    target_len: float,
    weight: float,
    max_step: float,
) -> None:
    threshold = target_len * 0.32
    for atom_id in sorted(selected):
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
            norm = math.hypot(dx, dy)
            if norm <= 1e-6:
                sx = x2 - x1
                sy = y2 - y1
                norm = math.hypot(sx, sy) or 1.0
                dx = -sy / norm
                dy = sx / norm
                norm = 1.0
            push = min(max_step, (threshold - dist) * max(0.0, weight))
            positions[atom_id] = (ax + dx / norm * push, ay + dy / norm * push)


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


def _rebuild_acyclic_components_if_needed(
    positions: dict[int, tuple[float, float]],
    graph: MolGraph,
    selected: set[int],
    bonds: list[Bond],
    target_len: float,
) -> dict[int, tuple[float, float]]:
    adjacency = _adjacency_for_bonds(bonds)
    bond_lookup = {
        tuple(sorted((bond.a1_id, bond.a2_id))): bond
        for bond in bonds
    }
    out = dict(positions)
    visited: set[int] = set()
    for start in sorted(selected):
        if start in visited:
            continue
        stack = [start]
        component: set[int] = set()
        while stack:
            atom_id = stack.pop()
            if atom_id in component:
                continue
            component.add(atom_id)
            visited.add(atom_id)
            stack.extend(neigh for neigh in adjacency.get(atom_id, []) if neigh in selected)
        comp_bonds = [
            bond for bond in bonds
            if bond.a1_id in component and bond.a2_id in component
        ]
        if not comp_bonds or _component_has_cycle(component, adjacency):
            continue
        if not _acyclic_component_needs_rebuild(out, graph, component, comp_bonds, target_len):
            continue
        before_center = _center({atom_id: out[atom_id] for atom_id in component})
        rebuilt = _layout_acyclic_component(out, graph, component, adjacency, bond_lookup, comp_bonds, target_len)
        after_center = _center(rebuilt)
        shift_x = before_center[0] - after_center[0]
        shift_y = before_center[1] - after_center[1]
        for atom_id, (x, y) in rebuilt.items():
            out[atom_id] = (x + shift_x, y + shift_y)
    return out


def _adjacency_for_bonds(bonds: list[Bond]) -> dict[int, list[int]]:
    adjacency: dict[int, list[int]] = {}
    for bond in bonds:
        adjacency.setdefault(bond.a1_id, []).append(bond.a2_id)
        adjacency.setdefault(bond.a2_id, []).append(bond.a1_id)
    return adjacency


def _component_has_cycle(component: set[int], adjacency: dict[int, list[int]]) -> bool:
    visited: set[int] = set()
    for start in component:
        if start in visited:
            continue
        stack: list[tuple[int, int]] = [(start, -1)]
        while stack:
            atom_id, parent_id = stack.pop()
            if atom_id in visited:
                continue
            visited.add(atom_id)
            for neigh in adjacency.get(atom_id, []):
                if neigh not in component or neigh == parent_id:
                    continue
                if neigh in visited:
                    return True
                stack.append((neigh, atom_id))
    return False


def _acyclic_component_needs_rebuild(
    positions: dict[int, tuple[float, float]],
    graph: MolGraph,
    component: set[int],
    bonds: list[Bond],
    target_len: float,
) -> bool:
    for bond in bonds:
        x1, y1 = positions[bond.a1_id]
        x2, y2 = positions[bond.a2_id]
        dist = math.hypot(x2 - x1, y2 - y1)
        desired = _target_length_for_bond(bond, target_len)
        if desired > 1e-6 and (dist < desired * 0.70 or dist > desired * 1.30):
            return True

    adjacency = _adjacency_for_bonds(bonds)
    for center_id in component:
        neighs = [neigh for neigh in adjacency.get(center_id, []) if neigh in component]
        if len(neighs) != 2:
            continue
        target = _target_angle_for_atom(graph, center_id, bonds)
        if target >= math.radians(165.0):
            continue
        left_id, right_id = neighs
        cx, cy = positions[center_id]
        lx, ly = positions[left_id]
        rx, ry = positions[right_id]
        a1 = math.atan2(ly - cy, lx - cx)
        a2 = math.atan2(ry - cy, rx - cx)
        current = abs(_angle_delta(a2 - a1))
        if current <= math.radians(12.0) or abs(math.pi - current) <= math.radians(12.0):
            return True
    return False


def _layout_acyclic_component(
    positions: dict[int, tuple[float, float]],
    graph: MolGraph,
    component: set[int],
    adjacency: dict[int, list[int]],
    bond_lookup: dict[tuple[int, int], Bond],
    bonds: list[Bond],
    target_len: float,
) -> dict[int, tuple[float, float]]:
    root = _choose_tree_root(positions, graph, component, adjacency)
    laid_out = {root: positions[root]}
    parent_of: dict[int, int | None] = {root: None}
    queue = [root]
    while queue:
        parent_id = queue.pop(0)
        px, py = laid_out[parent_id]
        existing_angles: list[float] = []
        parent_parent = parent_of.get(parent_id)
        if parent_parent is not None and parent_parent in laid_out:
            gx, gy = laid_out[parent_parent]
            existing_angles.append(math.atan2(gy - py, gx - px))

        neighbors = [
            neigh for neigh in adjacency.get(parent_id, [])
            if neigh in component and neigh not in laid_out
        ]
        neighbors.sort(key=lambda neigh: _neighbor_layout_priority(graph, neigh, component, adjacency))
        for neigh in neighbors:
            ox, oy = positions.get(neigh, (px + 1.0, py))
            old_angle = math.atan2(oy - positions[parent_id][1], ox - positions[parent_id][0])
            candidates = _layout_angle_candidates(graph, parent_id, bonds, existing_angles, old_angle)
            chosen = min(candidates, key=lambda angle: abs(_angle_delta(angle - old_angle)))
            bond = bond_lookup.get(tuple(sorted((parent_id, neigh))))
            length = _target_length_for_bond(bond, target_len) if bond is not None else target_len
            laid_out[neigh] = (
                px + math.cos(chosen) * length,
                py + math.sin(chosen) * length,
            )
            existing_angles.append(chosen)
            parent_of[neigh] = parent_id
            queue.append(neigh)
    return laid_out


def _choose_tree_root(
    positions: dict[int, tuple[float, float]],
    graph: MolGraph,
    component: set[int],
    adjacency: dict[int, list[int]],
) -> int:
    heavy = [
        atom_id for atom_id in component
        if graph.atoms.get(atom_id) is None or graph.atoms[atom_id].element != "H"
    ]
    candidates = heavy if heavy else list(component)
    leaves = [
        atom_id for atom_id in candidates
        if _non_h_degree(graph, atom_id, component, adjacency) <= 1
    ]
    if leaves:
        candidates = leaves
    return min(candidates, key=lambda atom_id: (positions[atom_id][0], positions[atom_id][1], atom_id))


def _neighbor_layout_priority(
    graph: MolGraph,
    atom_id: int,
    component: set[int],
    adjacency: dict[int, list[int]],
) -> tuple[int, int, int]:
    atom = graph.atoms.get(atom_id)
    is_h = 1 if atom is not None and atom.element == "H" else 0
    return is_h, -_non_h_degree(graph, atom_id, component, adjacency), atom_id


def _non_h_degree(
    graph: MolGraph,
    atom_id: int,
    component: set[int],
    adjacency: dict[int, list[int]],
) -> int:
    degree = 0
    for neigh in adjacency.get(atom_id, []):
        if neigh not in component:
            continue
        atom = graph.atoms.get(neigh)
        if atom is not None and atom.element == "H":
            continue
        degree += 1
    return degree


def _layout_angle_candidates(
    graph: MolGraph,
    atom_id: int,
    bonds: list[Bond],
    existing_angles: list[float],
    old_angle: float,
) -> list[float]:
    if not existing_angles:
        return [old_angle]

    target = _target_angle_for_atom(graph, atom_id, bonds)
    raw: list[float] = []
    if target >= math.radians(165.0):
        raw = [angle + math.pi for angle in existing_angles]
    else:
        for angle in existing_angles:
            raw.append(angle + target)
            raw.append(angle - target)

    filtered = [
        angle for angle in raw
        if all(abs(_angle_delta(angle - used)) >= math.radians(24.0) for used in existing_angles)
    ]
    return filtered or raw or [old_angle]


def _nudge_degenerate_angles(
    positions: dict[int, tuple[float, float]],
    graph: MolGraph,
    selected: set[int],
    bonds: list[Bond],
    target_len: float,
) -> None:
    neighbors: dict[int, list[int]] = {}
    for bond in bonds:
        neighbors.setdefault(bond.a1_id, []).append(bond.a2_id)
        neighbors.setdefault(bond.a2_id, []).append(bond.a1_id)

    for center_id in sorted(selected):
        neighs = sorted(neighbors.get(center_id, []))
        if len(neighs) != 2:
            continue
        target = _target_angle_for_atom(graph, center_id, bonds)
        if target >= math.radians(165.0):
            continue
        left_id, right_id = neighs
        cx, cy = positions[center_id]
        lx, ly = positions[left_id]
        rx, ry = positions[right_id]
        a1 = math.atan2(ly - cy, lx - cx)
        a2 = math.atan2(ry - cy, rx - cx)
        current = abs(_angle_delta(a2 - a1))
        is_degenerate = current <= math.radians(8.0) or abs(math.pi - current) <= math.radians(8.0)
        if not is_degenerate:
            continue
        axis_x = rx - lx
        axis_y = ry - ly
        norm = math.hypot(axis_x, axis_y)
        if norm <= 1e-6:
            axis_x, axis_y, norm = _deterministic_unit_vector(left_id, right_id)
        sign = -1.0 if int(center_id) % 2 else 1.0
        step = max(1.0, min(target_len * 0.16, target_len - 1.0))
        positions[center_id] = (
            cx + (-axis_y / norm) * step * sign,
            cy + (axis_x / norm) * step * sign,
        )


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
        tolerance = _angle_tolerance_for_atom(graph, center_id, bonds)
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
                if abs(abs_diff - target) <= tolerance:
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


def _angle_tolerance_for_atom(graph: MolGraph, atom_id: int, bonds: list[Bond]) -> float:
    has_triple = False
    has_double = False
    for bond in bonds:
        if atom_id not in {bond.a1_id, bond.a2_id}:
            continue
        has_triple = has_triple or int(getattr(bond, "order", 1) or 1) >= 3
        has_double = has_double or int(getattr(bond, "order", 1) or 1) == 2 or bool(getattr(bond, "is_aromatic", False))
    atom = graph.atoms.get(atom_id)
    if has_triple:
        return math.radians(6.0)
    if has_double or (atom is not None and bool(getattr(atom, "is_aromatic", False))):
        return math.radians(8.0)
    return math.radians(15.0)


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


def _apply_forces(
    positions: dict[int, tuple[float, float]],
    forces: dict[int, tuple[float, float]],
    *,
    damping: float,
    max_step: float,
) -> None:
    damping = max(0.0, min(1.0, float(damping)))
    max_step = max(0.1, float(max_step))
    for atom_id, (fx, fy) in forces.items():
        step_x = fx * damping
        step_y = fy * damping
        step_len = math.hypot(step_x, step_y)
        if step_len <= 1e-9:
            continue
        if step_len > max_step:
            scale = max_step / step_len
            step_x *= scale
            step_y *= scale
        x, y = positions[atom_id]
        positions[atom_id] = (x + step_x, y + step_y)


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


def _deterministic_unit_vector(a_id: int, b_id: int) -> tuple[float, float, float]:
    seed = ((int(a_id) * 1103515245) ^ (int(b_id) * 12345)) & 0xFFFF
    angle = (float(seed) / 65535.0) * math.tau
    return math.cos(angle), math.sin(angle), 1.0


def _max_position_delta(
    before: dict[int, tuple[float, float]],
    after: dict[int, tuple[float, float]],
) -> float:
    max_delta = 0.0
    for atom_id, (x0, y0) in before.items():
        x1, y1 = after.get(atom_id, (x0, y0))
        max_delta = max(max_delta, math.hypot(x1 - x0, y1 - y0))
    return max_delta


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
