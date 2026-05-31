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


def structure_preserving_length_polish(
    positions: dict[int, tuple[float, float]],
    bonds: list[Any],
    target_len: float,
    *,
    max_iterations: int = 4,
    max_displacement_per_atom: float | None = None,
    tolerance: float = 0.02,
) -> dict[int, tuple[float, float]]:
    """Normaliza enlaces acíclicos moviendo sólo el lado móvil del enlace.

    A diferencia de ``length_only_polish``, no reparte fuerzas por toda la
    molécula. Los anillos quedan anclados y las ramas/subcadenas se trasladan
    como bloques rígidos para preservar ángulos y conformación local.
    """
    if not positions or not bonds:
        return {}

    atom_ids = set(positions)
    clean_bonds = [
        bond for bond in bonds
        if _bond_a1_helper(bond) in atom_ids and _bond_a2_helper(bond) in atom_ids
    ]
    if not clean_bonds:
        return {}

    target_len = max(8.0, float(target_len))
    max_disp = (
        float(max_displacement_per_atom)
        if max_displacement_per_atom is not None
        else target_len * 8.0
    )
    ring_atoms = _ring_atoms_from_bonds(atom_ids, clean_bonds)
    adjacency = _adjacency(atom_ids, clean_bonds)
    out = {aid: (float(x), float(y)) for aid, (x, y) in positions.items()}

    for _iteration in range(max(1, int(max_iterations))):
        max_step = 0.0
        for bond in clean_bonds:
            a1 = _bond_a1_helper(bond)
            a2 = _bond_a2_helper(bond)
            side, anchor, root = _movable_side_for_bond(
                a1, a2, atom_ids, adjacency, ring_atoms,
            )
            if not side or anchor not in out or root not in out:
                continue

            ax, ay = out[anchor]
            rx, ry = out[root]
            dx = rx - ax
            dy = ry - ay
            dist = math.hypot(dx, dy)
            if dist < 1e-6:
                continue
            desired = _target_length_for_order(bond, target_len)
            if abs(dist - desired) <= max(0.25, desired * tolerance):
                continue

            new_rx = ax + dx / dist * desired
            new_ry = ay + dy / dist * desired
            shift_x = new_rx - rx
            shift_y = new_ry - ry
            step = math.hypot(shift_x, shift_y)
            if step <= 1e-6:
                continue
            for atom_id in side:
                ox, oy = out[atom_id]
                out[atom_id] = (ox + shift_x, oy + shift_y)
            max_step = max(max_step, step)
        if max_step < 0.1:
            break

    changed: dict[int, tuple[float, float]] = {}
    for aid, (x, y) in out.items():
        ox, oy = positions.get(aid, (x, y))
        dx = x - ox
        dy = y - oy
        displacement = math.hypot(dx, dy)
        if displacement <= 0.5:
            continue
        if displacement > max_disp:
            scale = max_disp / displacement
            x = ox + dx * scale
            y = oy + dy * scale
        changed[aid] = (x, y)
    return changed


def structure_preserving_geometry_polish(
    positions: dict[int, tuple[float, float]],
    bonds: list[Any],
    target_len: float,
    *,
    max_displacement_per_atom: float | None = None,
) -> dict[int, tuple[float, float]]:
    """Reconstruye cadenas/ramas externas manteniendo anillos como bloques.

    La función comprime cada sistema cíclico en un cuerpo rígido y reconstruye
    el árbol acíclico que conecta esos bloques con ángulos químicos simples
    (sp3≈109.5, sp2/aromático≈120, sp≈180) y longitudes por orden de enlace.
    """
    if not positions or not bonds:
        return {}

    atom_ids = set(positions)
    clean_bonds = [
        bond for bond in bonds
        if _bond_a1_helper(bond) in atom_ids and _bond_a2_helper(bond) in atom_ids
    ]
    if not clean_bonds:
        return {}

    target_len = max(8.0, float(target_len))
    max_disp = (
        float(max_displacement_per_atom)
        if max_displacement_per_atom is not None
        else target_len * 8.0
    )

    groups, atom_to_group = _rigid_groups_for_geometry(atom_ids, clean_bonds, positions)
    if len(groups) <= 1:
        return {}

    group_edges = _group_edges_for_geometry(clean_bonds, atom_to_group)
    if not group_edges:
        return {}

    group_adjacency = _group_adjacency(group_edges)
    out = {aid: (float(x), float(y)) for aid, (x, y) in positions.items()}

    visited_groups: set[int] = set()
    for start_gid in sorted(groups):
        if start_gid in visited_groups:
            continue
        component = _group_component(start_gid, group_adjacency)
        visited_groups.update(component)
        if len(component) <= 1:
            continue
        if _group_component_has_cycle(component, group_adjacency):
            continue
        rebuilt = _layout_rigid_group_tree(
            component,
            groups,
            group_adjacency,
            clean_bonds,
            target_len,
        )
        for atom_id, coord in rebuilt.items():
            out[atom_id] = coord

    changed: dict[int, tuple[float, float]] = {}
    for aid, (x, y) in out.items():
        ox, oy = positions.get(aid, (x, y))
        dx = x - ox
        dy = y - oy
        displacement = math.hypot(dx, dy)
        if displacement <= 0.5:
            continue
        if displacement > max_disp:
            scale = max_disp / displacement
            x = ox + dx * scale
            y = oy + dy * scale
        changed[aid] = (x, y)
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


def _target_angle_for_atom(atom_id: int, bonds: list[Any]) -> float:
    has_triple = False
    has_double = False
    has_aromatic = False
    for bond in bonds:
        a1 = _bond_a1_helper(bond)
        a2 = _bond_a2_helper(bond)
        if atom_id not in {a1, a2}:
            continue
        order = int(bond.order if hasattr(bond, "order") else bond.get("order", 1) or 1)
        has_triple = has_triple or order >= 3
        has_double = has_double or order == 2
        has_aromatic = has_aromatic or bool(
            bond.is_aromatic if hasattr(bond, "is_aromatic") else bond.get("is_aromatic", False)
        )
    if has_triple:
        return math.radians(180.0)
    if has_double or has_aromatic:
        return math.radians(120.0)
    return math.radians(109.5)


def _adjacency(atom_ids: set[int], bonds: list[Any]) -> dict[int, set[int]]:
    adjacency: dict[int, set[int]] = {atom_id: set() for atom_id in atom_ids}
    for bond in bonds:
        a1 = _bond_a1_helper(bond)
        a2 = _bond_a2_helper(bond)
        if a1 not in atom_ids or a2 not in atom_ids:
            continue
        adjacency.setdefault(a1, set()).add(a2)
        adjacency.setdefault(a2, set()).add(a1)
    return adjacency


def _cycle_bond_edges(atom_ids: set[int], bonds: list[Any]) -> set[tuple[int, int]]:
    adjacency = _adjacency(atom_ids, bonds)
    edges: set[tuple[int, int]] = set()
    for bond in bonds:
        a1 = _bond_a1_helper(bond)
        a2 = _bond_a2_helper(bond)
        side = _component_without_edge(a1, a1, a2, atom_ids, adjacency)
        if a2 in side:
            edges.add((min(a1, a2), max(a1, a2)))
    return edges


def _component_without_edge(
    start: int,
    a1: int,
    a2: int,
    atom_ids: set[int],
    adjacency: dict[int, set[int]],
) -> set[int]:
    if start not in atom_ids:
        return set()
    visited: set[int] = set()
    stack = [start]
    blocked = {a1, a2}
    while stack:
        atom_id = stack.pop()
        if atom_id in visited:
            continue
        visited.add(atom_id)
        for neigh in adjacency.get(atom_id, set()):
            if {atom_id, neigh} == blocked:
                continue
            if neigh in atom_ids and neigh not in visited:
                stack.append(neigh)
    return visited


def _ring_atoms_from_bonds(atom_ids: set[int], bonds: list[Any]) -> set[int]:
    ring_atoms: set[int] = set()
    for a1, a2 in _cycle_bond_edges(atom_ids, bonds):
        ring_atoms.update({a1, a2})
    return ring_atoms


def _movable_side_for_bond(
    a1: int,
    a2: int,
    atom_ids: set[int],
    adjacency: dict[int, set[int]],
    ring_atoms: set[int],
) -> tuple[set[int], int, int]:
    left = _component_without_edge(a1, a1, a2, atom_ids, adjacency)
    right = _component_without_edge(a2, a1, a2, atom_ids, adjacency)
    if not left or not right or a2 in left or a1 in right:
        return set(), a1, a2

    left_ring = len(left & ring_atoms)
    right_ring = len(right & ring_atoms)
    if left_ring == 0 and right_ring > 0:
        return left, a2, a1
    if right_ring == 0 and left_ring > 0:
        return right, a1, a2
    if left_ring != right_ring:
        return (left, a2, a1) if left_ring < right_ring else (right, a1, a2)
    if len(left) != len(right):
        return (left, a2, a1) if len(left) < len(right) else (right, a1, a2)
    return (left, a2, a1) if min(left) > min(right) else (right, a1, a2)


def _rigid_groups_for_geometry(
    atom_ids: set[int],
    bonds: list[Any],
    positions: dict[int, tuple[float, float]],
) -> tuple[dict[int, dict[str, Any]], dict[int, int]]:
    cycle_edges = _cycle_bond_edges(atom_ids, bonds)
    ring_adj: dict[int, set[int]] = {atom_id: set() for atom_id in atom_ids}
    for a1, a2 in cycle_edges:
        ring_adj.setdefault(a1, set()).add(a2)
        ring_adj.setdefault(a2, set()).add(a1)

    groups: dict[int, dict[str, Any]] = {}
    atom_to_group: dict[int, int] = {}
    visited: set[int] = set()
    next_gid = 1

    for atom_id in sorted(atom_ids):
        if atom_id in visited or not ring_adj.get(atom_id):
            continue
        stack = [atom_id]
        component: set[int] = set()
        while stack:
            current = stack.pop()
            if current in component:
                continue
            component.add(current)
            visited.add(current)
            stack.extend(ring_adj.get(current, set()))
        if len(component) < 3:
            continue
        center = _center_for_atoms(component, positions)
        groups[next_gid] = {
            "atoms": component,
            "is_ring": True,
            "center": center,
            "relative": {
                aid: (positions[aid][0] - center[0], positions[aid][1] - center[1])
                for aid in component
            },
        }
        for aid in component:
            atom_to_group[aid] = next_gid
        next_gid += 1

    for atom_id in sorted(atom_ids):
        if atom_id in atom_to_group:
            continue
        x, y = positions[atom_id]
        groups[next_gid] = {
            "atoms": {atom_id},
            "is_ring": False,
            "center": (float(x), float(y)),
            "relative": {atom_id: (0.0, 0.0)},
        }
        atom_to_group[atom_id] = next_gid
        next_gid += 1

    return groups, atom_to_group


def _group_edges_for_geometry(
    bonds: list[Any],
    atom_to_group: dict[int, int],
) -> list[dict[str, Any]]:
    edges: list[dict[str, Any]] = []
    for bond in bonds:
        a1 = _bond_a1_helper(bond)
        a2 = _bond_a2_helper(bond)
        g1 = atom_to_group.get(a1)
        g2 = atom_to_group.get(a2)
        if g1 is None or g2 is None or g1 == g2:
            continue
        edges.append({
            "g1": g1,
            "g2": g2,
            "a1": a1,
            "a2": a2,
            "bond": bond,
        })
    return edges


def _group_adjacency(group_edges: list[dict[str, Any]]) -> dict[int, list[dict[str, Any]]]:
    adjacency: dict[int, list[dict[str, Any]]] = {}
    for edge in group_edges:
        adjacency.setdefault(edge["g1"], []).append(edge)
        adjacency.setdefault(edge["g2"], []).append(edge)
    return adjacency


def _group_component(start_gid: int, adjacency: dict[int, list[dict[str, Any]]]) -> set[int]:
    component: set[int] = set()
    stack = [start_gid]
    while stack:
        gid = stack.pop()
        if gid in component:
            continue
        component.add(gid)
        for edge in adjacency.get(gid, []):
            stack.append(_other_group(edge, gid))
    return component


def _group_component_has_cycle(
    component: set[int],
    adjacency: dict[int, list[dict[str, Any]]],
) -> bool:
    edge_keys: set[tuple[int, int]] = set()
    for gid in component:
        for edge in adjacency.get(gid, []):
            other = _other_group(edge, gid)
            if other in component:
                edge_keys.add((min(gid, other), max(gid, other)))
    return len(edge_keys) >= len(component)


def _layout_rigid_group_tree(
    component: set[int],
    groups: dict[int, dict[str, Any]],
    adjacency: dict[int, list[dict[str, Any]]],
    bonds: list[Any],
    target_len: float,
) -> dict[int, tuple[float, float]]:
    root = _choose_geometry_root(component, groups)
    placed_centers: dict[int, tuple[float, float]] = {
        root: groups[root]["center"],
    }
    used_angles: dict[int, list[float]] = {root: []}
    parent: dict[int, int | None] = {root: None}
    queue = [root]

    while queue:
        gid = queue.pop(0)
        pending = [
            edge for edge in adjacency.get(gid, [])
            if _other_group(edge, gid) in component
            and _other_group(edge, gid) not in placed_centers
        ]
        pending.sort(key=lambda edge: _geometry_edge_priority(edge, gid, groups))
        for edge in pending:
            child = _other_group(edge, gid)
            parent_atom, child_atom = _edge_atoms_for_parent(edge, gid)
            parent_anchor = _group_atom_position(groups[gid], placed_centers[gid], parent_atom)
            old_parent_anchor = _group_atom_position(
                groups[gid],
                groups[gid]["center"],
                parent_atom,
            )
            old_child_anchor = _group_atom_position(
                groups[child],
                groups[child]["center"],
                child_atom,
            )
            old_angle = math.atan2(
                old_child_anchor[1] - old_parent_anchor[1],
                old_child_anchor[0] - old_parent_anchor[0],
            )
            angle = _choose_geometry_angle(
                gid,
                parent_atom,
                old_angle,
                groups,
                placed_centers,
                used_angles.setdefault(gid, []),
                bonds,
            )
            length = _target_length_for_order(edge["bond"], target_len)
            child_anchor = (
                parent_anchor[0] + math.cos(angle) * length,
                parent_anchor[1] + math.sin(angle) * length,
            )
            child_rel = groups[child]["relative"].get(child_atom, (0.0, 0.0))
            placed_centers[child] = (
                child_anchor[0] - child_rel[0],
                child_anchor[1] - child_rel[1],
            )
            used_angles[gid].append(angle)
            used_angles[child] = [(angle + math.pi) % math.tau]
            parent[child] = gid
            queue.append(child)

    rebuilt: dict[int, tuple[float, float]] = {}
    for gid in component:
        if gid not in placed_centers:
            continue
        center = placed_centers[gid]
        for atom_id, (rx, ry) in groups[gid]["relative"].items():
            rebuilt[atom_id] = (center[0] + rx, center[1] + ry)
    return rebuilt


def _choose_geometry_root(
    component: set[int],
    groups: dict[int, dict[str, Any]],
) -> int:
    ring_groups = [gid for gid in component if groups[gid]["is_ring"]]
    candidates = ring_groups if ring_groups else list(component)
    return min(
        candidates,
        key=lambda gid: (
            -len(groups[gid]["atoms"]),
            groups[gid]["center"][0],
            groups[gid]["center"][1],
            gid,
        ),
    )


def _geometry_edge_priority(
    edge: dict[str, Any],
    gid: int,
    groups: dict[int, dict[str, Any]],
) -> tuple[int, float, int]:
    child = _other_group(edge, gid)
    is_ring = 1 if groups[child]["is_ring"] else 0
    cx, cy = groups[child]["center"]
    return is_ring, cx + cy * 0.001, child


def _choose_geometry_angle(
    gid: int,
    atom_id: int,
    old_angle: float,
    groups: dict[int, dict[str, Any]],
    placed_centers: dict[int, tuple[float, float]],
    existing_angles: list[float],
    bonds: list[Any],
) -> float:
    group = groups[gid]
    if group["is_ring"]:
        center = placed_centers[gid]
        anchor = _group_atom_position(group, center, atom_id)
        outward = math.atan2(anchor[1] - center[1], anchor[0] - center[0])
        if math.hypot(anchor[0] - center[0], anchor[1] - center[1]) > 1e-6:
            candidates = [outward, outward + math.radians(30.0), outward - math.radians(30.0)]
        else:
            candidates = [old_angle]
    elif not existing_angles:
        candidates = [old_angle]
    else:
        target = _target_angle_for_atom(atom_id, bonds)
        candidates = []
        if target >= math.radians(165.0):
            candidates.extend(angle + math.pi for angle in existing_angles)
        else:
            for angle in existing_angles:
                candidates.append(angle + target)
                candidates.append(angle - target)

    filtered = [
        angle for angle in candidates
        if all(abs(_angle_delta(angle - used)) >= math.radians(24.0) for used in existing_angles)
    ]
    viable = filtered or candidates or [old_angle]
    return min(viable, key=lambda angle: abs(_angle_delta(angle - old_angle)))


def _group_atom_position(
    group: dict[str, Any],
    center: tuple[float, float],
    atom_id: int,
) -> tuple[float, float]:
    rx, ry = group["relative"].get(atom_id, (0.0, 0.0))
    return center[0] + rx, center[1] + ry


def _edge_atoms_for_parent(edge: dict[str, Any], gid: int) -> tuple[int, int]:
    if edge["g1"] == gid:
        return edge["a1"], edge["a2"]
    return edge["a2"], edge["a1"]


def _other_group(edge: dict[str, Any], gid: int) -> int:
    return edge["g2"] if edge["g1"] == gid else edge["g1"]


def _center_for_atoms(
    atom_ids: set[int],
    positions: dict[int, tuple[float, float]],
) -> tuple[float, float]:
    if not atom_ids:
        return 0.0, 0.0
    return (
        sum(positions[aid][0] for aid in atom_ids) / len(atom_ids),
        sum(positions[aid][1] for aid in atom_ids) / len(atom_ids),
    )


def _angle_delta(value: float) -> float:
    while value <= -math.pi:
        value += math.tau
    while value > math.pi:
        value -= math.tau
    return value
