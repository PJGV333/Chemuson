from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Iterable

from chemuson.chemio.depiction_candidates import block_donut_score, score_imported_depiction
from chemuson.clean2d.geometry import (
    adjacency_from_bonds,
    centroid,
    count_crossings,
    cycle_basis,
    finite_coords,
    graph_atom_coords,
    graph_with_coords,
    normalize_atom_ids,
)
from chemuson.clean2d.local_graph_cleaner import stereo_layout_signature
from chemuson.clean2d.safety import min_nonbonded_distance, ring_degeneracy_score
from chemuson.core.layers import BlockKind, MotifKind, build_multilayer_chemical_graph
from chemuson.core.model import Bond, MolGraph, bond_affects_valence


@dataclass(frozen=True)
class BlockUnwrapReport:
    ok: bool
    reason: str = ""
    source: str = "chemuson_block_unwrap"
    atom_count: int = 0
    ring_count: int = 0
    rigid_block_count: int = 0
    linker_count: int = 0
    macrocycle_count: int = 0
    cyclophane_count: int = 0
    aromatic_ring_count: int = 0
    donut_score_before: float = 0.0
    donut_score_after: float = 0.0
    score_before: float = 0.0
    score_after: float = 0.0
    metadata: dict[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class _RigidUnit:
    id: int
    atom_ids: frozenset[int]
    centroid: tuple[float, float]


@dataclass(frozen=True)
class _UnitEdge:
    left: int
    right: int
    path: tuple[int, ...]
    left_anchor: int
    right_anchor: int


def block_unwrap_layout(
    graph: MolGraph,
    atom_ids: Iterable[int] | None = None,
    *,
    target_bond_length: float = 40.0,
) -> tuple[dict[int, tuple[float, float]] | None, BlockUnwrapReport]:
    selected = normalize_atom_ids(graph, atom_ids)
    target = max(8.0, float(target_bond_length or 40.0))
    before = graph_atom_coords(graph, selected)
    bonds = _selected_bonds(graph, selected)
    score_before, _ = score_imported_depiction(graph_with_coords(graph, before), target_bond_length=target)
    donut_before, donut_meta = block_donut_score(graph_with_coords(graph, before), target_bond_length=target)
    try:
        layers = build_multilayer_chemical_graph(graph, selected)
    except Exception as exc:
        return None, _report(False, "layers_failed", selected, score_before, score_before, donut_before, donut_before, {"error": str(exc)})

    ring_motifs = [motif for motif in layers.motif_graph.motifs if motif.kind == MotifKind.RING]
    block_counts: dict[BlockKind, int] = {}
    for block in layers.block_graph.blocks:
        block_counts[block.kind] = block_counts.get(block.kind, 0) + 1
    rigid_block_count = sum(
        block_counts.get(kind, 0)
        for kind in (BlockKind.AROMATIC_RING, BlockKind.FUSED_SYSTEM, BlockKind.MACROCYCLE, BlockKind.CYCLOPHANE)
    )
    linker_count = block_counts.get(BlockKind.LINKER, 0) or _fallback_linker_count(graph, selected, ring_motifs)
    trigger = (
        block_counts.get(BlockKind.MACROCYCLE, 0) > 0
        or block_counts.get(BlockKind.CYCLOPHANE, 0) > 0
        or block_counts.get(BlockKind.AROMATIC_RING, 0) >= 4
        or (rigid_block_count >= 4 and linker_count >= 2)
        or (len(selected) >= 35 and len(ring_motifs) >= 4)
        or donut_before >= 4.0
    )
    base_metadata = {"donut_before_metadata": donut_meta, "linker_count_effective": linker_count}
    if not trigger:
        return None, _report(False, "not_block_rich_or_donut", selected, score_before, score_before, donut_before, donut_before, base_metadata, ring_count=len(ring_motifs), rigid_block_count=rigid_block_count, linker_count=linker_count, block_counts=block_counts)

    units = _rigid_units(graph, selected, ring_motifs, before)
    if len(units) < 4:
        return None, _report(False, "too_few_rigid_units", selected, score_before, score_before, donut_before, donut_before, base_metadata, ring_count=len(ring_motifs), rigid_block_count=len(units), linker_count=linker_count, block_counts=block_counts)
    edges = _unit_edges(selected, bonds, units)
    if not edges:
        return None, _report(False, "no_block_connectivity", selected, score_before, score_before, donut_before, donut_before, base_metadata, ring_count=len(ring_motifs), rigid_block_count=len(units), linker_count=linker_count, block_counts=block_counts)

    coords = _layout_units_and_linkers(graph, selected, bonds, before, units, edges, target)
    coords = _relax_nonrigid(selected, bonds, coords, _rigid_atoms(units), target)
    if set(coords) != selected:
        return None, _report(False, "missing_coordinates", selected, score_before, score_before, donut_before, donut_before, base_metadata, ring_count=len(ring_motifs), rigid_block_count=len(units), linker_count=linker_count, block_counts=block_counts)
    if not finite_coords(coords):
        return None, _report(False, "non_finite_coordinates", selected, score_before, score_before, donut_before, donut_before, base_metadata, ring_count=len(ring_motifs), rigid_block_count=len(units), linker_count=linker_count, block_counts=block_counts)

    after_graph = graph_with_coords(graph, coords)
    score_after, score_meta_after = score_imported_depiction(after_graph, target_bond_length=target)
    donut_after, donut_meta_after = block_donut_score(after_graph, target_bond_length=target)
    safety_reason = _safety_rejection_reason(graph, selected, bonds, before, coords, target, score_before, score_after, donut_before, donut_after)
    metadata = {**base_metadata, "score_after_metadata": score_meta_after, "donut_after_metadata": donut_meta_after, "unit_count": len(units), "edge_count": len(edges)}
    if safety_reason:
        return None, _report(False, safety_reason, selected, score_before, score_after, donut_before, donut_after, metadata, ring_count=len(ring_motifs), rigid_block_count=len(units), linker_count=linker_count, block_counts=block_counts)
    return coords, _report(True, "ok", selected, score_before, score_after, donut_before, donut_after, metadata, ring_count=len(ring_motifs), rigid_block_count=len(units), linker_count=linker_count, block_counts=block_counts)


def _report(ok: bool, reason: str, selected: set[int], score_before: float, score_after: float, donut_before: float, donut_after: float, metadata: dict[str, object], *, ring_count: int = 0, rigid_block_count: int = 0, linker_count: int = 0, block_counts: dict[BlockKind, int] | None = None) -> BlockUnwrapReport:
    counts = block_counts or {}
    return BlockUnwrapReport(
        ok=ok,
        reason=reason,
        atom_count=len(selected),
        ring_count=ring_count,
        rigid_block_count=rigid_block_count,
        linker_count=linker_count,
        macrocycle_count=counts.get(BlockKind.MACROCYCLE, 0),
        cyclophane_count=counts.get(BlockKind.CYCLOPHANE, 0),
        aromatic_ring_count=counts.get(BlockKind.AROMATIC_RING, 0),
        donut_score_before=donut_before,
        donut_score_after=donut_after,
        score_before=score_before,
        score_after=score_after,
        metadata=metadata,
    )


def _rigid_units(graph: MolGraph, selected: set[int], ring_motifs, coords: dict[int, tuple[float, float]]) -> list[_RigidUnit]:
    groups: list[set[int]] = [set(motif.atom_ids) & selected for motif in ring_motifs if len(set(motif.atom_ids) & selected) >= 3]
    changed = True
    while changed:
        changed = False
        merged: list[set[int]] = []
        while groups:
            current = groups.pop(0)
            rest: list[set[int]] = []
            for group in groups:
                if len(current & group) >= 2:
                    current |= group
                    changed = True
                else:
                    rest.append(group)
            groups = rest
            merged.append(current)
        groups = merged
    units: list[_RigidUnit] = []
    for idx, atoms in enumerate(sorted(groups, key=lambda item: (min(item), len(item))), start=1):
        unit_centroid = centroid(coords, atoms)
        if unit_centroid is not None:
            units.append(_RigidUnit(idx, frozenset(atoms), unit_centroid))
    return units


def _unit_edges(selected: set[int], bonds: list[Bond], units: list[_RigidUnit]) -> list[_UnitEdge]:
    adjacency = adjacency_from_bonds(selected, bonds)
    atom_to_unit: dict[int, int] = {atom_id: unit.id for unit in units for atom_id in unit.atom_ids}
    edges: dict[tuple[int, int], _UnitEdge] = {}
    for unit in units:
        for atom_id in unit.atom_ids:
            for neighbor in adjacency.get(atom_id, set()):
                if atom_to_unit.get(neighbor) == unit.id:
                    continue
                path = _path_to_other_unit(neighbor, unit.id, atom_to_unit, adjacency, max_len=10)
                if path is None:
                    continue
                other = atom_to_unit[path[-1]]
                key = (min(unit.id, other), max(unit.id, other))
                full_path = (atom_id, *path)
                edge = _UnitEdge(key[0], key[1], full_path if unit.id == key[0] else tuple(reversed(full_path)), full_path[0] if unit.id == key[0] else full_path[-1], full_path[-1] if unit.id == key[0] else full_path[0])
                previous = edges.get(key)
                if previous is None or len(edge.path) < len(previous.path):
                    edges[key] = edge
    return list(edges.values())


def _layout_units_and_linkers(graph: MolGraph, selected: set[int], bonds: list[Bond], before: dict[int, tuple[float, float]], units: list[_RigidUnit], edges: list[_UnitEdge], target: float) -> dict[int, tuple[float, float]]:
    skeleton = _longest_unit_path(units, edges)
    unit_positions: dict[int, tuple[float, float]] = {}
    for idx, unit_id in enumerate(skeleton):
        unit_positions[unit_id] = (idx * target * 3.3, (((idx % 2) * 2) - 1) * target * 1.35)
    adjacency_units: dict[int, set[int]] = {unit.id: set() for unit in units}
    for edge in edges:
        adjacency_units.setdefault(edge.left, set()).add(edge.right)
        adjacency_units.setdefault(edge.right, set()).add(edge.left)
    queue = list(skeleton)
    while queue:
        parent = queue.pop(0)
        for child in sorted(adjacency_units.get(parent, set())):
            if child in unit_positions:
                continue
            px, py = unit_positions[parent]
            side = -1 if len(unit_positions) % 2 else 1
            unit_positions[child] = (px + target * 1.4, py + side * target * 2.8)
            queue.append(child)
    coords: dict[int, tuple[float, float]] = {}
    for unit in units:
        ux, uy = unit_positions.get(unit.id, unit.centroid)
        cx, cy = unit.centroid
        for atom_id in unit.atom_ids:
            bx, by = before[atom_id]
            coords[atom_id] = (ux + bx - cx, uy + by - cy)
    for edge in edges:
        path = edge.path
        linker = [atom_id for atom_id in path[1:-1] if atom_id not in coords]
        if not linker or path[0] not in coords or path[-1] not in coords:
            continue
        endpoint_distance = math.hypot(coords[path[-1]][0] - coords[path[0]][0], coords[path[-1]][1] - coords[path[0]][1])
        arch = endpoint_distance > target * 5.0 and len(linker) <= 4
        _place_linker(coords, path[0], path[-1], linker, target, arch=arch)
    rigid = set(coords)
    adjacency = adjacency_from_bonds(selected, bonds)
    for atom_id in sorted(selected):
        if atom_id in coords:
            continue
        anchors = [neighbor for neighbor in adjacency.get(atom_id, set()) if neighbor in coords]
        if anchors:
            ax, ay = centroid(coords, anchors) or (0.0, 0.0)
            nearest_unit_center = _nearest_point((ax, ay), list(unit_positions.values())) or (ax - target, ay)
            vx, vy = ax - nearest_unit_center[0], ay - nearest_unit_center[1]
            length = math.hypot(vx, vy) or 1.0
            coords[atom_id] = (ax + vx / length * target, ay + vy / length * target)
        else:
            bx, by = before[atom_id]
            coords[atom_id] = (bx, by)
    return coords


def _place_linker(coords: dict[int, tuple[float, float]], left: int, right: int, linker: list[int], target: float, *, arch: bool = False) -> None:
    x1, y1 = coords[left]
    x2, y2 = coords[right]
    dx, dy = x2 - x1, y2 - y1
    length = math.hypot(dx, dy) or 1.0
    px, py = -dy / length, dx / length
    for idx, atom_id in enumerate(linker, start=1):
        t = idx / (len(linker) + 1)
        if arch:
            offset = target * (2.2 + 0.25 * math.sin(math.pi * t))
        else:
            offset = ((-1) ** idx) * target * (0.25 + 0.08 * (idx % 2))
        coords[atom_id] = (x1 + dx * t + px * offset, y1 + dy * t + py * offset)


def _relax_nonrigid(selected: set[int], bonds: list[Bond], coords: dict[int, tuple[float, float]], rigid_atoms: set[int], target: float) -> dict[int, tuple[float, float]]:
    out = dict(coords)
    for _ in range(18):
        shifts = {atom_id: [0.0, 0.0] for atom_id in selected if atom_id not in rigid_atoms}
        for bond in bonds:
            if bond.a1_id not in out or bond.a2_id not in out:
                continue
            x1, y1 = out[bond.a1_id]
            x2, y2 = out[bond.a2_id]
            dx, dy = x2 - x1, y2 - y1
            dist = math.hypot(dx, dy) or 1.0
            desired = target * (0.96 if bond.order == 2 else 0.98 if bond.is_aromatic else 1.0)
            delta = (dist - desired) * 0.18
            fx, fy = dx / dist * delta, dy / dist * delta
            if bond.a1_id in shifts:
                shifts[bond.a1_id][0] += fx
                shifts[bond.a1_id][1] += fy
            if bond.a2_id in shifts:
                shifts[bond.a2_id][0] -= fx
                shifts[bond.a2_id][1] -= fy
        for atom_id, (sx, sy) in shifts.items():
            step = math.hypot(sx, sy)
            if step > target * 0.20:
                scale = target * 0.20 / step
                sx *= scale
                sy *= scale
            x, y = out[atom_id]
            out[atom_id] = (x + sx, y + sy)
    return out


def _safety_rejection_reason(graph: MolGraph, selected: set[int], bonds: list[Bond], before: dict[int, tuple[float, float]], after: dict[int, tuple[float, float]], target: float, score_before: float, score_after: float, donut_before: float, donut_after: float) -> str:
    if count_crossings(after, bonds) > count_crossings(before, bonds):
        return "new_bond_crossings"
    before_nonbonded = min_nonbonded_distance(before, bonds, selected)
    after_nonbonded = min_nonbonded_distance(after, bonds, selected)
    if after_nonbonded < min(before_nonbonded * 0.70, target * 0.22):
        return "nonbonded_distance_worse"
    before_ring = _min_ring_degeneracy(selected, bonds, before)
    after_ring = _min_ring_degeneracy(selected, bonds, after)
    if before_ring < math.inf and after_ring + 0.05 < before_ring:
        return "ring_degeneracy_worse"
    try:
        if stereo_layout_signature(graph, before, selected) != stereo_layout_signature(graph, after, selected):
            return "stereo_signature_changed"
    except Exception:
        pass
    if not (score_after <= score_before * 0.85 or score_before - score_after > target * 8.0):
        return "score_not_improved"
    if donut_after >= donut_before:
        return "donut_score_not_improved"
    return ""


def _min_ring_degeneracy(selected: set[int], bonds: list[Bond], coords: dict[int, tuple[float, float]]) -> float:
    rings = cycle_basis(selected, bonds, max_size=18)
    if not rings:
        return math.inf
    return min(ring_degeneracy_score(coords, set(ring)) for ring in rings)


def _longest_unit_path(units: list[_RigidUnit], edges: list[_UnitEdge]) -> list[int]:
    ids = [unit.id for unit in units]
    adjacency: dict[int, set[int]] = {unit_id: set() for unit_id in ids}
    for edge in edges:
        adjacency.setdefault(edge.left, set()).add(edge.right)
        adjacency.setdefault(edge.right, set()).add(edge.left)
    best = [ids[0]] if ids else []
    for start in ids:
        queue = [(start, [start])]
        while queue:
            current, path = queue.pop(0)
            if len(path) > len(best):
                best = path
            for neighbor in sorted(adjacency.get(current, set())):
                if neighbor not in path:
                    queue.append((neighbor, path + [neighbor]))
    return best


def _path_to_other_unit(start: int, own_unit: int, atom_to_unit: dict[int, int], adjacency: dict[int, set[int]], max_len: int) -> tuple[int, ...] | None:
    queue = [(start, [start])]
    seen: set[int] = set()
    while queue:
        atom_id, path = queue.pop(0)
        if len(path) > max_len:
            continue
        unit = atom_to_unit.get(atom_id)
        if unit is not None and unit != own_unit:
            return tuple(path)
        if atom_id in seen:
            continue
        seen.add(atom_id)
        for neighbor in sorted(adjacency.get(atom_id, set())):
            if neighbor in path:
                continue
            neighbor_unit = atom_to_unit.get(neighbor)
            if neighbor_unit == own_unit:
                continue
            queue.append((neighbor, path + [neighbor]))
    return None


def _fallback_linker_count(graph: MolGraph, selected: set[int], ring_motifs) -> int:
    ring_atoms = set().union(*(set(motif.atom_ids) & selected for motif in ring_motifs)) if ring_motifs else set()
    if not ring_atoms:
        return 0
    non_ring = selected - ring_atoms
    adjacency = adjacency_from_bonds(selected, _selected_bonds(graph, selected))
    count = 0
    seen: set[int] = set()
    for start in sorted(non_ring):
        if start in seen:
            continue
        stack = [start]
        component: set[int] = set()
        ring_neighbors: set[int] = set()
        while stack:
            atom_id = stack.pop()
            if atom_id in component:
                continue
            component.add(atom_id)
            for neighbor in adjacency.get(atom_id, set()):
                if neighbor in non_ring:
                    stack.append(neighbor)
                elif neighbor in ring_atoms:
                    ring_neighbors.add(neighbor)
        seen.update(component)
        if len(ring_neighbors) >= 2:
            count += 1
    return count


def _selected_bonds(graph: MolGraph, selected: set[int]) -> list[Bond]:
    return [bond for bond in graph.bonds.values() if bond.a1_id in selected and bond.a2_id in selected and bond_affects_valence(bond)]


def _rigid_atoms(units: list[_RigidUnit]) -> set[int]:
    return {atom_id for unit in units for atom_id in unit.atom_ids}


def _nearest_point(point: tuple[float, float], points: list[tuple[float, float]]) -> tuple[float, float] | None:
    if not points:
        return None
    return min(points, key=lambda other: math.hypot(other[0] - point[0], other[1] - point[1]))

