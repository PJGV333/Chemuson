from __future__ import annotations

import copy
import math
from typing import Iterable

from chemuson.core.model import Bond, MolGraph


Point2D = tuple[float, float]
CoordMap = dict[int, Point2D]


def normalize_atom_ids(graph: MolGraph, atom_ids: Iterable[int] | None) -> set[int]:
    """Return existing atom ids, preserving the all-atoms default used by Clean2D."""
    if atom_ids is None:
        return set(graph.atoms)
    return {int(atom_id) for atom_id in atom_ids if int(atom_id) in graph.atoms}


def graph_atom_coords(graph: MolGraph, atom_ids: Iterable[int] | None = None) -> CoordMap:
    """Read 2D atom coordinates as floats for scoring/layout code."""
    selected = normalize_atom_ids(graph, atom_ids)
    return {
        atom_id: (float(graph.atoms[atom_id].x), float(graph.atoms[atom_id].y))
        for atom_id in selected
    }


def graph_with_coords(graph: MolGraph, coords: CoordMap) -> MolGraph:
    """Return a copy of graph with replacement coordinates for scoring only."""
    out = copy.deepcopy(graph)
    for atom_id, (x, y) in coords.items():
        if atom_id in out.atoms:
            out.atoms[atom_id].x = float(x)
            out.atoms[atom_id].y = float(y)
    return out


def apply_coords_in_place(graph: MolGraph, coords: CoordMap) -> None:
    """Apply coordinates to an already-copied graph."""
    for atom_id, (x, y) in coords.items():
        if atom_id in graph.atoms:
            graph.atoms[atom_id].x = float(x)
            graph.atoms[atom_id].y = float(y)


def finite_coords(coords: CoordMap) -> bool:
    return all(math.isfinite(x) and math.isfinite(y) for x, y in coords.values())


def centroid(coords: CoordMap, atom_ids: Iterable[int]) -> Point2D | None:
    points = [coords[atom_id] for atom_id in atom_ids if atom_id in coords]
    if not points:
        return None
    return sum(x for x, _y in points) / len(points), sum(y for _x, y in points) / len(points)


def adjacency_from_bonds(atom_ids: Iterable[int], bonds: Iterable[Bond]) -> dict[int, set[int]]:
    adjacency: dict[int, set[int]] = {int(atom_id): set() for atom_id in atom_ids}
    for bond in bonds:
        adjacency.setdefault(int(bond.a1_id), set()).add(int(bond.a2_id))
        adjacency.setdefault(int(bond.a2_id), set()).add(int(bond.a1_id))
    return adjacency


def count_crossings(coords: CoordMap, bonds: list[Bond]) -> int:
    crossings = 0
    for i, b1 in enumerate(bonds):
        if b1.a1_id not in coords or b1.a2_id not in coords:
            continue
        for b2 in bonds[i + 1:]:
            if {b1.a1_id, b1.a2_id} & {b2.a1_id, b2.a2_id}:
                continue
            if b2.a1_id not in coords or b2.a2_id not in coords:
                continue
            if segments_intersect(coords[b1.a1_id], coords[b1.a2_id], coords[b2.a1_id], coords[b2.a2_id]):
                crossings += 1
    return crossings


def segments_intersect(p1: Point2D, p2: Point2D, p3: Point2D, p4: Point2D) -> bool:
    def orient(a: Point2D, b: Point2D, c: Point2D) -> float:
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    o1 = orient(p1, p2, p3)
    o2 = orient(p1, p2, p4)
    o3 = orient(p3, p4, p1)
    o4 = orient(p3, p4, p2)
    if abs(o1) < 1e-9 or abs(o2) < 1e-9 or abs(o3) < 1e-9 or abs(o4) < 1e-9:
        return False
    return (o1 > 0) != (o2 > 0) and (o3 > 0) != (o4 > 0)


def cycle_basis(atom_ids: set[int], bonds: list[Bond], *, max_size: int) -> list[list[int]]:
    """Small cycle basis approximation used by geometric safety/scoring."""
    adjacency = adjacency_from_bonds(atom_ids, bonds)
    rings: list[list[int]] = []
    seen: set[tuple[int, ...]] = set()
    for bond in sorted(bonds, key=lambda item: item.id):
        path = shortest_path_without_edge(adjacency, bond.a1_id, bond.a2_id)
        if path is None or not (3 <= len(path) <= max_size):
            continue
        key = tuple(sorted(path))
        if key in seen:
            continue
        seen.add(key)
        rings.append(path)
    return rings


def shortest_path_without_edge(
    adjacency: dict[int, set[int]],
    start: int,
    end: int,
) -> list[int] | None:
    blocked = {start, end}
    queue: list[tuple[int, list[int]]] = [(start, [start])]
    seen: set[int] = set()
    while queue:
        atom_id, path = queue.pop(0)
        if atom_id == end and len(path) >= 3:
            return path
        if atom_id in seen:
            continue
        seen.add(atom_id)
        for neighbor in sorted(adjacency.get(atom_id, set())):
            if {atom_id, neighbor} == blocked or neighbor in path:
                continue
            queue.append((neighbor, path + [neighbor]))
    return None
