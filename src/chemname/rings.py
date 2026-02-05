from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Set, Tuple

from chemcalc.valence import implicit_h_count

from .molview import MolView


def find_rings_simple(view: MolView) -> List[frozenset[int]]:
    """Find simple rings using shortest-path edge removal.

    Returns a list of ring atom sets. Best effort for non-fused rings.
    """
    adjacency = _build_adjacency(view)
    edges: Set[Tuple[int, int]] = set()
    for a, nbrs in adjacency.items():
        for b in nbrs:
            if a < b:
                edges.add((a, b))

    rings: Dict[frozenset[int], List[int]] = {}
    for a, b in edges:
        path = _shortest_path_excluding_edge(adjacency, a, b, blocked_edge=(a, b))
        if not path:
            continue
        ring_nodes = frozenset(path)
        if len(ring_nodes) < 3:
            continue
        if ring_nodes not in rings:
            rings[ring_nodes] = path

    return _sorted_rings(rings.keys())


def has_ring(view: MolView) -> bool:
    return bool(find_rings_simple(view))


def is_simple_ring(view: MolView, ring_nodes: Iterable[int]) -> bool:
    ring_set = set(ring_nodes)
    if len(ring_set) < 3:
        return False
    for atom_id in ring_set:
        neighbors = [nbr for nbr in view.neighbors(atom_id) if nbr in ring_set]
        if len(neighbors) != 2:
            return False
    return True


def ring_bonds(view: MolView, ring_nodes: Iterable[int]) -> Set[frozenset[int]]:
    ring_set = set(ring_nodes)
    bonds: Set[frozenset[int]] = set()
    for atom_id in ring_set:
        for nbr in view.neighbors(atom_id):
            if nbr in ring_set:
                bonds.add(frozenset({atom_id, nbr}))
    return bonds


def ring_order(view: MolView, ring_nodes: Iterable[int]) -> List[int]:
    ring_set = set(ring_nodes)
    if not is_simple_ring(view, ring_set):
        return []
    adjacency = {node: [nbr for nbr in view.neighbors(node) if nbr in ring_set] for node in ring_set}
    start = min(ring_set)
    next_node = min(adjacency[start])
    order = [start, next_node]
    prev = start
    current = next_node
    while True:
        nbrs = adjacency[current]
        next_candidate = nbrs[0] if nbrs[1] == prev else nbrs[1]
        if next_candidate == start:
            break
        order.append(next_candidate)
        prev, current = current, next_candidate
        if len(order) > len(ring_set):
            return []
    if len(order) != len(ring_set):
        return []
    return order


def perceive_aromaticity_basic(view: MolView, rings: Iterable[frozenset[int]]) -> Set[frozenset[int]]:
    """Mark benzene-like rings as aromatic (conservative)."""
    aromatic: Set[frozenset[int]] = set()
    for ring_nodes in rings:
        if len(ring_nodes) != 6:
            continue
        if not is_simple_ring(view, ring_nodes):
            continue
        if any(view.element(atom_id) != "C" for atom_id in ring_nodes):
            continue
        if _ring_aromatic_basic(view, ring_nodes):
            aromatic.add(ring_nodes)
    return aromatic


def classify_aromatic_ring(view: MolView, ring_nodes: Iterable[int]) -> Optional[dict]:
    """Classify simple aromatic rings into retained-name kinds."""
    ring_set = set(ring_nodes)
    if not is_simple_ring(view, ring_set):
        return None
    order = ring_order(view, ring_set)
    if not order:
        return None
    size = len(order)
    if size not in {5, 6}:
        return None
    if not _ring_aromatic_basic(view, ring_set):
        return None

    hetero_atoms = [atom_id for atom_id in order if view.element(atom_id) != "C"]
    hetero_positions = [idx + 1 for idx, atom_id in enumerate(order) if atom_id in hetero_atoms]
    hetero_elements = [view.element(atom_id) for atom_id in hetero_atoms]

    kind = None
    if size == 6:
        if not hetero_atoms:
            kind = "benzene"
        elif len(hetero_atoms) == 1 and hetero_elements[0] == "N":
            kind = "pyridine"
        elif len(hetero_atoms) == 2 and all(elem == "N" for elem in hetero_elements):
            positions = sorted(hetero_positions)
            delta = positions[1] - positions[0]
            distance = min(delta, size - delta)
            if distance == 1:
                kind = "pyridazine"
            elif distance == 2:
                kind = "pyrimidine"
            elif distance == 3:
                kind = "pyrazine"
        elif len(hetero_atoms) == 3 and all(elem == "N" for elem in hetero_elements):
            positions = sorted(hetero_positions)
            deltas = []
            for idx in range(len(positions)):
                current = positions[idx]
                nxt = positions[(idx + 1) % len(positions)]
                deltas.append((nxt - current) % size)
            if sorted(deltas) == [2, 2, 2]:
                kind = "triazine"
    else:
        if len(hetero_atoms) == 1:
            elem = hetero_elements[0]
            if elem == "O":
                kind = "furan"
            elif elem == "S":
                kind = "thiophene"
            elif elem == "N":
                if _hetero_has_h(view, hetero_atoms[0]):
                    kind = "pyrrole"
        elif len(hetero_atoms) == 2:
            elements = sorted(hetero_elements)
            if elements == ["N", "N"]:
                if any(_hetero_has_h(view, atom_id) for atom_id in hetero_atoms):
                    kind = "imidazole"
            elif elements == ["N", "O"]:
                kind = "oxazole"
            elif elements == ["N", "S"]:
                kind = "thiazole"
        elif len(hetero_atoms) == 3:
            if all(elem == "N" for elem in hetero_elements):
                kind = "triazole"

    if kind is None:
        return None
    preferred_start = None
    if kind == "triazole":
        preferred_start = [atom_id for atom_id in hetero_atoms if _hetero_has_h(view, atom_id)]
    return {
        "kind": kind,
        "order": order,
        "hetero_atoms": hetero_atoms,
        "hetero_positions": hetero_positions,
        "preferred_start": preferred_start,
    }


def detect_naphthalene(
    view: MolView, rings: Iterable[frozenset[int]]
) -> Optional[dict]:
    ring_list = [ring for ring in rings if len(ring) == 6]
    if len(ring_list) < 2:
        return None
    candidates = []
    for i in range(len(ring_list)):
        for j in range(i + 1, len(ring_list)):
            r1 = ring_list[i]
            r2 = ring_list[j]
            if not _ring_aromatic_basic(view, r1) or not _ring_aromatic_basic(view, r2):
                continue
            if any(view.element(atom_id) != "C" for atom_id in r1 | r2):
                continue
            shared = r1 & r2
            if len(shared) != 2:
                continue
            a, b = tuple(shared)
            if view.bond_order_between(a, b) <= 0:
                continue
            union = r1 | r2
            if len(union) != 10:
                continue
            candidates.append((tuple(sorted(shared)), union, (r1, r2)))
    if not candidates:
        return None
    candidates.sort(key=lambda item: (tuple(sorted(item[1])), item[0]))
    fusion_atoms, union, rings_pair = candidates[0]
    return {
        "atoms": union,
        "fusion_atoms": fusion_atoms,
        "rings": rings_pair,
    }


def ring_type_basic(view: MolView, ring_nodes: Iterable[int]) -> Optional[str]:
    ring_set = set(ring_nodes)
    if not is_simple_ring(view, ring_set):
        return None
    size = len(ring_set)
    if size == 6 and all(view.element(atom_id) == "C" for atom_id in ring_set):
        if _ring_aromatic_basic(view, ring_set):
            return "benzene"
        if all(
            view.bond_order_between(
                ring_order(view, ring_set)[idx],
                ring_order(view, ring_set)[(idx + 1) % size],
            )
            == 1
            for idx in range(size)
        ):
            return "cyclohexane"
    return None


@dataclass(frozen=True)
class RingContext:
    rings: List[frozenset[int]]
    ring_types: Dict[frozenset[int], str]
    atom_rings: Dict[int, List[frozenset[int]]]


def build_ring_context(view: MolView) -> RingContext:
    rings = find_rings_simple(view)
    ring_types: Dict[frozenset[int], str] = {}
    atom_rings: Dict[int, List[frozenset[int]]] = {}
    for ring in rings:
        rtype = ring_type_basic(view, ring)
        if rtype:
            ring_types[ring] = rtype
        for atom_id in ring:
            atom_rings.setdefault(atom_id, []).append(ring)
    return RingContext(rings=rings, ring_types=ring_types, atom_rings=atom_rings)


def _build_adjacency(view: MolView) -> Dict[int, List[int]]:
    adjacency: Dict[int, List[int]] = {}
    for atom_id in view.atoms():
        adjacency[atom_id] = list(view.neighbors(atom_id))
    return adjacency


def _ring_aromatic_basic(view: MolView, ring_nodes: Iterable[int]) -> bool:
    ring_set = set(ring_nodes)
    bond_pairs = list(ring_bonds(view, ring_set))
    if not bond_pairs:
        return False
    if all(view.bond_is_aromatic(*tuple(pair)) for pair in bond_pairs):
        return True

    size = len(ring_set)
    double_bonds = 0
    double_bond_count_per_atom: Dict[int, int] = {atom_id: 0 for atom_id in ring_set}
    for pair in bond_pairs:
        a, b = tuple(pair)
        order = view.bond_order_between(a, b)
        if order not in {1, 2}:
            return False
        if order == 2:
            double_bonds += 1
            double_bond_count_per_atom[a] += 1
            double_bond_count_per_atom[b] += 1

    if size == 6:
        return double_bonds == 3 and all(count == 1 for count in double_bond_count_per_atom.values())
    if size == 5:
        return double_bonds == 2 and all(count <= 1 for count in double_bond_count_per_atom.values())
    return False


def _hetero_has_h(view: MolView, atom_id: int) -> bool:
    return implicit_h_count(view, atom_id) + view.explicit_h(atom_id) > 0


def _shortest_path_excluding_edge(
    adjacency: Dict[int, List[int]],
    start: int,
    goal: int,
    blocked_edge: Tuple[int, int],
) -> List[int]:
    blocked = {blocked_edge, (blocked_edge[1], blocked_edge[0])}
    queue: deque[int] = deque([start])
    parent: Dict[int, Optional[int]] = {start: None}
    while queue:
        node = queue.popleft()
        if node == goal:
            break
        for nbr in adjacency.get(node, []):
            if (node, nbr) in blocked:
                continue
            if nbr in parent:
                continue
            parent[nbr] = node
            queue.append(nbr)

    if goal not in parent:
        return []
    path: List[int] = []
    current: Optional[int] = goal
    while current is not None:
        path.append(current)
        current = parent.get(current)
    path.reverse()
    return path


def _sorted_rings(rings: Iterable[frozenset[int]]) -> List[frozenset[int]]:
    return sorted(rings, key=lambda ring: (len(ring), tuple(sorted(ring))))
