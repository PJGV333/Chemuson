from __future__ import annotations

from collections import deque
from typing import Dict, List, Optional, Set, Tuple

from .molview import MolView


def carbon_skeleton_nodes(view: MolView) -> Set[int]:
    """Return atom ids that are carbon atoms."""
    return {atom_id for atom_id in view.atoms() if view.element(atom_id) == "C"}


def longest_carbon_chain(view: MolView) -> List[int]:
    """Return the longest carbon chain in an acyclic graph.

    Uses a deterministic tie-break: lexicographically smallest atom id tuple.
    """
    carbon_nodes = sorted(carbon_skeleton_nodes(view))
    if not carbon_nodes:
        return []

    adjacency: Dict[int, List[int]] = {node: [] for node in carbon_nodes}
    carbon_set = set(carbon_nodes)
    for node in carbon_nodes:
        neighbors = [nbr for nbr in view.neighbors(node) if nbr in carbon_set]
        adjacency[node] = sorted(neighbors)

    best_len = -1
    best_tuple: Optional[Tuple[int, ...]] = None
    best_path: List[int] = []

    for start in carbon_nodes:
        dist, parent = _bfs(start, adjacency)
        for end, d in dist.items():
            if d < best_len:
                continue
            path = _reconstruct_path(parent, start, end)
            if not path:
                continue
            path_tuple = tuple(path)
            rev_tuple = tuple(reversed(path))
            canon = min(path_tuple, rev_tuple)
            if d > best_len or best_tuple is None or canon < best_tuple:
                best_len = d
                best_tuple = canon
                best_path = list(canon)

    return best_path


def longest_chain_in_subset(view: MolView, allowed_nodes: Set[int]) -> List[int]:
    carbon_nodes = sorted(node for node in allowed_nodes if view.element(node) == "C")
    if not carbon_nodes:
        return []
    adjacency: Dict[int, List[int]] = {node: [] for node in carbon_nodes}
    carbon_set = set(carbon_nodes)
    for node in carbon_nodes:
        neighbors = [nbr for nbr in view.neighbors(node) if nbr in carbon_set]
        adjacency[node] = sorted(neighbors)

    best_len = -1
    best_tuple: Optional[Tuple[int, ...]] = None
    best_path: List[int] = []

    for start in carbon_nodes:
        dist, parent = _bfs(start, adjacency)
        for end, d in dist.items():
            if d < best_len:
                continue
            path = _reconstruct_path(parent, start, end)
            if not path:
                continue
            path_tuple = tuple(path)
            rev_tuple = tuple(reversed(path))
            canon = min(path_tuple, rev_tuple)
            if d > best_len or best_tuple is None or canon < best_tuple:
                best_len = d
                best_tuple = canon
                best_path = list(canon)
    return best_path


def _bfs(start: int, adjacency: Dict[int, List[int]]) -> Tuple[Dict[int, int], Dict[int, Optional[int]]]:
    dist: Dict[int, int] = {start: 0}
    parent: Dict[int, Optional[int]] = {start: None}
    queue: deque[int] = deque([start])
    while queue:
        node = queue.popleft()
        for nbr in adjacency.get(node, []):
            if nbr in dist:
                continue
            dist[nbr] = dist[node] + 1
            parent[nbr] = node
            queue.append(nbr)
    return dist, parent


def _reconstruct_path(
    parent: Dict[int, Optional[int]], start: int, end: int
) -> List[int]:
    if end not in parent:
        return []
    path: List[int] = []
    current: Optional[int] = end
    while current is not None:
        path.append(current)
        if current == start:
            break
        current = parent.get(current)
    if not path or path[-1] != start:
        return []
    return list(reversed(path))
