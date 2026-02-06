"""Selección de cadena principal para nomenclatura IUPAC-lite.

Incluye funciones para identificar cadenas de carbono y elegir la más larga
con criterios deterministas de desempate.
"""

from __future__ import annotations

from collections import deque
from typing import Dict, List, Optional, Set, Tuple

from .molview import MolView


def carbon_skeleton_nodes(view: MolView) -> Set[int]:
    """Devuelve los IDs de átomos que son carbono.

    Args:
        view: Vista del grafo molecular.

    Returns:
        Conjunto de IDs de átomos de carbono.
    """
    return {atom_id for atom_id in view.atoms() if view.element(atom_id) == "C"}


def longest_carbon_chain(view: MolView) -> List[int]:
    """Encuentra la cadena de carbono más larga en un grafo acíclico.

    En caso de empate, se elige la tupla de IDs lexicográficamente menor
    para garantizar determinismo.

    Args:
        view: Vista del grafo molecular.

    Returns:
        Lista de IDs de átomos que forman la cadena principal.
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
        # BFS desde cada nodo para obtener el diámetro de la subred de carbonos.
        dist, parent = _bfs(start, adjacency)
        for end, d in dist.items():
            if d < best_len:
                continue
            path = _reconstruct_path(parent, start, end)
            if not path:
                continue
            path_tuple = tuple(path)
            rev_tuple = tuple(reversed(path))
            # Canonicalizamos la cadena para desempates estables.
            canon = min(path_tuple, rev_tuple)
            if d > best_len or best_tuple is None or canon < best_tuple:
                best_len = d
                best_tuple = canon
                best_path = list(canon)

    return best_path


def longest_chain_in_subset(view: MolView, allowed_nodes: Set[int]) -> List[int]:
    """Devuelve la cadena de carbono más larga dentro de un subconjunto.

    Args:
        view: Vista del grafo molecular.
        allowed_nodes: Conjunto de IDs permitidos para la cadena.

    Returns:
        Lista de IDs de la cadena más larga encontrada.
    """
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
    """Búsqueda BFS para calcular distancias y padres.

    Args:
        start: Nodo inicial.
        adjacency: Lista de adyacencia.

    Returns:
        Diccionarios de distancia y padre por nodo.
    """
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
    """Reconstruye el camino desde `start` hasta `end` usando padres.

    Args:
        parent: Diccionario de padres generado por BFS.
        start: Nodo inicial.
        end: Nodo final.

    Returns:
        Lista con el camino reconstruido o lista vacía si no existe.
    """
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
