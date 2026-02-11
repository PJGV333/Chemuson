"""Utilidades de numeración visual para átomos y estructuras."""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from typing import Dict, List

from core.model import MolGraph


@dataclass(frozen=True)
class NumberedStructure:
    """Resultado de numeración para una estructura desconectada."""

    structure_id: int
    number: int
    x: float
    y: float
    atom_ids: tuple[int, ...]


def compute_atom_numbers(model: MolGraph, strategy: str = "xy") -> Dict[int, int]:
    """Calcula numeración de átomos estable y reproducible.

    Estrategia implementada:
    - ``xy`` por estructura: para cada componente conexa, ordena por ``x``
      ascendente y luego ``y`` descendente, reiniciando numeración en 1.

    TODO:
    - Implementar estrategia basada en conectividad (BFS/DFS) cuando se
      necesite una numeración más "química".
    """

    if strategy != "xy":
        raise NotImplementedError(f"Atom numbering strategy not implemented: {strategy}")

    atom_numbers: dict[int, int] = {}
    for component in _connected_components(model):
        atoms = sorted(
            (model.atoms[atom_id] for atom_id in component),
            key=lambda atom: (float(atom.x), -float(atom.y), int(atom.id)),
        )
        for idx, atom in enumerate(atoms, start=1):
            atom_numbers[int(atom.id)] = idx
    return atom_numbers


def compute_structure_numbers(model: MolGraph) -> List[NumberedStructure]:
    """Numeración de estructuras por componentes conexas del grafo."""

    sortable: list[tuple[float, float, int, tuple[int, ...]]] = []
    for component in _connected_components(model):
        atoms = [model.atoms[atom_id] for atom_id in component]
        center_x = sum(float(atom.x) for atom in atoms) / len(atoms)
        center_y = sum(float(atom.y) for atom in atoms) / len(atoms)
        atom_ids_tuple = tuple(component)
        structure_id = min(atom_ids_tuple)
        sortable.append((center_x, -center_y, structure_id, atom_ids_tuple))

    sortable.sort(key=lambda entry: (entry[0], entry[1], entry[2]))
    result: list[NumberedStructure] = []
    for idx, (center_x, neg_center_y, structure_id, atom_ids_tuple) in enumerate(sortable, start=1):
        result.append(
            NumberedStructure(
                structure_id=structure_id,
                number=idx,
                x=center_x,
                y=-neg_center_y,
                atom_ids=atom_ids_tuple,
            )
        )
    return result


def _connected_components(model: MolGraph) -> list[list[int]]:
    """Devuelve componentes conexas sobre los enlaces del modelo."""
    atom_ids = sorted(int(atom_id) for atom_id in model.atoms.keys())
    if not atom_ids:
        return []

    adjacency: dict[int, set[int]] = {atom_id: set() for atom_id in atom_ids}
    for bond in model.bonds.values():
        a1 = int(bond.a1_id)
        a2 = int(bond.a2_id)
        if a1 not in adjacency or a2 not in adjacency:
            continue
        adjacency[a1].add(a2)
        adjacency[a2].add(a1)

    components: list[list[int]] = []
    visited: set[int] = set()
    for start in atom_ids:
        if start in visited:
            continue
        queue = deque([start])
        component: list[int] = []
        while queue:
            atom_id = queue.popleft()
            if atom_id in visited:
                continue
            visited.add(atom_id)
            component.append(atom_id)
            for neighbor in sorted(adjacency.get(atom_id, set())):
                if neighbor not in visited:
                    queue.append(neighbor)
        components.append(sorted(component))
    return components
