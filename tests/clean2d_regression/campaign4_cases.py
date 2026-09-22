from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

from chemuson.core.model import MolGraph


@dataclass(frozen=True)
class Campaign4Case:
    name: str
    family: str
    builder: Callable[[], MolGraph]
    tags: tuple[str, ...]


def _regular_ring(size: int, *, aromatic: bool, atom_start: int = 1, bond_start: int = 1, radius: float = 42.0) -> MolGraph:
    graph = MolGraph()
    for index in range(size):
        angle = math.radians(index * 360.0 / size - 90.0)
        graph.add_atom(
            "C",
            math.cos(angle) * radius,
            math.sin(angle) * radius,
            atom_id=atom_start + index,
        )
    for index in range(size):
        graph.add_bond(
            atom_start + index,
            atom_start + ((index + 1) % size),
            bond_id=bond_start + index,
            is_aromatic=aromatic,
        )
    return graph


def _fused_bicyclic(*, aromatic: bool) -> MolGraph:
    graph = MolGraph()
    coords = {
        1: (0.0, 24.0), 2: (36.0, 0.0), 3: (72.0, 24.0), 4: (72.0, 66.0),
        5: (36.0, 90.0), 6: (0.0, 66.0), 7: (108.0, 0.0), 8: (144.0, 24.0),
        9: (144.0, 66.0), 10: (108.0, 90.0),
    }
    for atom_id, (x, y) in coords.items():
        graph.add_atom("C", x, y, atom_id=atom_id)
    bonds = ((1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1), (3, 7), (7, 8), (8, 9), (9, 10), (10, 4))
    for bond_id, (left, right) in enumerate(bonds, start=1):
        graph.add_bond(left, right, bond_id=bond_id, is_aromatic=aromatic)
    return graph


def _three_fused_ring_polycyclic() -> MolGraph:
    graph = MolGraph()
    coords = {
        1: (0.0, 24.0), 2: (36.0, 0.0), 3: (72.0, 24.0), 4: (72.0, 66.0),
        5: (36.0, 90.0), 6: (0.0, 66.0), 7: (108.0, 0.0), 8: (144.0, 24.0),
        9: (144.0, 66.0), 10: (108.0, 90.0), 11: (180.0, 0.0), 12: (216.0, 24.0),
        13: (216.0, 66.0), 14: (180.0, 90.0),
    }
    for atom_id, (x, y) in coords.items():
        graph.add_atom("C", x, y, atom_id=atom_id)
    bonds = (
        (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1),
        (3, 7), (7, 8), (8, 9), (9, 10), (10, 4),
        (9, 11), (11, 12), (12, 13), (13, 14), (14, 10),
    )
    for bond_id, (left, right) in enumerate(bonds, start=1):
        graph.add_bond(left, right, bond_id=bond_id, is_aromatic=True)
    return graph


def _spiro_bicyclic_graph() -> MolGraph:
    graph = MolGraph()
    graph.add_atom("C", 0.0, 0.0, atom_id=1)
    next_atom = 2
    next_bond = 1
    for center_angle in (0.0, 180.0):
        ring_ids = list(range(next_atom, next_atom + 5))
        next_atom += 5
        for index, atom_id in enumerate(ring_ids):
            angle = math.radians(center_angle + 30.0 + index * 72.0)
            graph.add_atom("C", math.cos(angle) * 42.0, math.sin(angle) * 42.0, atom_id=atom_id)
        path = [1, *ring_ids]
        for left, right in zip(path, path[1:]):
            graph.add_bond(left, right, bond_id=next_bond)
            next_bond += 1
        graph.add_bond(ring_ids[-1], 1, bond_id=next_bond)
        next_bond += 1
    return graph


def _bridged_bicyclic_graph() -> MolGraph:
    graph = MolGraph()
    coords = {1: (-42.0, 0.0), 2: (-12.0, 42.0), 3: (12.0, 42.0), 4: (42.0, 0.0), 5: (0.0, -42.0), 6: (0.0, 4.0)}
    for atom_id, (x, y) in coords.items():
        graph.add_atom("C", x, y, atom_id=atom_id)
    for bond_id, pair in enumerate(((1, 2), (2, 3), (3, 4), (1, 5), (5, 4), (1, 6), (6, 4)), start=1):
        graph.add_bond(*pair, bond_id=bond_id)
    return graph


def _polycyclic_graph() -> MolGraph:
    return _fused_bicyclic(aromatic=True)


def _fused_one_substituent() -> MolGraph:
    graph = _fused_bicyclic(aromatic=True)
    graph.add_atom("O", 36.0, 132.0, atom_id=11)
    graph.add_bond(5, 11, bond_id=12)
    return graph


def _fused_multiple_substituents() -> MolGraph:
    graph = _fused_bicyclic(aromatic=True)
    graph.add_atom("O", 36.0, 132.0, atom_id=11)
    graph.add_atom("N", 180.0, 24.0, atom_id=12)
    graph.add_bond(5, 11, bond_id=12)
    graph.add_bond(8, 12, bond_id=13)
    return graph


def _spiro_substituent() -> MolGraph:
    graph = _spiro_bicyclic_graph()
    graph.add_atom("O", 0.0, 84.0, atom_id=12)
    graph.add_bond(3, 12, bond_id=13)
    return graph


def _two_rigid_blocks_linker() -> MolGraph:
    left = _regular_ring(6, aromatic=True)
    right = _regular_ring(6, aromatic=True, atom_start=7, bond_start=7)
    for atom_id, atom in right.atoms.items():
        atom.x += 118.0
    for atom_id, atom in right.atoms.items():
        left.atoms[atom_id] = atom
    for bond_id, bond in right.bonds.items():
        left.bonds[bond_id] = bond
    left.add_atom("C", 118.0, 0.0, atom_id=13)
    left.add_atom("C", 160.0, 0.0, atom_id=14)
    left.add_bond(1, 13, bond_id=13)
    left.add_bond(13, 14, bond_id=14)
    left.add_bond(14, 7, bond_id=15)
    return left


def _congested_attachments() -> MolGraph:
    graph = _regular_ring(6, aromatic=True)
    graph.add_atom("C", -42.0, -30.0, atom_id=7)
    graph.add_atom("O", 42.0, -30.0, atom_id=8)
    graph.add_atom("N", 0.0, 78.0, atom_id=9)
    graph.add_bond(1, 7, bond_id=7)
    graph.add_bond(2, 8, bond_id=8)
    graph.add_bond(3, 9, bond_id=9)
    return graph


CAMPAIGN4_CASES: tuple[Campaign4Case, ...] = (
    Campaign4Case("regular_benzene_control", "monocycle", lambda: _regular_ring(6, aromatic=True), ("control", "aromatic")),
    Campaign4Case("simple_monocycle_control", "monocycle", lambda: _regular_ring(6, aromatic=False), ("control",)),
    Campaign4Case("fused_aromatic_bicyclic", "fused", lambda: _fused_bicyclic(aromatic=True), ("fused", "aromatic")),
    Campaign4Case("fused_non_aromatic_bicyclic", "fused", lambda: _fused_bicyclic(aromatic=False), ("fused",)),
    Campaign4Case("spiro_bicyclic", "spiro", _spiro_bicyclic_graph, ("spiro",)),
    Campaign4Case("bridged_bicyclic", "bridged", _bridged_bicyclic_graph, ("bridged",)),
    Campaign4Case("larger_polycyclic", "polycyclic", _polycyclic_graph, ("polycyclic", "fused")),
    Campaign4Case("three_fused_ring_polycyclic", "polycyclic", _three_fused_ring_polycyclic, ("polycyclic", "fused")),
    Campaign4Case("fused_one_substituent", "fused-substitution", _fused_one_substituent, ("fused", "congested")),
    Campaign4Case("fused_multiple_substituents", "fused-substitution", _fused_multiple_substituents, ("fused", "congested")),
    Campaign4Case("spiro_substituent", "spiro-substitution", _spiro_substituent, ("spiro", "congested")),
    Campaign4Case("two_rigid_blocks_linker", "multiple-rigid-blocks", _two_rigid_blocks_linker, ("multiblock",)),
    Campaign4Case("congested_ring_attachments", "congested", _congested_attachments, ("congested",)),
)

CASES_BY_NAME = {case.name: case for case in CAMPAIGN4_CASES}
