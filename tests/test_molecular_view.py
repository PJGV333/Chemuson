"""Caracterización de la vista molecular compartida."""

from __future__ import annotations

from dataclasses import dataclass

import pytest

from chemuson.chemcalc import implicit_h_count, molecular_formula
from chemuson.chemname.errors import ChemNameNotSupported
from chemuson.chemname.molview import MolView
from chemuson.core import MolecularViewNotSupported
from chemuson.core.model import MolGraph


@dataclass
class _Atom:
    id: int
    element: str
    charge: int = 0
    isotope: int | None = None
    radical_electrons: int = 0
    explicit_h: int | None = None
    no_implicit: bool = False


class _Graph:
    def __init__(self, atoms, bonds=()) -> None:
        self.atoms = atoms
        self.bonds = bonds


def test_molgraph_queries_and_formula_do_not_mutate_graph() -> None:
    graph = MolGraph()
    carbon = graph.add_atom("C", 0.0, 0.0)
    oxygen = graph.add_atom("O", 1.0, 0.0, isotope=18, radical_electrons=1)
    graph.add_bond(carbon.id, oxygen.id, order=2, is_aromatic=True)
    atoms_before = dict(graph.atoms)
    bonds_before = dict(graph.bonds)

    view = MolView(graph)

    assert view.atoms() == [carbon.id, oxygen.id]
    assert view.neighbors(carbon.id) == [oxygen.id]
    assert view.bond_order_between(carbon.id, oxygen.id) == 2
    assert view.bond_is_aromatic(carbon.id, oxygen.id) is True
    assert view.isotope(oxygen.id) == 18
    assert view.has_radical(oxygen.id) is True
    assert molecular_formula(graph) == {"C": 1, "O": 1}
    assert graph.atoms == atoms_before
    assert graph.bonds == bonds_before


def test_mapping_atoms_and_tuple_bonds_are_supported() -> None:
    graph = _Graph(
        {
            1: {"element": "C", "explicit_h": 1},
            2: {"symbol": "N", "formal_charge": 1},
        },
        [(1, 2, 1)],
    )

    view = MolView(graph)

    assert view.atoms() == [1, 2]
    assert view.element(2) == "N"
    assert view.explicit_h(1) == 1
    assert view.atom_charge(2) == 1
    assert view.bonds() == [(1, 2, 1)]
    assert molecular_formula(graph) == {"C": 1, "N": 1, "H": 5}


def test_iterable_atoms_and_metadata_are_supported() -> None:
    first = _Atom(1, "C", explicit_h=2)
    second = _Atom(2, "O", charge=-1)
    graph = _Graph([first, second], [(1, 2, 1)])
    view = MolView(graph)

    assert view.atoms() == [1, 2]
    assert view.formal_charge(2) == -1
    assert view.explicit_h(1) == 2
    assert view.radical_electrons(1) == 0
    assert implicit_h_count(view, 2) == 0


def test_nodes_edges_and_neighbor_only_graphs_are_supported() -> None:
    class NodeGraph:
        nodes = {1: {"element": "C"}, 2: {"element": "C"}}
        edges = [(1, 2, {"order": 2, "is_aromatic": True})]

    class NeighborGraph:
        atoms = {1: {"element": "C"}, 2: {"element": "C"}}
        bonds = ()

        @staticmethod
        def neighbors(atom_id: int) -> list[int]:
            return [2] if atom_id == 1 else [1]

    edge_view = MolView(NodeGraph())
    neighbor_view = MolView(NeighborGraph())

    assert edge_view.bond_order_between(1, 2) == 2
    assert edge_view.bond_is_aromatic(1, 2) is True
    assert neighbor_view.neighbors(1) == [2]
    assert neighbor_view.bond_order_between(1, 2) == 1


def test_acyclicity_and_structural_bonds_are_preserved() -> None:
    graph = MolGraph()
    first = graph.add_atom("C", 0.0, 0.0)
    second = graph.add_atom("C", 1.0, 0.0)
    third = graph.add_atom("C", 0.5, 1.0)
    graph.add_bond(first.id, second.id)
    graph.add_bond(second.id, third.id)
    graph.add_bond(third.id, first.id)

    assert MolView(graph).is_acyclic() is False


def test_historical_errors_are_chemname_not_supported() -> None:
    class EmptyGraph:
        pass

    class MissingElementGraph:
        atoms = {1: {}}
        bonds = ()

    with pytest.raises(ChemNameNotSupported, match="Cannot read atoms from graph"):
        MolView(EmptyGraph()).atoms()
    with pytest.raises(ChemNameNotSupported, match="Cannot resolve atom element"):
        MolView(MissingElementGraph()).element(1)


def test_historical_error_is_the_neutral_core_error() -> None:
    assert ChemNameNotSupported is MolecularViewNotSupported
