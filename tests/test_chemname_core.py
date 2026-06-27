"""Core ChemName and MolView behavior."""

from chemuson.chemname import iupac_name
from chemuson.chemname.molview import MolView
from chemuson.core.model import MolGraph


def _build_linear_chain(graph: MolGraph, length: int) -> list[int]:
    ids = []
    previous_id = None
    for index in range(length):
        atom = graph.add_atom("C", float(index), 0.0)
        ids.append(atom.id)
        if previous_id is not None:
            graph.add_bond(previous_id, atom.id, order=1)
        previous_id = atom.id
    return ids


def test_iupac_name_returns_not_available_for_empty_input() -> None:
    assert iupac_name(None) == "N/D"


def test_molview_detects_acyclic_chain() -> None:
    graph = MolGraph()
    _build_linear_chain(graph, 3)

    assert MolView(graph).is_acyclic()


def test_molview_detects_cycle() -> None:
    graph = MolGraph()
    a1 = graph.add_atom("C", 0.0, 0.0)
    a2 = graph.add_atom("C", 1.0, 0.0)
    a3 = graph.add_atom("C", 0.5, 1.0)
    graph.add_bond(a1.id, a2.id, order=1)
    graph.add_bond(a2.id, a3.id, order=1)
    graph.add_bond(a3.id, a1.id, order=1)

    assert not MolView(graph).is_acyclic()
