"""Molecular formula and mass calculations."""

from chemuson.chemcalc import format_formula, molecular_formula, molecular_weight
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


def test_methane_formula_from_single_carbon() -> None:
    graph = MolGraph()
    graph.add_atom("C", 0.0, 0.0)

    assert format_formula(molecular_formula(graph)) == "CH4"


def test_chloromethane_formula() -> None:
    graph = MolGraph()
    carbon = graph.add_atom("C", 0.0, 0.0)
    chlorine = graph.add_atom("Cl", 1.0, 0.0)
    graph.add_bond(carbon.id, chlorine.id, order=1)

    assert format_formula(molecular_formula(graph)) == "CH3Cl"


def test_bromo_chloro_dodecane_formula_and_weight() -> None:
    graph = MolGraph()
    chain_ids = _build_linear_chain(graph, 12)
    chlorine = graph.add_atom("Cl", -1.0, 0.0)
    graph.add_bond(chain_ids[0], chlorine.id, order=1)
    bromine = graph.add_atom("Br", 13.0, 0.0)
    graph.add_bond(chain_ids[-1], bromine.id, order=1)

    formula = molecular_formula(graph)

    assert format_formula(formula) == "C12H24BrCl"
    assert molecular_weight(formula) > 0.0
