"""ChemName template matching behavior."""

from chemuson.chemname.molview import MolView
from chemuson.chemname.template import load_template
from chemuson.chemname.template_match import (
    _mapping_substituent_locants,
    match_template_exact,
    select_template_mapping,
)
from chemuson.core.model import MolGraph
from chemuson.utils.resources import open_resource_path


def _build_benzene(graph: MolGraph) -> list[int]:
    atoms = [graph.add_atom("C", float(index), 0.0) for index in range(6)]
    for index in range(6):
        a1 = atoms[index].id
        a2 = atoms[(index + 1) % 6].id
        if graph.find_bond_between(a1, a2) is None:
            graph.add_bond(a1, a2, order=1, is_aromatic=True)
    return [atom.id for atom in atoms]


def test_benzene_template_matches_aromatic_ring() -> None:
    with open_resource_path("chemname", "templates", "simple", "benzene.mol") as template_path:
        template = load_template(template_path)
    graph = MolGraph()
    ring = _build_benzene(graph)

    mappings = match_template_exact(template, MolView(graph), atom_ids=ring)

    assert len(mappings) >= 1


def test_template_selector_is_deterministic_and_prefers_lowest_locant() -> None:
    with open_resource_path("chemname", "templates", "simple", "benzene.mol") as template_path:
        template = load_template(template_path)
    graph = MolGraph()
    ring = _build_benzene(graph)
    chlorine = graph.add_atom("Cl", -1.0, 0.0)
    graph.add_bond(ring[2], chlorine.id, order=1)
    view = MolView(graph)
    mappings = match_template_exact(template, view, atom_ids=ring)

    chosen = select_template_mapping(template, view, mappings)

    assert chosen is not None
    assert sorted(_mapping_substituent_locants(template, view, chosen)) == [1]
    assert select_template_mapping(template, view, mappings) == chosen
