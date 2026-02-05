import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import MolGraph
from chemname.molview import MolView
from chemname.template import load_template
from chemname.template_match import match_template_exact, select_template_mapping, _mapping_substituent_locants


def build_benzene(graph: MolGraph) -> list[int]:
    atoms = [graph.add_atom("C", float(i), 0.0) for i in range(6)]
    for i in range(6):
        a1 = atoms[i].id
        a2 = atoms[(i + 1) % 6].id
        if graph.find_bond_between(a1, a2) is None:
            graph.add_bond(a1, a2, order=1, is_aromatic=True)
    return [atom.id for atom in atoms]


class TemplateMatchPRATest(unittest.TestCase):
    def test_benzene_template_match(self):
        template_path = os.path.join(
            os.path.dirname(__file__),
            "..",
            "src",
            "chemname",
            "templates",
            "simple",
            "benzene.mol",
        )
        template = load_template(template_path)
        graph = MolGraph()
        ring = build_benzene(graph)
        view = MolView(graph)
        mappings = match_template_exact(template, view, atom_ids=ring)
        self.assertGreaterEqual(len(mappings), 1)

    def test_selector_lowest_locant(self):
        template_path = os.path.join(
            os.path.dirname(__file__),
            "..",
            "src",
            "chemname",
            "templates",
            "simple",
            "benzene.mol",
        )
        template = load_template(template_path)
        graph = MolGraph()
        ring = build_benzene(graph)
        cl = graph.add_atom("Cl", -1.0, 0.0)
        graph.add_bond(ring[2], cl.id, order=1)
        view = MolView(graph)
        mappings = match_template_exact(template, view, atom_ids=ring)
        chosen = select_template_mapping(template, view, mappings)
        self.assertIsNotNone(chosen)
        locants = _mapping_substituent_locants(template, view, chosen)
        self.assertEqual(sorted(locants), [1])
        chosen_again = select_template_mapping(template, view, mappings)
        self.assertEqual(chosen, chosen_again)


if __name__ == "__main__":
    unittest.main()
