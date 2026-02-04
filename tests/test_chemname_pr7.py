import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import MolGraph
from chemname import iupac_name


def build_linear_chain(graph: MolGraph, length: int) -> list[int]:
    ids = []
    prev_id = None
    for i in range(length):
        atom = graph.add_atom("C", float(i), 0.0)
        ids.append(atom.id)
        if prev_id is not None:
            graph.add_bond(prev_id, atom.id, order=1)
        prev_id = atom.id
    return ids


class ChemNamePR7Test(unittest.TestCase):
    def test_propene_orientation(self):
        graph = MolGraph()
        c1 = graph.add_atom("C", 0.0, 0.0)
        c2 = graph.add_atom("C", 1.0, 0.0)
        c3 = graph.add_atom("C", 2.0, 0.0)
        graph.add_bond(c1.id, c2.id, order=2)
        graph.add_bond(c2.id, c3.id, order=1)
        self.assertEqual(iupac_name(graph), "prop-1-ene")

    def test_propene_reverse_orientation(self):
        graph = MolGraph()
        c1 = graph.add_atom("C", 0.0, 0.0)
        c2 = graph.add_atom("C", 1.0, 0.0)
        c3 = graph.add_atom("C", 2.0, 0.0)
        graph.add_bond(c1.id, c2.id, order=1)
        graph.add_bond(c2.id, c3.id, order=2)
        self.assertEqual(iupac_name(graph), "prop-1-ene")

    def test_bromopropene(self):
        graph = MolGraph()
        c1 = graph.add_atom("C", 0.0, 0.0)
        c2 = graph.add_atom("C", 1.0, 0.0)
        c3 = graph.add_atom("C", 2.0, 0.0)
        graph.add_bond(c1.id, c2.id, order=1)
        graph.add_bond(c2.id, c3.id, order=2)
        br = graph.add_atom("Br", -1.0, 0.0)
        graph.add_bond(c1.id, br.id, order=1)
        self.assertEqual(iupac_name(graph), "3-bromoprop-1-ene")

    def test_propan_1_ol(self):
        graph = MolGraph()
        chain = build_linear_chain(graph, 3)
        o = graph.add_atom("O", 3.0, 0.0)
        graph.add_bond(chain[-1], o.id, order=1)
        self.assertEqual(iupac_name(graph), "propan-1-ol")

    def test_propan_2_ol(self):
        graph = MolGraph()
        chain = build_linear_chain(graph, 3)
        o = graph.add_atom("O", 1.0, 1.0)
        graph.add_bond(chain[1], o.id, order=1)
        self.assertEqual(iupac_name(graph), "propan-2-ol")

    def test_propan_1_amine(self):
        graph = MolGraph()
        chain = build_linear_chain(graph, 3)
        n = graph.add_atom("N", 3.0, 0.0)
        graph.add_bond(chain[-1], n.id, order=1)
        self.assertEqual(iupac_name(graph), "propan-1-amine")


if __name__ == "__main__":
    unittest.main()
