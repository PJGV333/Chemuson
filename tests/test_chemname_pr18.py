import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import MolGraph
from chemname import iupac_name


def build_benzene_kekule(graph: MolGraph) -> list[int]:
    atoms = [graph.add_atom("C", float(i), 0.0) for i in range(6)]
    for i in range(6):
        order = 2 if i % 2 == 0 else 1
        graph.add_bond(atoms[i].id, atoms[(i + 1) % 6].id, order=order)
    return [atom.id for atom in atoms]


def build_cyclohexane(graph: MolGraph) -> list[int]:
    atoms = [graph.add_atom("C", float(i), 0.0) for i in range(6)]
    for i in range(6):
        graph.add_bond(atoms[i].id, atoms[(i + 1) % 6].id, order=1)
    return [atom.id for atom in atoms]


def build_linear_chain(graph: MolGraph, length: int) -> list[int]:
    ids = []
    prev = None
    for i in range(length):
        atom = graph.add_atom("C", float(i), 0.0)
        ids.append(atom.id)
        if prev is not None:
            graph.add_bond(prev, atom.id, order=1)
        prev = atom.id
    return ids


class ChemNamePR18Test(unittest.TestCase):
    def test_ethylbenzene_parent(self):
        graph = MolGraph()
        ring = build_benzene_kekule(graph)
        ethyl = graph.add_atom("C", -1.0, 0.0)
        ethyl2 = graph.add_atom("C", -2.0, 0.0)
        graph.add_bond(ring[0], ethyl.id, order=1)
        graph.add_bond(ethyl.id, ethyl2.id, order=1)
        self.assertEqual(iupac_name(graph), "ethylbenzene")

    def test_propylcyclohexane_parent(self):
        graph = MolGraph()
        ring = build_cyclohexane(graph)
        p1 = graph.add_atom("C", -1.0, 0.0)
        p2 = graph.add_atom("C", -2.0, 0.0)
        p3 = graph.add_atom("C", -3.0, 0.0)
        graph.add_bond(ring[0], p1.id, order=1)
        graph.add_bond(p1.id, p2.id, order=1)
        graph.add_bond(p2.id, p3.id, order=1)
        self.assertEqual(iupac_name(graph), "propylcyclohexane")

    def test_phenethyl_alcohol(self):
        graph = MolGraph()
        chain = build_linear_chain(graph, 2)
        o = graph.add_atom("O", 2.0, 0.0)
        ring = build_benzene_kekule(graph)
        graph.add_bond(chain[0], o.id, order=1)
        graph.add_bond(chain[1], ring[0], order=1)
        self.assertEqual(iupac_name(graph), "2-phenylethan-1-ol")


if __name__ == "__main__":
    unittest.main()
