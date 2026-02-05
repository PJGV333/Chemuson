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


def build_benzene_kekule(graph: MolGraph) -> list[int]:
    atoms = [graph.add_atom("C", float(i), 0.0) for i in range(6)]
    for i in range(6):
        order = 2 if i % 2 == 0 else 1
        graph.add_bond(atoms[i].id, atoms[(i + 1) % 6].id, order=order)
    return [atom.id for atom in atoms]


class ChemNamePR24Test(unittest.TestCase):
    def test_methoxypropane(self):
        graph = MolGraph()
        chain = build_linear_chain(graph, 3)
        o = graph.add_atom("O", 1.0, 1.0)
        me = graph.add_atom("C", 1.0, 2.0)
        graph.add_bond(chain[1], o.id, order=1)
        graph.add_bond(o.id, me.id, order=1)
        self.assertEqual(iupac_name(graph), "2-methoxypropane")

    def test_chloro_methoxypropane(self):
        graph = MolGraph()
        chain = build_linear_chain(graph, 3)
        o = graph.add_atom("O", -1.0, 0.0)
        me = graph.add_atom("C", -2.0, 0.0)
        cl = graph.add_atom("Cl", 3.0, 0.0)
        graph.add_bond(chain[0], o.id, order=1)
        graph.add_bond(o.id, me.id, order=1)
        graph.add_bond(chain[2], cl.id, order=1)
        self.assertEqual(iupac_name(graph), "1-chloro-3-methoxypropane")

    def test_nitroethane(self):
        graph = MolGraph()
        chain = build_linear_chain(graph, 2)
        n = graph.add_atom("N", -1.0, 0.0)
        o1 = graph.add_atom("O", -2.0, 0.0)
        o2 = graph.add_atom("O", -2.0, 1.0)
        graph.add_bond(chain[0], n.id, order=1)
        graph.add_bond(n.id, o1.id, order=2)
        graph.add_bond(n.id, o2.id, order=1)
        self.assertEqual(iupac_name(graph), "1-nitroethane")

    def test_chloronitrobenzene(self):
        graph = MolGraph()
        ring = build_benzene_kekule(graph)
        n = graph.add_atom("N", -1.0, 0.0)
        o1 = graph.add_atom("O", -2.0, 0.0)
        o2 = graph.add_atom("O", -2.0, 1.0)
        cl = graph.add_atom("Cl", 3.0, 1.0)
        graph.add_bond(ring[0], n.id, order=1)
        graph.add_bond(n.id, o1.id, order=2)
        graph.add_bond(n.id, o2.id, order=1)
        graph.add_bond(ring[3], cl.id, order=1)
        self.assertEqual(iupac_name(graph), "1-chloro-4-nitrobenzene")

    def test_aminobenzene(self):
        graph = MolGraph()
        ring = build_benzene_kekule(graph)
        n = graph.add_atom("N", -1.0, 0.0)
        graph.add_bond(ring[0], n.id, order=1)
        self.assertEqual(iupac_name(graph), "aminobenzene")


if __name__ == "__main__":
    unittest.main()
