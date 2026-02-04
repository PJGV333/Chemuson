import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import MolGraph
from chemname import iupac_name


def build_ring(graph: MolGraph, size: int) -> list[int]:
    atoms = [graph.add_atom("C", float(i), 0.0) for i in range(size)]
    for i in range(size):
        graph.add_bond(atoms[i].id, atoms[(i + 1) % size].id, order=1)
    return [atom.id for atom in atoms]


class ChemNamePR10Test(unittest.TestCase):
    def test_dimethylcyclohexane(self):
        graph = MolGraph()
        ring = build_ring(graph, 6)
        m1 = graph.add_atom("C", -1.0, 0.0)
        m2 = graph.add_atom("C", 3.0, 1.0)
        graph.add_bond(ring[0], m1.id, order=1)
        graph.add_bond(ring[3], m2.id, order=1)
        self.assertEqual(iupac_name(graph), "1,4-dimethylcyclohexane")

    def test_bromochlorocyclohexane(self):
        graph = MolGraph()
        ring = build_ring(graph, 6)
        br = graph.add_atom("Br", -1.0, 0.0)
        cl = graph.add_atom("Cl", 3.0, 1.0)
        graph.add_bond(ring[0], br.id, order=1)
        graph.add_bond(ring[3], cl.id, order=1)
        self.assertEqual(iupac_name(graph), "1-bromo-4-chlorocyclohexane")


if __name__ == "__main__":
    unittest.main()
