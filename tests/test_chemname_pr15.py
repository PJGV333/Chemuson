import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import MolGraph
from chemname import iupac_name


def build_naphthalene(graph: MolGraph) -> list[int]:
    atoms = [graph.add_atom("C", float(i), 0.0) for i in range(10)]
    ring1 = [0, 1, 2, 3, 4, 5]
    ring2 = [3, 4, 6, 7, 8, 9]

    for ring in (ring1, ring2):
        for i in range(len(ring)):
            a1 = atoms[ring[i]].id
            a2 = atoms[ring[(i + 1) % len(ring)]].id
            if graph.find_bond_between(a1, a2) is None:
                graph.add_bond(a1, a2, order=1, is_aromatic=True)
    return [atom.id for atom in atoms]


class ChemNamePR15Test(unittest.TestCase):
    def test_naphthalene(self):
        graph = MolGraph()
        build_naphthalene(graph)
        self.assertEqual(iupac_name(graph), "naphthalene")

    def test_1_chloronaphthalene(self):
        graph = MolGraph()
        atoms = build_naphthalene(graph)
        cl = graph.add_atom("Cl", -1.0, 0.0)
        graph.add_bond(atoms[2], cl.id, order=1)
        self.assertEqual(iupac_name(graph), "1-chloronaphthalene")

    def test_2_chloronaphthalene(self):
        graph = MolGraph()
        atoms = build_naphthalene(graph)
        cl = graph.add_atom("Cl", -1.0, 0.0)
        graph.add_bond(atoms[0], cl.id, order=1)
        self.assertEqual(iupac_name(graph), "2-chloronaphthalene")

    def test_1_methylnaphthalene(self):
        graph = MolGraph()
        atoms = build_naphthalene(graph)
        me = graph.add_atom("C", -1.0, 0.0)
        graph.add_bond(atoms[5], me.id, order=1)
        self.assertEqual(iupac_name(graph), "1-methylnaphthalene")


if __name__ == "__main__":
    unittest.main()
