import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import MolGraph
from chemname import iupac_name


def build_anthracene(graph: MolGraph) -> list[int]:
    atoms = [graph.add_atom("C", float(i), 0.0) for i in range(14)]
    ring1 = [0, 1, 2, 3, 4, 5]
    ring2 = [3, 4, 6, 7, 8, 9]
    ring3 = [8, 9, 10, 11, 12, 13]
    for ring in (ring1, ring2, ring3):
        for i in range(len(ring)):
            a1 = atoms[ring[i]].id
            a2 = atoms[ring[(i + 1) % len(ring)]].id
            if graph.find_bond_between(a1, a2) is None:
                graph.add_bond(a1, a2, order=1, is_aromatic=True)
    return [atom.id for atom in atoms]


def build_phenanthrene(graph: MolGraph) -> list[int]:
    atoms = [graph.add_atom("C", float(i), 0.0) for i in range(14)]
    ring1 = [0, 1, 2, 3, 4, 5]
    ring2 = [3, 4, 6, 7, 8, 9]
    ring3 = [4, 6, 10, 11, 12, 13]
    for ring in (ring1, ring2, ring3):
        for i in range(len(ring)):
            a1 = atoms[ring[i]].id
            a2 = atoms[ring[(i + 1) % len(ring)]].id
            if graph.find_bond_between(a1, a2) is None:
                graph.add_bond(a1, a2, order=1, is_aromatic=True)
    return [atom.id for atom in atoms]


class ChemNamePR20Test(unittest.TestCase):
    def test_anthracene(self):
        graph = MolGraph()
        build_anthracene(graph)
        self.assertEqual(iupac_name(graph), "anthracene")

    def test_chloroanthracene(self):
        graph = MolGraph()
        atoms = build_anthracene(graph)
        cl = graph.add_atom("Cl", -1.0, 0.0)
        graph.add_bond(atoms[2], cl.id, order=1)
        self.assertEqual(iupac_name(graph), "1-chloroanthracene")

    def test_phenanthrene(self):
        graph = MolGraph()
        build_phenanthrene(graph)
        self.assertEqual(iupac_name(graph), "phenanthrene")

    def test_chlorophenanthrene(self):
        graph = MolGraph()
        atoms = build_phenanthrene(graph)
        cl = graph.add_atom("Cl", -1.0, 0.0)
        graph.add_bond(atoms[2], cl.id, order=1)
        self.assertEqual(iupac_name(graph), "1-chlorophenanthrene")


if __name__ == "__main__":
    unittest.main()
