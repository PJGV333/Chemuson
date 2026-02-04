import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import MolGraph
from chemname import iupac_name


def add_aromatic_ring(graph: MolGraph, atoms: list, ring: list[int]) -> None:
    for i in range(len(ring)):
        a1 = atoms[ring[i]].id
        a2 = atoms[ring[(i + 1) % len(ring)]].id
        if graph.find_bond_between(a1, a2) is None:
            graph.add_bond(a1, a2, order=1, is_aromatic=True)


def build_quinoline(graph: MolGraph, n_index: int) -> list[int]:
    elements = ["C"] * 10
    elements[n_index] = "N"
    atoms = [graph.add_atom(elem, float(i), 0.0) for i, elem in enumerate(elements)]
    ring1 = [0, 1, 2, 3, 4, 5]
    ring2 = [3, 4, 6, 7, 8, 9]
    add_aromatic_ring(graph, atoms, ring1)
    add_aromatic_ring(graph, atoms, ring2)
    return [atom.id for atom in atoms]


def build_indole(graph: MolGraph) -> list[int]:
    elements = ["C"] * 9
    elements[6] = "N"
    atoms = [graph.add_atom(elem, float(i), 0.0) for i, elem in enumerate(elements)]
    ring1 = [0, 1, 2, 3, 4, 5]
    ring2 = [3, 4, 6, 7, 8]
    add_aromatic_ring(graph, atoms, ring1)
    add_aromatic_ring(graph, atoms, ring2)
    return [atom.id for atom in atoms]


class ChemNamePR21Test(unittest.TestCase):
    def test_quinoline(self):
        graph = MolGraph()
        build_quinoline(graph, n_index=6)
        self.assertEqual(iupac_name(graph), "quinoline")

    def test_chloroquinoline(self):
        graph = MolGraph()
        atoms = build_quinoline(graph, n_index=6)
        cl = graph.add_atom("Cl", -1.0, 0.0)
        graph.add_bond(atoms[7], cl.id, order=1)
        self.assertEqual(iupac_name(graph), "2-chloroquinoline")

    def test_isoquinoline(self):
        graph = MolGraph()
        build_quinoline(graph, n_index=7)
        self.assertEqual(iupac_name(graph), "isoquinoline")

    def test_chloroisoquinoline(self):
        graph = MolGraph()
        atoms = build_quinoline(graph, n_index=7)
        cl = graph.add_atom("Cl", -1.0, 0.0)
        graph.add_bond(atoms[6], cl.id, order=1)
        self.assertEqual(iupac_name(graph), "2-chloroisoquinoline")

    def test_indole(self):
        graph = MolGraph()
        build_indole(graph)
        self.assertEqual(iupac_name(graph), "indole")

    def test_methylindole(self):
        graph = MolGraph()
        atoms = build_indole(graph)
        methyl = graph.add_atom("C", -1.0, 0.0)
        graph.add_bond(atoms[7], methyl.id, order=1)
        self.assertEqual(iupac_name(graph), "2-methylindole")


if __name__ == "__main__":
    unittest.main()
