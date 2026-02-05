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


def build_ring(graph: MolGraph, size: int) -> list[int]:
    atoms = [graph.add_atom("C", float(i), 0.0) for i in range(size)]
    for i in range(size):
        a1 = atoms[i].id
        a2 = atoms[(i + 1) % size].id
        graph.add_bond(a1, a2, order=1)
    return [atom.id for atom in atoms]


class ChemNamePR22Test(unittest.TestCase):
    def test_butadiene(self):
        graph = MolGraph()
        chain = build_linear_chain(graph, 4)
        graph.find_bond_between(chain[0], chain[1]).order = 2
        graph.find_bond_between(chain[2], chain[3]).order = 2
        self.assertEqual(iupac_name(graph), "buta-1,3-diene")

    def test_hexadiyne(self):
        graph = MolGraph()
        chain = build_linear_chain(graph, 6)
        graph.find_bond_between(chain[0], chain[1]).order = 3
        graph.find_bond_between(chain[4], chain[5]).order = 3
        self.assertEqual(iupac_name(graph), "hexa-1,5-diyne")

    def test_hexatriene(self):
        graph = MolGraph()
        chain = build_linear_chain(graph, 6)
        graph.find_bond_between(chain[0], chain[1]).order = 2
        graph.find_bond_between(chain[2], chain[3]).order = 2
        graph.find_bond_between(chain[4], chain[5]).order = 2
        self.assertEqual(iupac_name(graph), "hexa-1,3,5-triene")

    def test_hex_en_yne(self):
        graph = MolGraph()
        chain = build_linear_chain(graph, 6)
        graph.find_bond_between(chain[0], chain[1]).order = 2
        graph.find_bond_between(chain[4], chain[5]).order = 3
        self.assertEqual(iupac_name(graph), "hex-1-en-5-yne")

    def test_cyclohexadiene(self):
        graph = MolGraph()
        ring = build_ring(graph, 6)
        graph.find_bond_between(ring[0], ring[1]).order = 2
        graph.find_bond_between(ring[2], ring[3]).order = 2
        self.assertEqual(iupac_name(graph), "cyclohex-1,3-diene")


if __name__ == "__main__":
    unittest.main()
