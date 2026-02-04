import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import MolGraph
from chemname.molview import MolView
from chemname.rings import find_rings_simple, is_simple_ring, perceive_aromaticity_basic


def build_ring(graph: MolGraph, size: int, aromatic: bool = False, kekule: bool = False):
    atoms = [graph.add_atom("C", float(i), 0.0) for i in range(size)]
    for i in range(size):
        order = 1
        if kekule and i % 2 == 0:
            order = 2
        graph.add_bond(
            atoms[i].id,
            atoms[(i + 1) % size].id,
            order=order,
            is_aromatic=aromatic,
        )
    return atoms


class ChemNamePR8Test(unittest.TestCase):
    def test_cyclohexane_ring_non_aromatic(self):
        graph = MolGraph()
        build_ring(graph, 6)
        view = MolView(graph)
        rings = find_rings_simple(view)
        self.assertEqual(len(rings), 1)
        self.assertEqual(len(next(iter(rings))), 6)
        self.assertTrue(is_simple_ring(view, next(iter(rings))))
        aromatic = perceive_aromaticity_basic(view, rings)
        self.assertFalse(aromatic)

    def test_benzene_aromatic_flagged(self):
        graph = MolGraph()
        build_ring(graph, 6, aromatic=True)
        view = MolView(graph)
        rings = find_rings_simple(view)
        aromatic = perceive_aromaticity_basic(view, rings)
        self.assertEqual(len(aromatic), 1)

    def test_benzene_kekule_aromatic(self):
        graph = MolGraph()
        build_ring(graph, 6, kekule=True)
        view = MolView(graph)
        rings = find_rings_simple(view)
        aromatic = perceive_aromaticity_basic(view, rings)
        self.assertEqual(len(aromatic), 1)

    def test_fused_ring_no_crash(self):
        graph = MolGraph()
        atoms = [graph.add_atom("C", float(i), 0.0) for i in range(10)]
        # ring 1: 1-2-3-4-5-6-1
        ring1 = [0, 1, 2, 3, 4, 5]
        for i in range(len(ring1)):
            a1 = atoms[ring1[i]].id
            a2 = atoms[ring1[(i + 1) % len(ring1)]].id
            graph.add_bond(a1, a2, order=1)
        # ring 2: 4-5-7-8-9-10-4 (fused at 4-5)
        ring2 = [3, 4, 6, 7, 8, 9]
        for i in range(len(ring2)):
            a1 = atoms[ring2[i]].id
            a2 = atoms[ring2[(i + 1) % len(ring2)]].id
            if graph.find_bond_between(a1, a2) is None:
                graph.add_bond(a1, a2, order=1)

        view = MolView(graph)
        rings = find_rings_simple(view)
        self.assertTrue(len(rings) >= 1)
        aromatic = perceive_aromaticity_basic(view, rings)
        self.assertFalse(aromatic)


if __name__ == "__main__":
    unittest.main()
