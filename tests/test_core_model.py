import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import BondStyle, BondStereo, MolGraph


class MolGraphTest(unittest.TestCase):
    def test_add_atom_and_bond(self):
        graph = MolGraph()
        a1 = graph.add_atom("C", 0.0, 0.0)
        a2 = graph.add_atom("O", 1.0, 0.0)
        bond = graph.add_bond(a1.id, a2.id, order=2, style=BondStyle.PLAIN, stereo=BondStereo.NONE)

        self.assertEqual(len(graph.atoms), 2)
        self.assertEqual(len(graph.bonds), 1)
        self.assertEqual(bond.order, 2)

    def test_update_bond(self):
        graph = MolGraph()
        a1 = graph.add_atom("C", 0.0, 0.0)
        a2 = graph.add_atom("C", 1.0, 0.0)
        bond = graph.add_bond(a1.id, a2.id, order=1, style=BondStyle.PLAIN, stereo=BondStereo.NONE)

        graph.update_bond(bond.id, order=3, style=BondStyle.WAVY)
        updated = graph.get_bond(bond.id)
        self.assertEqual(updated.order, 3)
        self.assertEqual(updated.style, BondStyle.WAVY)


if __name__ == "__main__":
    unittest.main()
