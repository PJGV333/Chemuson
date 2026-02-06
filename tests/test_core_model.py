"""Pruebas unitarias para test_core_model."""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import BondStyle, BondStereo, MolGraph


class MolGraphTest(unittest.TestCase):
    """Casos de prueba para MolGraphTest."""
    def test_add_atom_and_bond(self):
        """Verifica add atom and bond.

        Returns:
            None.

        """
        graph = MolGraph()
        a1 = graph.add_atom("C", 0.0, 0.0)
        a2 = graph.add_atom("O", 1.0, 0.0)
        bond = graph.add_bond(a1.id, a2.id, order=2, style=BondStyle.PLAIN, stereo=BondStereo.NONE)

        self.assertEqual(len(graph.atoms), 2)
        self.assertEqual(len(graph.bonds), 1)
        self.assertEqual(bond.order, 2)

    def test_update_bond(self):
        """Verifica update bond.

        Returns:
            None.

        """
        graph = MolGraph()
        a1 = graph.add_atom("C", 0.0, 0.0)
        a2 = graph.add_atom("C", 1.0, 0.0)
        bond = graph.add_bond(a1.id, a2.id, order=1, style=BondStyle.PLAIN, stereo=BondStereo.NONE)

        graph.update_bond(bond.id, order=3, style=BondStyle.WAVY)
        updated = graph.get_bond(bond.id)
        self.assertEqual(updated.order, 3)
        self.assertEqual(updated.style, BondStyle.WAVY)

    def test_validate_allows_hypervalent_phosphorus(self):
        """Verifica validate allows hypervalent phosphorus.

        Returns:
            None.

        """
        graph = MolGraph()
        p = graph.add_atom("P", 0.0, 0.0)
        o_dbl = graph.add_atom("O", 0.0, 1.0)
        o1 = graph.add_atom("O", -1.0, 0.0)
        o2 = graph.add_atom("O", -1.0, -1.0)
        c = graph.add_atom("C", 1.0, 0.0)

        graph.add_bond(p.id, o_dbl.id, order=2)
        graph.add_bond(p.id, o1.id, order=1)
        graph.add_bond(p.id, o2.id, order=1)
        graph.add_bond(p.id, c.id, order=1)

        self.assertNotIn(p.id, graph.validate())

    def test_validate_allows_sf6(self):
        """Verifica validate allows sf6.

        Returns:
            None.

        """
        graph = MolGraph()
        s = graph.add_atom("S", 0.0, 0.0)
        fs = [graph.add_atom("F", float(i), 1.0) for i in range(6)]
        for f in fs:
            graph.add_bond(s.id, f.id, order=1)
        self.assertNotIn(s.id, graph.validate())

    def test_validate_allows_if7(self):
        """Verifica validate allows if7.

        Returns:
            None.

        """
        graph = MolGraph()
        i = graph.add_atom("I", 0.0, 0.0)
        fs = [graph.add_atom("F", float(n), 1.0) for n in range(7)]
        for f in fs:
            graph.add_bond(i.id, f.id, order=1)
        self.assertNotIn(i.id, graph.validate())

    def test_validate_allows_ammonium(self):
        """Verifica validate allows ammonium.

        Returns:
            None.

        """
        graph = MolGraph()
        n = graph.add_atom("N", 0.0, 0.0)
        hs = [graph.add_atom("H", float(k), 1.0) for k in range(4)]
        for h in hs:
            graph.add_bond(n.id, h.id, order=1)
        self.assertNotIn(n.id, graph.validate())

    def test_validate_still_flags_overvalent_carbon(self):
        """Verifica validate still flags overvalent carbon.

        Returns:
            None.

        """
        graph = MolGraph()
        c = graph.add_atom("C", 0.0, 0.0)
        hs = [graph.add_atom("H", float(k), 1.0) for k in range(5)]
        for h in hs:
            graph.add_bond(c.id, h.id, order=1)
        self.assertIn(c.id, graph.validate())


if __name__ == "__main__":
    unittest.main()
