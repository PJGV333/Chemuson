"""Pruebas unitarias para test_chemname_pr16."""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import MolGraph
from chemname import iupac_name


def build_benzene_kekule(graph: MolGraph) -> list[int]:
    """Función de prueba auxiliar para build benzene kekule.

    Args:
        graph: Descripción del parámetro.

    Returns:
        None.

    """
    atoms = [graph.add_atom("C", float(i), 0.0) for i in range(6)]
    for i in range(6):
        order = 2 if i % 2 == 0 else 1
        graph.add_bond(atoms[i].id, atoms[(i + 1) % 6].id, order=order)
    return [atom.id for atom in atoms]


def build_cyclohexane(graph: MolGraph) -> list[int]:
    """Función de prueba auxiliar para build cyclohexane.

    Args:
        graph: Descripción del parámetro.

    Returns:
        None.

    """
    atoms = [graph.add_atom("C", float(i), 0.0) for i in range(6)]
    for i in range(6):
        graph.add_bond(atoms[i].id, atoms[(i + 1) % 6].id, order=1)
    return [atom.id for atom in atoms]


class ChemNamePR16Test(unittest.TestCase):
    """Casos de prueba para ChemNamePR16Test."""
    def test_isopropylbenzene(self):
        """Verifica isopropylbenzene.

        Returns:
            None.

        """
        graph = MolGraph()
        ring = build_benzene_kekule(graph)
        root = graph.add_atom("C", -1.0, 0.0)
        m1 = graph.add_atom("C", -2.0, 1.0)
        m2 = graph.add_atom("C", -2.0, -1.0)
        graph.add_bond(ring[0], root.id, order=1)
        graph.add_bond(root.id, m1.id, order=1)
        graph.add_bond(root.id, m2.id, order=1)
        self.assertEqual(iupac_name(graph), "isopropylbenzene")

    def test_tert_butylbenzene(self):
        """Verifica tert butylbenzene.

        Returns:
            None.

        """
        graph = MolGraph()
        ring = build_benzene_kekule(graph)
        root = graph.add_atom("C", -1.0, 0.0)
        m1 = graph.add_atom("C", -2.0, 1.0)
        m2 = graph.add_atom("C", -2.0, -1.0)
        m3 = graph.add_atom("C", -2.0, 0.0)
        graph.add_bond(ring[0], root.id, order=1)
        graph.add_bond(root.id, m1.id, order=1)
        graph.add_bond(root.id, m2.id, order=1)
        graph.add_bond(root.id, m3.id, order=1)
        self.assertEqual(iupac_name(graph), "tert-butylbenzene")

    def test_sec_butylcyclohexane(self):
        """Verifica sec butylcyclohexane.

        Returns:
            None.

        """
        graph = MolGraph()
        ring = build_cyclohexane(graph)
        root = graph.add_atom("C", -1.0, 0.0)
        m1 = graph.add_atom("C", -2.0, 1.0)
        e1 = graph.add_atom("C", -2.0, -1.0)
        e2 = graph.add_atom("C", -3.0, -1.0)
        graph.add_bond(ring[0], root.id, order=1)
        graph.add_bond(root.id, m1.id, order=1)
        graph.add_bond(root.id, e1.id, order=1)
        graph.add_bond(e1.id, e2.id, order=1)
        self.assertEqual(iupac_name(graph), "sec-butylcyclohexane")


if __name__ == "__main__":
    unittest.main()
