"""Pruebas unitarias para test_chemname_pr9."""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import MolGraph
from chemname import iupac_name


def build_ring(graph: MolGraph, size: int) -> list[int]:
    """Función de prueba auxiliar para build ring.

    Args:
        graph: Descripción del parámetro.
        size: Descripción del parámetro.

    Returns:
        None.

    """
    atoms = [graph.add_atom("C", float(i), 0.0) for i in range(size)]
    for i in range(size):
        graph.add_bond(atoms[i].id, atoms[(i + 1) % size].id, order=1)
    return [atom.id for atom in atoms]


class ChemNamePR9Test(unittest.TestCase):
    """Casos de prueba para ChemNamePR9Test."""
    def test_cyclohexane(self):
        """Verifica cyclohexane.

        Returns:
            None.

        """
        graph = MolGraph()
        build_ring(graph, 6)
        self.assertEqual(iupac_name(graph), "cyclohexane")

    def test_methylcyclohexane(self):
        """Verifica methylcyclohexane.

        Returns:
            None.

        """
        graph = MolGraph()
        ring = build_ring(graph, 6)
        methyl = graph.add_atom("C", -1.0, 0.0)
        graph.add_bond(ring[0], methyl.id, order=1)
        self.assertEqual(iupac_name(graph), "methylcyclohexane")

    def test_chlorocyclohexane(self):
        """Verifica chlorocyclohexane.

        Returns:
            None.

        """
        graph = MolGraph()
        ring = build_ring(graph, 6)
        cl = graph.add_atom("Cl", -1.0, 0.0)
        graph.add_bond(ring[0], cl.id, order=1)
        self.assertEqual(iupac_name(graph), "chlorocyclohexane")


if __name__ == "__main__":
    unittest.main()
