"""Pruebas unitarias para test_chemname_pr25."""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import MolGraph
from chemname import iupac_name


def build_ring(graph: MolGraph, elements: list[str], double_bonds: set[int]) -> list[int]:
    """Función de prueba auxiliar para build ring.

    Args:
        graph: Descripción del parámetro.
        elements: Descripción del parámetro.
        double_bonds: Descripción del parámetro.

    Returns:
        None.

    """
    atoms = [graph.add_atom(elem, float(i), 0.0) for i, elem in enumerate(elements)]
    size = len(atoms)
    for i in range(size):
        order = 2 if i in double_bonds else 1
        graph.add_bond(atoms[i].id, atoms[(i + 1) % size].id, order=order)
    return [atom.id for atom in atoms]


class ChemNamePR25Test(unittest.TestCase):
    """Casos de prueba para ChemNamePR25Test."""
    def test_triazine(self):
        """Verifica triazine.

        Returns:
            None.

        """
        graph = MolGraph()
        elements = ["N", "C", "N", "C", "N", "C"]
        ring = build_ring(graph, elements, {0, 2, 4})
        self.assertEqual(iupac_name(graph), "triazine")

    def test_chlorotriazine(self):
        """Verifica chlorotriazine.

        Returns:
            None.

        """
        graph = MolGraph()
        elements = ["N", "C", "N", "C", "N", "C"]
        ring = build_ring(graph, elements, {0, 2, 4})
        cl = graph.add_atom("Cl", -1.0, 0.0)
        graph.add_bond(ring[1], cl.id, order=1)
        self.assertEqual(iupac_name(graph), "2-chlorotriazine")

    def test_triazole(self):
        """Verifica triazole.

        Returns:
            None.

        """
        graph = MolGraph()
        elements = ["N", "N", "N", "C", "C"]
        atoms = [graph.add_atom(elem, float(i), 0.0, explicit_h=1 if i == 0 else None) for i, elem in enumerate(elements)]
        for i in range(len(atoms)):
            order = 2 if i in {0, 2} else 1
            graph.add_bond(atoms[i].id, atoms[(i + 1) % len(atoms)].id, order=order)
        self.assertEqual(iupac_name(graph), "triazole")


if __name__ == "__main__":
    unittest.main()
