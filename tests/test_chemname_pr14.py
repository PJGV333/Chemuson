"""Pruebas unitarias para test_chemname_pr14."""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import MolGraph
from chemname import iupac_name


def build_kekule_ring(graph: MolGraph, elements, explicit_h=None):
    """Función de prueba auxiliar para build kekule ring.

    Args:
        graph: Descripción del parámetro.
        elements: Descripción del parámetro.
        explicit_h: Descripción del parámetro.

    Returns:
        None.

    """
    explicit_h = explicit_h or {}
    atoms = []
    for i, elem in enumerate(elements):
        atoms.append(graph.add_atom(elem, float(i), 0.0, explicit_h=explicit_h.get(i)))
    size = len(elements)
    for i in range(size):
        if size == 5:
            order = 2 if i in {0, 2} else 1
        else:
            order = 2 if i % 2 == 0 else 1
        graph.add_bond(atoms[i].id, atoms[(i + 1) % size].id, order=order)
    return [atom.id for atom in atoms]


class ChemNamePR14Test(unittest.TestCase):
    """Casos de prueba para ChemNamePR14Test."""
    def test_pyridazine(self):
        """Verifica pyridazine.

        Returns:
            None.

        """
        graph = MolGraph()
        build_kekule_ring(graph, ["N", "N", "C", "C", "C", "C"])
        self.assertEqual(iupac_name(graph), "pyridazine")

    def test_pyrimidine(self):
        """Verifica pyrimidine.

        Returns:
            None.

        """
        graph = MolGraph()
        ring = build_kekule_ring(graph, ["N", "C", "N", "C", "C", "C"])
        self.assertEqual(iupac_name(graph), "pyrimidine")
        cl = graph.add_atom("Cl", -1.0, 0.0)
        graph.add_bond(ring[1], cl.id, order=1)
        self.assertEqual(iupac_name(graph), "2-chloropyrimidine")

    def test_pyrazine(self):
        """Verifica pyrazine.

        Returns:
            None.

        """
        graph = MolGraph()
        build_kekule_ring(graph, ["N", "C", "C", "N", "C", "C"])
        self.assertEqual(iupac_name(graph), "pyrazine")

    def test_imidazole(self):
        """Verifica imidazole.

        Returns:
            None.

        """
        graph = MolGraph()
        build_kekule_ring(graph, ["N", "C", "N", "C", "C"], explicit_h={0: 1})
        self.assertEqual(iupac_name(graph), "imidazole")

    def test_oxazole(self):
        """Verifica oxazole.

        Returns:
            None.

        """
        graph = MolGraph()
        build_kekule_ring(graph, ["O", "C", "N", "C", "C"])
        self.assertEqual(iupac_name(graph), "oxazole")

    def test_thiazole(self):
        """Verifica thiazole.

        Returns:
            None.

        """
        graph = MolGraph()
        build_kekule_ring(graph, ["S", "C", "N", "C", "C"])
        self.assertEqual(iupac_name(graph), "thiazole")


if __name__ == "__main__":
    unittest.main()
