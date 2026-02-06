"""Pruebas unitarias para test_chemname_pr13."""

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


class ChemNamePR13Test(unittest.TestCase):
    """Casos de prueba para ChemNamePR13Test."""
    def test_pyridine(self):
        """Verifica pyridine.

        Returns:
            None.

        """
        graph = MolGraph()
        build_kekule_ring(graph, ["N", "C", "C", "C", "C", "C"])
        self.assertEqual(iupac_name(graph), "pyridine")

    def test_chloropyridine_2(self):
        """Verifica chloropyridine 2.

        Returns:
            None.

        """
        graph = MolGraph()
        ring = build_kekule_ring(graph, ["N", "C", "C", "C", "C", "C"])
        cl = graph.add_atom("Cl", -1.0, 0.0)
        graph.add_bond(ring[1], cl.id, order=1)
        self.assertEqual(iupac_name(graph), "2-chloropyridine")

    def test_methylpyridine_3(self):
        """Verifica methylpyridine 3.

        Returns:
            None.

        """
        graph = MolGraph()
        ring = build_kekule_ring(graph, ["N", "C", "C", "C", "C", "C"])
        me = graph.add_atom("C", -1.0, 0.0)
        graph.add_bond(ring[2], me.id, order=1)
        self.assertEqual(iupac_name(graph), "3-methylpyridine")

    def test_furan(self):
        """Verifica furan.

        Returns:
            None.

        """
        graph = MolGraph()
        build_kekule_ring(graph, ["O", "C", "C", "C", "C"])
        self.assertEqual(iupac_name(graph), "furan")

    def test_thiophene(self):
        """Verifica thiophene.

        Returns:
            None.

        """
        graph = MolGraph()
        build_kekule_ring(graph, ["S", "C", "C", "C", "C"])
        self.assertEqual(iupac_name(graph), "thiophene")

    def test_pyrrole(self):
        """Verifica pyrrole.

        Returns:
            None.

        """
        graph = MolGraph()
        build_kekule_ring(graph, ["N", "C", "C", "C", "C"], explicit_h={0: 1})
        self.assertEqual(iupac_name(graph), "pyrrole")


if __name__ == "__main__":
    unittest.main()
