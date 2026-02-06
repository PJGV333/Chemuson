"""Pruebas unitarias para test_chemname_pr1."""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import MolGraph
from chemname.molview import MolView


def build_linear_chain(graph: MolGraph, length: int) -> list[int]:
    """Función de prueba auxiliar para build linear chain.

    Args:
        graph: Descripción del parámetro.
        length: Descripción del parámetro.

    Returns:
        None.

    """
    ids = []
    prev_id = None
    for i in range(length):
        atom = graph.add_atom("C", float(i), 0.0)
        ids.append(atom.id)
        if prev_id is not None:
            graph.add_bond(prev_id, atom.id, order=1)
        prev_id = atom.id
    return ids


class ChemNamePR1Test(unittest.TestCase):
    """Casos de prueba para ChemNamePR1Test."""
    def test_is_acyclic_true(self):
        """Verifica is acyclic true.

        Returns:
            None.

        """
        graph = MolGraph()
        build_linear_chain(graph, 3)
        view = MolView(graph)
        self.assertTrue(view.is_acyclic())

    def test_is_acyclic_false(self):
        """Verifica is acyclic false.

        Returns:
            None.

        """
        graph = MolGraph()
        a1 = graph.add_atom("C", 0.0, 0.0)
        a2 = graph.add_atom("C", 1.0, 0.0)
        a3 = graph.add_atom("C", 0.5, 1.0)
        graph.add_bond(a1.id, a2.id, order=1)
        graph.add_bond(a2.id, a3.id, order=1)
        graph.add_bond(a3.id, a1.id, order=1)

        view = MolView(graph)
        self.assertFalse(view.is_acyclic())


if __name__ == "__main__":
    unittest.main()
