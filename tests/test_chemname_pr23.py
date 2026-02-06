"""Pruebas unitarias para test_chemname_pr23."""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import MolGraph
from chemname import iupac_name


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


class ChemNamePR23Test(unittest.TestCase):
    """Casos de prueba para ChemNamePR23Test."""
    def test_propanal(self):
        """Verifica propanal.

        Returns:
            None.

        """
        graph = MolGraph()
        chain = build_linear_chain(graph, 3)
        o = graph.add_atom("O", 3.0, 0.0)
        graph.add_bond(chain[-1], o.id, order=2)
        self.assertEqual(iupac_name(graph), "propanal")

    def test_propan_2_one(self):
        """Verifica propan 2 one.

        Returns:
            None.

        """
        graph = MolGraph()
        chain = build_linear_chain(graph, 3)
        o = graph.add_atom("O", 1.0, 1.0)
        graph.add_bond(chain[1], o.id, order=2)
        self.assertEqual(iupac_name(graph), "propan-2-one")

    def test_butanoic_acid(self):
        """Verifica butanoic acid.

        Returns:
            None.

        """
        graph = MolGraph()
        chain = build_linear_chain(graph, 4)
        o1 = graph.add_atom("O", 4.0, 0.0)
        graph.add_bond(chain[-1], o1.id, order=2)
        o2 = graph.add_atom("O", 4.0, 1.0, explicit_h=1)
        graph.add_bond(chain[-1], o2.id, order=1)
        self.assertEqual(iupac_name(graph), "butanoic acid")

    def test_propanenitrile(self):
        """Verifica propanenitrile.

        Returns:
            None.

        """
        graph = MolGraph()
        chain = build_linear_chain(graph, 3)
        n = graph.add_atom("N", 3.0, 0.0)
        graph.add_bond(chain[-1], n.id, order=3)
        self.assertEqual(iupac_name(graph), "propanenitrile")


if __name__ == "__main__":
    unittest.main()
