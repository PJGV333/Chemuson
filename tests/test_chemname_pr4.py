"""Pruebas unitarias para test_chemname_pr4."""

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


class ChemNamePR4Test(unittest.TestCase):
    """Casos de prueba para ChemNamePR4Test."""
    def test_bromo_chloro_dodecane(self):
        """Verifica bromo chloro dodecane.

        Returns:
            None.

        """
        graph = MolGraph()
        chain_ids = build_linear_chain(graph, 12)

        cl = graph.add_atom("Cl", -1.0, 0.0)
        graph.add_bond(chain_ids[0], cl.id, order=1)
        br = graph.add_atom("Br", 13.0, 0.0)
        graph.add_bond(chain_ids[-1], br.id, order=1)

        self.assertEqual(iupac_name(graph), "1-bromo-12-chlorododecane")

    def test_chloromethane(self):
        """Verifica chloromethane.

        Returns:
            None.

        """
        graph = MolGraph()
        c1 = graph.add_atom("C", 0.0, 0.0)
        cl = graph.add_atom("Cl", 1.0, 0.0)
        graph.add_bond(c1.id, cl.id, order=1)

        self.assertEqual(iupac_name(graph), "1-chloromethane")

    def test_dichloropropane(self):
        """Verifica dichloropropane.

        Returns:
            None.

        """
        graph = MolGraph()
        chain_ids = build_linear_chain(graph, 3)

        cl1 = graph.add_atom("Cl", -1.0, 0.0)
        graph.add_bond(chain_ids[0], cl1.id, order=1)
        cl2 = graph.add_atom("Cl", 3.0, 0.0)
        graph.add_bond(chain_ids[-1], cl2.id, order=1)

        self.assertEqual(iupac_name(graph), "1,3-dichloropropane")


if __name__ == "__main__":
    unittest.main()
