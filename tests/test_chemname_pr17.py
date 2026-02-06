"""Pruebas unitarias para test_chemname_pr17."""

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


def build_linear_chain(graph: MolGraph, length: int) -> list[int]:
    """Función de prueba auxiliar para build linear chain.

    Args:
        graph: Descripción del parámetro.
        length: Descripción del parámetro.

    Returns:
        None.

    """
    ids = []
    prev = None
    for i in range(length):
        atom = graph.add_atom("C", float(i), 0.0)
        ids.append(atom.id)
        if prev is not None:
            graph.add_bond(prev, atom.id, order=1)
        prev = atom.id
    return ids


class ChemNamePR17Test(unittest.TestCase):
    """Casos de prueba para ChemNamePR17Test."""
    def test_phenylheptane(self):
        """Verifica phenylheptane.

        Returns:
            None.

        """
        graph = MolGraph()
        chain = build_linear_chain(graph, 7)
        ring = build_benzene_kekule(graph)
        graph.add_bond(chain[0], ring[0], order=1)
        self.assertEqual(iupac_name(graph), "1-phenylheptane")

    def test_benzylheptane(self):
        """Verifica benzylheptane.

        Returns:
            None.

        """
        graph = MolGraph()
        chain = build_linear_chain(graph, 7)
        benzyl = graph.add_atom("C", -1.0, 0.0)
        ring = build_benzene_kekule(graph)
        graph.add_bond(chain[0], benzyl.id, order=1)
        graph.add_bond(benzyl.id, ring[0], order=1)
        self.assertEqual(iupac_name(graph), "1-benzylheptane")

    def test_cyclohexylheptane(self):
        """Verifica cyclohexylheptane.

        Returns:
            None.

        """
        graph = MolGraph()
        chain = build_linear_chain(graph, 7)
        ring = build_cyclohexane(graph)
        graph.add_bond(chain[0], ring[0], order=1)
        self.assertEqual(iupac_name(graph), "1-cyclohexylheptane")


if __name__ == "__main__":
    unittest.main()
