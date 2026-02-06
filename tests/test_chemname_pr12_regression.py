"""Pruebas unitarias para test_chemname_pr12_regression."""

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


class ChemNamePR12Regression(unittest.TestCase):
    """Casos de prueba para ChemNamePR12Regression."""
    def test_linear_regressions(self):
        """Verifica linear regressions.

        Returns:
            None.

        """
        graph = MolGraph()
        chain = build_linear_chain(graph, 12)
        cl = graph.add_atom("Cl", -1.0, 0.0)
        br = graph.add_atom("Br", 13.0, 0.0)
        graph.add_bond(chain[0], cl.id, order=1)
        graph.add_bond(chain[-1], br.id, order=1)
        self.assertEqual(iupac_name(graph), "1-bromo-12-chlorododecane")

        graph = MolGraph()
        c1 = graph.add_atom("C", 0.0, 0.0)
        cl = graph.add_atom("Cl", 1.0, 0.0)
        graph.add_bond(c1.id, cl.id, order=1)
        self.assertEqual(iupac_name(graph), "1-chloromethane")

        graph = MolGraph()
        chain = build_linear_chain(graph, 3)
        cl1 = graph.add_atom("Cl", -1.0, 0.0)
        cl2 = graph.add_atom("Cl", 3.0, 0.0)
        graph.add_bond(chain[0], cl1.id, order=1)
        graph.add_bond(chain[-1], cl2.id, order=1)
        self.assertEqual(iupac_name(graph), "1,3-dichloropropane")

        graph = MolGraph()
        c1 = graph.add_atom("C", 0.0, 0.0)
        c2 = graph.add_atom("C", 1.0, 0.0)
        c3 = graph.add_atom("C", 2.0, 0.0)
        graph.add_bond(c1.id, c2.id, order=2)
        graph.add_bond(c2.id, c3.id, order=1)
        self.assertEqual(iupac_name(graph), "prop-1-ene")

        graph = MolGraph()
        chain = build_linear_chain(graph, 3)
        o = graph.add_atom("O", 3.0, 0.0)
        graph.add_bond(chain[-1], o.id, order=1)
        self.assertEqual(iupac_name(graph), "propan-1-ol")

    def test_ring_regressions(self):
        """Verifica ring regressions.

        Returns:
            None.

        """
        graph = MolGraph()
        build_ring(graph, 6)
        self.assertEqual(iupac_name(graph), "cyclohexane")

        graph = MolGraph()
        ring = build_ring(graph, 6)
        me = graph.add_atom("C", -1.0, 0.0)
        graph.add_bond(ring[0], me.id, order=1)
        self.assertEqual(iupac_name(graph), "methylcyclohexane")

        graph = MolGraph()
        ring = build_ring(graph, 6)
        cl = graph.add_atom("Cl", -1.0, 0.0)
        graph.add_bond(ring[0], cl.id, order=1)
        self.assertEqual(iupac_name(graph), "chlorocyclohexane")

        graph = MolGraph()
        ring = build_ring(graph, 6)
        br = graph.add_atom("Br", -1.0, 0.0)
        cl = graph.add_atom("Cl", 3.0, 1.0)
        graph.add_bond(ring[0], br.id, order=1)
        graph.add_bond(ring[3], cl.id, order=1)
        self.assertEqual(iupac_name(graph), "1-bromo-4-chlorocyclohexane")

        graph = MolGraph()
        ring = build_benzene_kekule(graph)
        self.assertEqual(iupac_name(graph), "benzene")

        graph = MolGraph()
        ring = build_benzene_kekule(graph)
        cl = graph.add_atom("Cl", -1.0, 0.0)
        graph.add_bond(ring[0], cl.id, order=1)
        self.assertEqual(iupac_name(graph), "chlorobenzene")

        graph = MolGraph()
        ring = build_benzene_kekule(graph)
        me = graph.add_atom("C", -1.0, 0.0)
        graph.add_bond(ring[0], me.id, order=1)
        self.assertEqual(iupac_name(graph), "methylbenzene")


if __name__ == "__main__":
    unittest.main()
