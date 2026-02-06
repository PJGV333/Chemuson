"""Pruebas unitarias para test_chemname_pr11."""

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


class ChemNamePR11Test(unittest.TestCase):
    """Casos de prueba para ChemNamePR11Test."""
    def test_benzene(self):
        """Verifica benzene.

        Returns:
            None.

        """
        graph = MolGraph()
        build_benzene_kekule(graph)
        self.assertEqual(iupac_name(graph), "benzene")

    def test_chlorobenzene(self):
        """Verifica chlorobenzene.

        Returns:
            None.

        """
        graph = MolGraph()
        ring = build_benzene_kekule(graph)
        cl = graph.add_atom("Cl", -1.0, 0.0)
        graph.add_bond(ring[0], cl.id, order=1)
        self.assertEqual(iupac_name(graph), "chlorobenzene")

    def test_methylbenzene(self):
        """Verifica methylbenzene.

        Returns:
            None.

        """
        graph = MolGraph()
        ring = build_benzene_kekule(graph)
        me = graph.add_atom("C", -1.0, 0.0)
        graph.add_bond(ring[0], me.id, order=1)
        self.assertEqual(iupac_name(graph), "methylbenzene")

    def test_bromochlorobenzene(self):
        """Verifica bromochlorobenzene.

        Returns:
            None.

        """
        graph = MolGraph()
        ring = build_benzene_kekule(graph)
        br = graph.add_atom("Br", -1.0, 0.0)
        cl = graph.add_atom("Cl", 3.0, 1.0)
        graph.add_bond(ring[0], br.id, order=1)
        graph.add_bond(ring[3], cl.id, order=1)
        self.assertEqual(iupac_name(graph), "1-bromo-4-chlorobenzene")

    def test_hydroxybenzene(self):
        """Verifica hydroxybenzene.

        Returns:
            None.

        """
        graph = MolGraph()
        ring = build_benzene_kekule(graph)
        o = graph.add_atom("O", -1.0, 0.0)
        graph.add_bond(ring[0], o.id, order=1)
        self.assertEqual(iupac_name(graph), "hydroxybenzene")

    def test_nitrobenzene(self):
        """Verifica nitrobenzene.

        Returns:
            None.

        """
        graph = MolGraph()
        ring = build_benzene_kekule(graph)
        n = graph.add_atom("N", -1.0, 0.0)
        o1 = graph.add_atom("O", -2.0, 0.0)
        o2 = graph.add_atom("O", -2.0, 1.0)
        graph.add_bond(ring[0], n.id, order=1)
        graph.add_bond(n.id, o1.id, order=2)
        graph.add_bond(n.id, o2.id, order=1)
        self.assertEqual(iupac_name(graph), "nitrobenzene")


if __name__ == "__main__":
    unittest.main()
