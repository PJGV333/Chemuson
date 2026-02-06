"""Pruebas unitarias para test_chemname_pr19."""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import MolGraph
from chemname import iupac_name


def build_naphthalene(graph: MolGraph) -> list[int]:
    """Función de prueba auxiliar para build naphthalene.

    Args:
        graph: Descripción del parámetro.

    Returns:
        None.

    """
    atoms = [graph.add_atom("C", float(i), 0.0) for i in range(10)]
    ring1 = [0, 1, 2, 3, 4, 5]
    ring2 = [3, 4, 6, 7, 8, 9]

    for ring in (ring1, ring2):
        for i in range(len(ring)):
            a1 = atoms[ring[i]].id
            a2 = atoms[ring[(i + 1) % len(ring)]].id
            if graph.find_bond_between(a1, a2) is None:
                graph.add_bond(a1, a2, order=1, is_aromatic=True)
    return [atom.id for atom in atoms]


class ChemNamePR19Test(unittest.TestCase):
    """Casos de prueba para ChemNamePR19Test."""
    def test_1_2_dichloronaphthalene(self):
        """Verifica 1 2 dichloronaphthalene.

        Returns:
            None.

        """
        graph = MolGraph()
        atoms = build_naphthalene(graph)
        cl1 = graph.add_atom("Cl", -1.0, 0.0)
        cl2 = graph.add_atom("Cl", -2.0, 0.0)
        graph.add_bond(atoms[5], cl1.id, order=1)
        graph.add_bond(atoms[0], cl2.id, order=1)
        self.assertEqual(iupac_name(graph), "1,2-dichloronaphthalene")

    def test_1_5_dibromonaphthalene(self):
        """Verifica 1 5 dibromonaphthalene.

        Returns:
            None.

        """
        graph = MolGraph()
        atoms = build_naphthalene(graph)
        br1 = graph.add_atom("Br", -1.0, 0.0)
        br2 = graph.add_atom("Br", -2.0, 0.0)
        graph.add_bond(atoms[5], br1.id, order=1)
        graph.add_bond(atoms[9], br2.id, order=1)
        self.assertEqual(iupac_name(graph), "1,5-dibromonaphthalene")

    def test_2_6_dimethylnaphthalene(self):
        """Verifica 2 6 dimethylnaphthalene.

        Returns:
            None.

        """
        graph = MolGraph()
        atoms = build_naphthalene(graph)
        m1 = graph.add_atom("C", -1.0, 0.0)
        m2 = graph.add_atom("C", -2.0, 0.0)
        graph.add_bond(atoms[0], m1.id, order=1)
        graph.add_bond(atoms[8], m2.id, order=1)
        self.assertEqual(iupac_name(graph), "2,6-dimethylnaphthalene")


if __name__ == "__main__":
    unittest.main()
