"""Pruebas unitarias para test_chemname_pr3."""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import MolGraph
from chemname.locants import choose_orientation, substituents_on_chain
from chemname.molview import MolView
from chemname.options import NameOptions
from chemname.parent_chain import longest_carbon_chain


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


class ChemNamePR3Test(unittest.TestCase):
    """Casos de prueba para ChemNamePR3Test."""
    def test_orientation_prefers_bromo_at_one(self):
        """Verifica orientation prefers bromo at one.

        Returns:
            None.

        """
        graph = MolGraph()
        chain_ids = build_linear_chain(graph, 12)

        cl = graph.add_atom("Cl", -1.0, 0.0)
        graph.add_bond(chain_ids[0], cl.id, order=1)
        br = graph.add_atom("Br", 13.0, 0.0)
        graph.add_bond(chain_ids[-1], br.id, order=1)

        view = MolView(graph)
        chain = longest_carbon_chain(view)
        oriented = choose_orientation(view, chain, NameOptions())
        subs = substituents_on_chain(view, oriented)

        locants = {sub.name: sub.locant for sub in subs}
        self.assertEqual(locants.get("bromo"), 1)
        self.assertEqual(locants.get("chloro"), 12)


if __name__ == "__main__":
    unittest.main()
