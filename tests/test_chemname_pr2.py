import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import MolGraph
from chemname.molview import MolView
from chemname.parent_chain import longest_carbon_chain


def build_linear_chain(graph: MolGraph, length: int) -> list[int]:
    ids = []
    prev_id = None
    for i in range(length):
        atom = graph.add_atom("C", float(i), 0.0)
        ids.append(atom.id)
        if prev_id is not None:
            graph.add_bond(prev_id, atom.id, order=1)
        prev_id = atom.id
    return ids


def is_valid_chain(view: MolView, chain: list[int]) -> bool:
    for i in range(len(chain) - 1):
        if chain[i + 1] not in view.neighbors(chain[i]):
            return False
    return True


class ChemNamePR2Test(unittest.TestCase):
    def test_longest_chain_dodecane(self):
        graph = MolGraph()
        build_linear_chain(graph, 12)
        view = MolView(graph)
        chain = longest_carbon_chain(view)
        self.assertEqual(len(chain), 12)
        self.assertTrue(is_valid_chain(view, chain))

    def test_longest_chain_branch(self):
        graph = MolGraph()
        main = build_linear_chain(graph, 4)
        b1 = graph.add_atom("C", 1.0, 1.0)
        graph.add_bond(main[1], b1.id, order=1)
        b2 = graph.add_atom("C", 1.0, 2.0)
        graph.add_bond(b1.id, b2.id, order=1)

        view = MolView(graph)
        chain = longest_carbon_chain(view)
        self.assertEqual(len(chain), 5)
        self.assertTrue(is_valid_chain(view, chain))


if __name__ == "__main__":
    unittest.main()
