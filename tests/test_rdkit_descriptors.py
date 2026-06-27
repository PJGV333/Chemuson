"""Descriptores químicos calculados por worker RDKit aislado."""

from __future__ import annotations

import unittest


from chemuson.chemio.rdkit_safe import molecular_descriptors_isolated
from chemuson.core.model import MolGraph


class RdkitDescriptorWorkerTest(unittest.TestCase):
    def test_descriptor_worker_returns_lipinski_fields(self):
        graph = MolGraph()
        c1 = graph.add_atom("C", 0.0, 0.0)
        c2 = graph.add_atom("C", 40.0, 0.0)
        o = graph.add_atom("O", 80.0, 0.0)
        graph.add_bond(c1.id, c2.id, order=1)
        graph.add_bond(c2.id, o.id, order=1)

        descriptors, error = molecular_descriptors_isolated(graph, timeout_s=5.0)
        if error == "rdkit_unavailable":
            self.skipTest("RDKit no disponible en el worker")

        self.assertIsNone(error)
        self.assertIsNotNone(descriptors)
        assert descriptors is not None
        self.assertIn("logp", descriptors)
        self.assertIn("tpsa", descriptors)
        self.assertIn("hbd", descriptors)
        self.assertIn("hba", descriptors)
        self.assertIn("rotatable_bonds", descriptors)
        self.assertIn("lipinski_violations", descriptors)
        self.assertGreaterEqual(float(descriptors["tpsa"]), 0.0)
        self.assertGreaterEqual(int(descriptors["hbd"]), 1)


if __name__ == "__main__":
    unittest.main()
