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

        self.assertIsNone(error, msg=f"RDKit isolated worker failed: {error}")
        self.assertIsNotNone(descriptors)
        assert descriptors is not None
        self.assertIn("logp", descriptors)
        self.assertIn("tpsa", descriptors)
        self.assertIn("hbd", descriptors)
        self.assertIn("hba", descriptors)
        self.assertIn("rotatable_bonds", descriptors)
        self.assertIn("lipinski_violations", descriptors)
        self.assertAlmostEqual(float(descriptors["logp"]), -0.0014, delta=1e-4)
        self.assertAlmostEqual(float(descriptors["tpsa"]), 20.23, delta=1e-6)
        self.assertEqual(int(descriptors["hbd"]), 1)
        self.assertEqual(int(descriptors["hba"]), 1)
        self.assertAlmostEqual(float(descriptors["molecular_weight"]), 46.069, delta=1e-6)


if __name__ == "__main__":
    unittest.main()
