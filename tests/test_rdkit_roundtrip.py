"""Pruebas unitarias para test_rdkit_roundtrip."""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.core.model import MolGraph
from chemuson.chemio.rdkit_io import (
    molfile_to_molgraph,
    molgraph_to_molfile,
    molgraph_to_smiles,
    molgraph_to_rdkit_with_map,
)

try:
    from rdkit import Chem  # noqa: F401
    RDKit_AVAILABLE = True
except Exception:
    RDKit_AVAILABLE = False


class RdkitRoundtripTest(unittest.TestCase):
    """Casos de prueba para RdkitRoundtripTest."""
    @unittest.skipIf(not RDKit_AVAILABLE, "RDKit no disponible")
    def test_molgraph_roundtrip_smiles(self):
        """Verifica molgraph roundtrip smiles.

        Returns:
            None.

        """
        graph = MolGraph()
        a1 = graph.add_atom("C", 0.0, 0.0)
        a2 = graph.add_atom("C", 1.5, 0.0)
        graph.add_bond(a1.id, a2.id, order=1)

        molfile = molgraph_to_molfile(graph)
        graph2 = molfile_to_molgraph(molfile)
        smiles = molgraph_to_smiles(graph2)

        self.assertIn(smiles, {"CC", "C-C"})

    @unittest.skipIf(not RDKit_AVAILABLE, "RDKit no disponible")
    def test_duplicate_bonds_do_not_crash_rdkit_conversion(self):
        graph = MolGraph()
        a1 = graph.add_atom("C", 0.0, 0.0)
        a2 = graph.add_atom("C", 1.5, 0.0)
        graph.add_bond(a1.id, a2.id, order=1)
        graph.add_bond(a1.id, a2.id, order=2)

        mol, id_map = molgraph_to_rdkit_with_map(graph)
        self.assertIsNotNone(mol)
        self.assertEqual(len(id_map), 2)
        self.assertEqual(mol.GetNumBonds(), 1)


if __name__ == "__main__":
    unittest.main()
