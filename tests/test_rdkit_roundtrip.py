import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import MolGraph
from chemio.rdkit_io import molfile_to_molgraph, molgraph_to_molfile, molgraph_to_smiles

try:
    from rdkit import Chem  # noqa: F401
    RDKit_AVAILABLE = True
except Exception:
    RDKit_AVAILABLE = False


class RdkitRoundtripTest(unittest.TestCase):
    @unittest.skipIf(not RDKit_AVAILABLE, "RDKit no disponible")
    def test_molgraph_roundtrip_smiles(self):
        graph = MolGraph()
        a1 = graph.add_atom("C", 0.0, 0.0)
        a2 = graph.add_atom("C", 1.5, 0.0)
        graph.add_bond(a1.id, a2.id, order=1)

        molfile = molgraph_to_molfile(graph)
        graph2 = molfile_to_molgraph(molfile)
        smiles = molgraph_to_smiles(graph2)

        self.assertIn(smiles, {"CC", "C-C"})


if __name__ == "__main__":
    unittest.main()
