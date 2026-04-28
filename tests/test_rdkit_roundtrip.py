"""Pruebas unitarias para test_rdkit_roundtrip."""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.core.model import MolGraph
from chemuson.chemio import rdkit_io, rdkit_safe
from chemuson.chemio.rdkit_io import (
    kekulize_display_orders,
    molfile_to_molgraph,
    molgraph_to_molfile,
    molgraph_to_smiles,
    molgraph_to_rdkit_with_map,
    smiles_to_molgraph,
)

try:
    from rdkit import Chem  # noqa: F401
    RDKit_AVAILABLE = True
except Exception:
    RDKit_AVAILABLE = False


class RdkitRoundtripTest(unittest.TestCase):
    """Casos de prueba para RdkitRoundtripTest."""

    def test_fallback_smiles_keeps_carbonyl_branch_on_parent_atom(self):
        graph = MolGraph()
        n1 = graph.add_atom("N", 0.0, 0.0)
        carbonyl = graph.add_atom("C", 1.5, 0.0)
        oxygen = graph.add_atom("O", 3.0, 0.0)
        n2 = graph.add_atom("N", 1.5, -1.5)
        graph.add_bond(n1.id, carbonyl.id, order=1)
        graph.add_bond(carbonyl.id, oxygen.id, order=2)
        graph.add_bond(carbonyl.id, n2.id, order=1)

        original = rdkit_io._rdkit_available
        rdkit_io._rdkit_available = lambda: False
        try:
            smiles = molgraph_to_smiles(graph)
        finally:
            rdkit_io._rdkit_available = original

        self.assertEqual(smiles, "NC(=O)N")

    def test_smiles_export_uses_isolated_rdkit_worker_when_available(self):
        graph = MolGraph()
        a1 = graph.add_atom("C", 0.0, 0.0)
        a2 = graph.add_atom("C", 1.5, 0.0)
        graph.add_bond(a1.id, a2.id, order=1)

        original_available = rdkit_io._rdkit_available
        original_isolated = rdkit_safe.molgraph_to_smiles_isolated
        rdkit_io._rdkit_available = lambda: True
        rdkit_safe.molgraph_to_smiles_isolated = lambda _graph: ("worker-smiles", None)
        try:
            smiles = rdkit_io.molgraph_to_smiles(graph)
        finally:
            rdkit_io._rdkit_available = original_available
            rdkit_safe.molgraph_to_smiles_isolated = original_isolated

        self.assertEqual(smiles, "worker-smiles")

    def test_smiles_export_falls_back_when_isolated_worker_fails(self):
        graph = MolGraph()
        a1 = graph.add_atom("C", 0.0, 0.0)
        a2 = graph.add_atom("C", 1.5, 0.0)
        graph.add_bond(a1.id, a2.id, order=1)

        original_available = rdkit_io._rdkit_available
        original_isolated = rdkit_safe.molgraph_to_smiles_isolated
        rdkit_io._rdkit_available = lambda: True
        rdkit_safe.molgraph_to_smiles_isolated = lambda _graph: (None, "timeout")
        try:
            smiles = rdkit_io.molgraph_to_smiles(graph)
        finally:
            rdkit_io._rdkit_available = original_available
            rdkit_safe.molgraph_to_smiles_isolated = original_isolated

        self.assertIn(smiles, {"CC", "C-C"})

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
    def test_smiles_import_regression_for_naphthyl_urea(self):
        smiles = "O=C(NCCCl)Nc1cccc2ccccc12"
        expected = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), canonical=True)

        graph = smiles_to_molgraph(smiles)

        self.assertGreater(len(graph.atoms), 0)
        self.assertEqual(molgraph_to_smiles(graph), expected)

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

    @unittest.skipIf(not RDKit_AVAILABLE, "RDKit no disponible")
    def test_pseudoatoms_use_dummy_atoms_in_rdkit_conversion(self):
        graph = MolGraph()
        a1 = graph.add_atom("OH", 0.0, 0.0, is_explicit=True)
        a2 = graph.add_atom("CH2OH", 1.5, 0.0, is_explicit=True)
        graph.add_bond(a1.id, a2.id, order=1)

        mol, id_map = molgraph_to_rdkit_with_map(graph)
        self.assertEqual(len(id_map), 2)
        self.assertEqual(mol.GetNumAtoms(), 2)
        self.assertEqual(mol.GetAtomWithIdx(id_map[a1.id]).GetAtomicNum(), 0)
        self.assertEqual(mol.GetAtomWithIdx(id_map[a2.id]).GetAtomicNum(), 0)

    def test_kekulize_display_orders_returns_none_for_pseudoatoms(self):
        graph = MolGraph()
        a1 = graph.add_atom("OH", 0.0, 0.0, is_explicit=True)
        a2 = graph.add_atom("C", 1.5, 0.0, is_explicit=True)
        graph.add_bond(a1.id, a2.id, order=1, is_aromatic=True)
        self.assertIsNone(kekulize_display_orders(graph))

    def test_molfile_duplicate_bond_lines_are_deduped(self):
        molfile = (
            "Chemuson\n"
            "Chemuson\n"
            "\n"
            "  2  2  0  0  0  0  0  0  0 0999 V2000\n"
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "  1  2  1  0  0  0  0\n"
            "  1  2  2  0  0  0  0\n"
            "M  END\n"
        )
        graph = molfile_to_molgraph(molfile)
        self.assertEqual(len(graph.atoms), 2)
        self.assertEqual(len(graph.bonds), 1)
        only_bond = next(iter(graph.bonds.values()))
        self.assertEqual(only_bond.order, 2)


if __name__ == "__main__":
    unittest.main()
