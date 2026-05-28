"""Pruebas unitarias para test_rdkit_roundtrip."""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.core.model import MolGraph
from chemuson.chemio import rdkit_io, rdkit_safe
from chemuson.chemio.rdkit_io import (
    expand_abbreviations_for_calculation,
    kekulize_display_orders,
    molfile_to_molgraph,
    molgraph_to_molfile,
    molgraph_to_smiles,
    molgraph_to_rdkit_with_map,
    smiles_to_molgraph,
)

RDKit_AVAILABLE = rdkit_io._rdkit_available()


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

    def test_smiles_export_expands_methoxy_abbreviation_for_fallback(self):
        graph = MolGraph()
        a1 = graph.add_atom("C", 0.0, 0.0)
        a2 = graph.add_atom("OMe", 1.5, 0.0, is_explicit=True)
        graph.add_bond(a1.id, a2.id, order=1)

        original_available = rdkit_io._rdkit_available
        rdkit_io._rdkit_available = lambda: False
        try:
            smiles = molgraph_to_smiles(graph)
        finally:
            rdkit_io._rdkit_available = original_available

        self.assertEqual(smiles, "COC")
        self.assertNotIn("*", smiles)
        self.assertNotIn("OMe", smiles)
        self.assertEqual(len(graph.atoms), 2)
        self.assertEqual(graph.get_atom(a2.id).element, "OMe")

    def test_isolated_smiles_export_expands_methoxy_before_worker(self):
        graph = MolGraph()
        a1 = graph.add_atom("C", 0.0, 0.0)
        a2 = graph.add_atom("OMe", 1.5, 0.0, is_explicit=True)
        graph.add_bond(a1.id, a2.id, order=1)

        captured = {}
        original_isolated = rdkit_safe.molgraph_to_smiles_isolated

        def _fake_isolated(expanded_graph, timeout_s=5.0):
            captured["elements"] = sorted(atom.element for atom in expanded_graph.atoms.values())
            captured["bond_count"] = len(expanded_graph.bonds)
            captured["timeout_s"] = timeout_s
            return "COC", None

        rdkit_safe.molgraph_to_smiles_isolated = _fake_isolated
        try:
            smiles = rdkit_io.molgraph_to_smiles_isolated_or_error(graph, timeout_s=1.5)
        finally:
            rdkit_safe.molgraph_to_smiles_isolated = original_isolated

        self.assertEqual(smiles, "COC")
        self.assertEqual(captured["elements"], ["C", "C", "O"])
        self.assertEqual(captured["bond_count"], 2)
        self.assertEqual(captured["timeout_s"], 1.5)

    def test_superatom_expansion_supports_publication_abbreviations(self):
        graph = MolGraph()
        root = graph.add_atom("N", 0.0, 0.0)
        ac = graph.add_atom("Ac", 1.5, 0.0, is_explicit=True)
        boc = graph.add_atom("Boc", 3.0, 0.0, is_explicit=True)
        ts = graph.add_atom("Ts", 4.5, 0.0, is_explicit=True)
        graph.add_bond(root.id, ac.id, order=1)
        graph.add_bond(root.id, boc.id, order=1)
        graph.add_bond(root.id, ts.id, order=1)

        expanded = expand_abbreviations_for_calculation(graph)
        elements = [atom.element for atom in expanded.atoms.values()]

        self.assertNotIn("Ac", elements)
        self.assertNotIn("Boc", elements)
        self.assertNotIn("Ts", elements)
        self.assertEqual(elements.count("S"), 1)
        self.assertGreaterEqual(elements.count("O"), 5)
        self.assertGreaterEqual(elements.count("C"), 12)
        self.assertEqual(graph.get_atom(ac.id).element, "Ac")

    def test_smiles_export_merges_nearly_duplicate_ring_closure_atoms(self):
        graph = MolGraph()
        c1 = graph.add_atom("C", 0.0, 0.0)
        n = graph.add_atom("N", 40.0, 0.0, is_explicit=True)
        c2 = graph.add_atom("C", 80.0, 0.0)
        duplicate_c1 = graph.add_atom("C", 0.5, 0.5)
        graph.add_bond(c1.id, n.id, order=1)
        graph.add_bond(n.id, c2.id, order=1)
        graph.add_bond(c2.id, duplicate_c1.id, order=1)

        original_available = rdkit_io._rdkit_available
        rdkit_io._rdkit_available = lambda: False
        try:
            smiles = molgraph_to_smiles(graph)
        finally:
            rdkit_io._rdkit_available = original_available

        self.assertEqual(smiles, "C1NC1")
        self.assertEqual(len(graph.atoms), 4)
        self.assertNotEqual(c1.id, duplicate_c1.id)

    def test_isolated_smiles_export_merges_duplicates_before_worker(self):
        graph = MolGraph()
        c1 = graph.add_atom("C", 0.0, 0.0)
        n = graph.add_atom("N", 40.0, 0.0, is_explicit=True)
        c2 = graph.add_atom("C", 80.0, 0.0)
        duplicate_c1 = graph.add_atom("C", 0.5, 0.5)
        graph.add_bond(c1.id, n.id, order=1)
        graph.add_bond(n.id, c2.id, order=1)
        graph.add_bond(c2.id, duplicate_c1.id, order=1)

        captured = {}
        original_isolated = rdkit_safe.molgraph_to_smiles_isolated

        def _fake_isolated(export_graph, timeout_s=5.0):
            captured["atom_count"] = len(export_graph.atoms)
            captured["bond_pairs"] = sorted(
                tuple(sorted((bond.a1_id, bond.a2_id)))
                for bond in export_graph.bonds.values()
            )
            return "C1NC1", None

        rdkit_safe.molgraph_to_smiles_isolated = _fake_isolated
        try:
            smiles = rdkit_io.molgraph_to_smiles_isolated_or_error(graph, timeout_s=1.5)
        finally:
            rdkit_safe.molgraph_to_smiles_isolated = original_isolated

        self.assertEqual(smiles, "C1NC1")
        self.assertEqual(captured["atom_count"], 3)
        self.assertEqual(captured["bond_pairs"], [(1, 2), (1, 3), (2, 3)])

    def test_smiles_export_does_not_merge_different_elements_by_distance(self):
        graph = MolGraph()
        c1 = graph.add_atom("C", 0.0, 0.0)
        n = graph.add_atom("N", 40.0, 0.0, is_explicit=True)
        c2 = graph.add_atom("C", 80.0, 0.0)
        oxygen = graph.add_atom("O", 0.5, 0.5, is_explicit=True)
        graph.add_bond(c1.id, n.id, order=1)
        graph.add_bond(n.id, c2.id, order=1)
        graph.add_bond(c2.id, oxygen.id, order=1)

        original_available = rdkit_io._rdkit_available
        rdkit_io._rdkit_available = lambda: False
        try:
            smiles = molgraph_to_smiles(graph)
        finally:
            rdkit_io._rdkit_available = original_available

        self.assertEqual(smiles, "CNCO")

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
        from rdkit import Chem

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
