"""Pruebas para validación/normalización estricta con RDKit."""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.core.model import BondStyle, MolGraph
from chemuson.chemio.rdkit_io import molgraph_to_rdkit_with_map, strict_validate_and_normalize

try:
    from rdkit import Chem

    RDKit_AVAILABLE = True
except Exception:
    RDKit_AVAILABLE = False


class StrictRdkitValidationTest(unittest.TestCase):
    """Cobertura de chequeos estrictos para iones/complejos."""

    @unittest.skipIf(not RDKit_AVAILABLE, "RDKit no disponible")
    def test_strict_validation_accepts_tetramethylammonium(self):
        """N tetravalente con carga positiva debe validarse."""
        graph = MolGraph()
        n = graph.add_atom("N", 0.0, 0.0, charge=1)
        carbons = [graph.add_atom("C", float(i + 1), 0.0) for i in range(4)]
        for carbon in carbons:
            graph.add_bond(n.id, carbon.id, order=1)

        result = strict_validate_and_normalize(graph)

        self.assertTrue(result.is_valid)
        self.assertEqual(result.errors, [])
        self.assertIsNotNone(result.normalized_smiles)
        self.assertIn("[N+]", result.normalized_smiles)
        self.assertIsNotNone(result.normalized_graph)

    @unittest.skipIf(not RDKit_AVAILABLE, "RDKit no disponible")
    def test_strict_validation_rejects_overvalent_neutral_carbon(self):
        """Un C neutro pentavalente debe fallar sanitización estricta."""
        graph = MolGraph()
        c = graph.add_atom("C", 0.0, 0.0)
        hydrogens = [graph.add_atom("H", float(i + 1), 0.0) for i in range(5)]
        for h in hydrogens:
            graph.add_bond(c.id, h.id, order=1)

        result = strict_validate_and_normalize(graph)

        self.assertFalse(result.is_valid)
        self.assertTrue(result.errors)

    @unittest.skipIf(not RDKit_AVAILABLE, "RDKit no disponible")
    def test_coordination_bond_is_exported_as_dative(self):
        """Los enlaces de coordinación se exportan a RDKit como dativos."""
        graph = MolGraph()
        donor = graph.add_atom("N", 0.0, 0.0)
        metal = graph.add_atom("Pd", 1.5, 0.0, is_coordination_center=True)
        graph.add_bond(
            donor.id,
            metal.id,
            order=1,
            style=BondStyle.COORDINATION,
            donor_atom_id=donor.id,
        )

        mol, id_map = molgraph_to_rdkit_with_map(graph)
        bond = mol.GetBondBetweenAtoms(id_map[donor.id], id_map[metal.id])

        self.assertIsNotNone(bond)
        self.assertEqual(bond.GetBondType(), Chem.BondType.DATIVE)


if __name__ == "__main__":
    unittest.main()
