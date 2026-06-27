"""Pruebas de nomenclatura de coordinación (MVP experimental)."""

import unittest


from chemuson.core.model import BondStyle, MolGraph
from chemuson.chemname import NameOptions, iupac_name


class ChemNameCoordinationTest(unittest.TestCase):
    """Cobertura mínima de coordinación: carbonyl, ammine/halo, Cp."""

    def test_nickel_tetracarbonyl(self):
        graph = MolGraph()
        ni = graph.add_atom("Ni", 0.0, 0.0, is_coordination_center=True)
        for idx in range(4):
            c = graph.add_atom("C", float(idx + 1), 0.0)
            o = graph.add_atom("O", float(idx + 1), 1.0)
            graph.add_bond(ni.id, c.id, style=BondStyle.COORDINATION, donor_atom_id=c.id)
            graph.add_bond(c.id, o.id, order=3)
        self.assertEqual(
            iupac_name(graph, NameOptions(allow_coordination=True)),
            "tetracarbonylnickel(0)",
        )

    def test_cis_and_trans_diamminedichloroplatinum(self):
        graph = MolGraph()
        pt = graph.add_atom("Pt", 0.0, 0.0, is_coordination_center=True)
        n1 = graph.add_atom("N", 1.0, 0.0, explicit_h=3)
        n2 = graph.add_atom("N", 0.0, 1.0, explicit_h=3)
        cl1 = graph.add_atom("Cl", -1.0, 0.0, charge=-1)
        cl2 = graph.add_atom("Cl", 0.0, -1.0, charge=-1)
        graph.add_bond(pt.id, n1.id, style=BondStyle.COORDINATION, donor_atom_id=n1.id)
        graph.add_bond(pt.id, n2.id, style=BondStyle.COORDINATION, donor_atom_id=n2.id)
        graph.add_bond(pt.id, cl1.id, style=BondStyle.COORDINATION, donor_atom_id=cl1.id)
        graph.add_bond(pt.id, cl2.id, style=BondStyle.COORDINATION, donor_atom_id=cl2.id)
        self.assertEqual(
            iupac_name(graph, NameOptions(allow_coordination=True)),
            "cis-diamminedichloroplatinum(II)",
        )

        graph = MolGraph()
        pt = graph.add_atom("Pt", 0.0, 0.0, is_coordination_center=True)
        n1 = graph.add_atom("N", 1.0, 0.0, explicit_h=3)
        n2 = graph.add_atom("N", -1.0, 0.0, explicit_h=3)
        cl1 = graph.add_atom("Cl", 0.0, 1.0, charge=-1)
        cl2 = graph.add_atom("Cl", 0.0, -1.0, charge=-1)
        graph.add_bond(pt.id, n1.id, style=BondStyle.COORDINATION, donor_atom_id=n1.id)
        graph.add_bond(pt.id, n2.id, style=BondStyle.COORDINATION, donor_atom_id=n2.id)
        graph.add_bond(pt.id, cl1.id, style=BondStyle.COORDINATION, donor_atom_id=cl1.id)
        graph.add_bond(pt.id, cl2.id, style=BondStyle.COORDINATION, donor_atom_id=cl2.id)
        self.assertEqual(
            iupac_name(graph, NameOptions(allow_coordination=True)),
            "trans-diamminedichloroplatinum(II)",
        )

    def test_ferrocene_eta5(self):
        graph = MolGraph()
        fe = graph.add_atom("Fe", 0.0, 0.0, is_coordination_center=True)
        ring_a = [graph.add_atom("C", float(i), 2.0) for i in range(5)]
        ring_b = [graph.add_atom("C", float(i), -2.0) for i in range(5)]
        for ring in (ring_a, ring_b):
            for i in range(5):
                graph.add_bond(ring[i].id, ring[(i + 1) % 5].id, order=1, is_aromatic=True)
                graph.add_bond(
                    fe.id,
                    ring[i].id,
                    style=BondStyle.COORDINATION,
                    donor_atom_id=ring[i].id,
                )
        self.assertEqual(
            iupac_name(graph, NameOptions(allow_coordination=True)),
            "bis(η5-cyclopentadienyl)iron(II)",
        )


if __name__ == "__main__":
    unittest.main()
