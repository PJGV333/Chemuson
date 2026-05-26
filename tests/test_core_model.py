"""Pruebas unitarias para test_core_model."""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.core.model import BondStyle, BondStereo, MolGraph
from chemuson.chemname.molview import MolView


class MolGraphTest(unittest.TestCase):
    """Casos de prueba para MolGraphTest."""
    def test_add_atom_and_bond(self):
        """Verifica add atom and bond.

        Returns:
            None.

        """
        graph = MolGraph()
        a1 = graph.add_atom("C", 0.0, 0.0)
        a2 = graph.add_atom("O", 1.0, 0.0)
        bond = graph.add_bond(a1.id, a2.id, order=2, style=BondStyle.PLAIN, stereo=BondStereo.NONE)

        self.assertEqual(len(graph.atoms), 2)
        self.assertEqual(len(graph.bonds), 1)
        self.assertEqual(bond.order, 2)

    def test_update_bond(self):
        """Verifica update bond.

        Returns:
            None.

        """
        graph = MolGraph()
        a1 = graph.add_atom("C", 0.0, 0.0)
        a2 = graph.add_atom("C", 1.0, 0.0)
        bond = graph.add_bond(a1.id, a2.id, order=1, style=BondStyle.PLAIN, stereo=BondStereo.NONE)

        graph.update_bond(bond.id, order=3, style=BondStyle.WAVY)
        updated = graph.get_bond(bond.id)
        self.assertEqual(updated.order, 3)
        self.assertEqual(updated.style, BondStyle.WAVY)

    def test_update_bond_to_flexible_style(self):
        """Verifica que el estilo flexible esté disponible en el modelo."""
        graph = MolGraph()
        a1 = graph.add_atom("C", 0.0, 0.0)
        a2 = graph.add_atom("C", 1.0, 0.0)
        bond = graph.add_bond(a1.id, a2.id, order=1, style=BondStyle.PLAIN, stereo=BondStereo.NONE)

        graph.update_bond(bond.id, style=BondStyle.FLEX)
        updated = graph.get_bond(bond.id)
        self.assertEqual(updated.style, BondStyle.FLEX)
        self.assertEqual(updated.order, 1)

    def test_flexible_bond_curves_reset_when_style_changes(self):
        """Al salir de FLEX, las curvas almacenadas deben limpiarse."""
        graph = MolGraph()
        a1 = graph.add_atom("C", 0.0, 0.0)
        a2 = graph.add_atom("C", 1.0, 0.0)
        bond = graph.add_bond(
            a1.id,
            a2.id,
            order=1,
            style=BondStyle.FLEX,
            stereo=BondStereo.NONE,
            flex_curve_1=0.3,
            flex_curve_2=-0.25,
        )
        self.assertEqual(bond.flex_curve_1, 0.3)
        self.assertEqual(bond.flex_curve_2, -0.25)

        graph.update_bond(bond.id, style=BondStyle.PLAIN)
        updated = graph.get_bond(bond.id)
        self.assertEqual(updated.style, BondStyle.PLAIN)
        self.assertIsNone(updated.flex_curve_1)
        self.assertIsNone(updated.flex_curve_2)

    def test_validate_allows_hypervalent_phosphorus(self):
        """Verifica validate allows hypervalent phosphorus.

        Returns:
            None.

        """
        graph = MolGraph()
        p = graph.add_atom("P", 0.0, 0.0)
        o_dbl = graph.add_atom("O", 0.0, 1.0)
        o1 = graph.add_atom("O", -1.0, 0.0)
        o2 = graph.add_atom("O", -1.0, -1.0)
        c = graph.add_atom("C", 1.0, 0.0)

        graph.add_bond(p.id, o_dbl.id, order=2)
        graph.add_bond(p.id, o1.id, order=1)
        graph.add_bond(p.id, o2.id, order=1)
        graph.add_bond(p.id, c.id, order=1)

        self.assertNotIn(p.id, graph.validate())

    def test_validate_allows_sf6(self):
        """Verifica validate allows sf6.

        Returns:
            None.

        """
        graph = MolGraph()
        s = graph.add_atom("S", 0.0, 0.0)
        fs = [graph.add_atom("F", float(i), 1.0) for i in range(6)]
        for f in fs:
            graph.add_bond(s.id, f.id, order=1)
        self.assertNotIn(s.id, graph.validate())

    def test_validate_allows_if7(self):
        """Verifica validate allows if7.

        Returns:
            None.

        """
        graph = MolGraph()
        i = graph.add_atom("I", 0.0, 0.0)
        fs = [graph.add_atom("F", float(n), 1.0) for n in range(7)]
        for f in fs:
            graph.add_bond(i.id, f.id, order=1)
        self.assertNotIn(i.id, graph.validate())

    def test_validate_flags_neutral_ammonium(self):
        """N neutro tetravalente debe marcarse como error de valencia."""

        graph = MolGraph()
        n = graph.add_atom("N", 0.0, 0.0)
        hs = [graph.add_atom("H", float(k), 1.0) for k in range(4)]
        for h in hs:
            graph.add_bond(n.id, h.id, order=1)
        self.assertIn(n.id, graph.validate())

    def test_validate_detailed_reports_interactive_valence_issue(self):
        """El diagnóstico detallado expone código, severidad, hint y matemática."""
        graph = MolGraph()
        n = graph.add_atom("N", 0.0, 0.0)
        hs = [graph.add_atom("H", float(k), 1.0) for k in range(4)]
        for h in hs:
            graph.add_bond(n.id, h.id, order=1)

        issues = graph.validate_detailed()
        issue = issues[n.id]

        self.assertEqual(issue.code, "VALENCE_EXCEEDED")
        self.assertEqual(issue.severity, "error")
        self.assertEqual(issue.target_type, "atom")
        self.assertIn("enlaces 4.00", issue.message)
        self.assertIn("carga formal", issue.hint)
        self.assertEqual(issue.allowed_valences, (3, 5))
        self.assertEqual(issue.observed_valence, 4.0)
        self.assertEqual(issue.bond_order_sum, 4.0)
        payload = issue.as_dict()
        self.assertEqual(payload["atom_id"], n.id)
        self.assertEqual(payload["code"], "VALENCE_EXCEEDED")
        self.assertIn("Sugerencia:", issue.tooltip_text())

    def test_validate_allows_tetramethylammonium_with_positive_charge(self):
        """Verifica que una carga positiva permita valencia tetravalente en N."""
        graph = MolGraph()
        n = graph.add_atom("N", 0.0, 0.0, charge=1)
        carbons = [graph.add_atom("C", float(k), 1.0) for k in range(4)]
        for carbon in carbons:
            graph.add_bond(n.id, carbon.id, order=1)
        self.assertNotIn(n.id, graph.validate())

    def test_validate_allows_neutral_amine_with_drawn_hydrogen(self):
        """Un N con dos enlaces pesados y un H explícito sigue siendo trivalente."""
        graph = MolGraph()
        n = graph.add_atom("N", 0.0, 0.0)
        c1 = graph.add_atom("C", -1.0, 0.0)
        c2 = graph.add_atom("C", 1.0, 0.0)
        h = graph.add_atom("H", 0.0, 1.0, is_explicit=True)

        graph.add_bond(n.id, c1.id, order=1)
        graph.add_bond(n.id, c2.id, order=1)
        graph.add_bond(n.id, h.id, order=1)

        self.assertNotIn(n.id, graph.validate())

    def test_implicit_h_count_does_not_double_count_drawn_hydrogen(self):
        """Un H explícito enlazado no debe impedir el H implícito restante en una amina."""
        graph = MolGraph()
        n = graph.add_atom("N", 0.0, 0.0)
        c = graph.add_atom("C", -1.0, 0.0)
        h = graph.add_atom("H", 1.0, 0.0, is_explicit=True)

        graph.add_bond(n.id, c.id, order=1)
        graph.add_bond(n.id, h.id, order=1)

        self.assertEqual(graph.implicit_h_count(n.id), 1)
        self.assertNotIn(n.id, graph.validate())

    def test_implicit_h_count_keeps_pyridine_like_aromatic_n_neutral(self):
        """El N aromático de un anillo de seis miembros no debe recibir H implícito."""
        graph = MolGraph()
        atoms = [graph.add_atom("N" if idx == 0 else "C", float(idx), 0.0) for idx in range(6)]
        for idx in range(6):
            graph.add_bond(
                atoms[idx].id,
                atoms[(idx + 1) % 6].id,
                order=1,
                is_aromatic=True,
            )

        self.assertEqual(graph.implicit_h_count(atoms[0].id), 0)
        self.assertEqual(graph.implicit_h_count(atoms[1].id), 1)

    def test_implicit_h_count_preserves_pyrrolic_default_in_five_member_ring(self):
        """El ajuste piridínico no debe eliminar el H implícito de un anillo de cinco miembros."""
        graph = MolGraph()
        atoms = [graph.add_atom("N" if idx == 0 else "C", float(idx), 0.0) for idx in range(5)]
        for idx in range(5):
            graph.add_bond(
                atoms[idx].id,
                atoms[(idx + 1) % 5].id,
                order=1,
                is_aromatic=True,
            )

        self.assertEqual(graph.implicit_h_count(atoms[0].id), 1)

    def test_validate_flags_pentavalent_carbon_when_charged(self):
        """C+ pentavalente sigue siendo inválido con regla iso-electrónica."""
        graph = MolGraph()
        c = graph.add_atom("C", 0.0, 0.0, charge=1)
        hs = [graph.add_atom("H", float(k), 1.0) for k in range(5)]
        for hydrogen in hs:
            graph.add_bond(c.id, hydrogen.id, order=1)
        self.assertIn(c.id, graph.validate())

    def test_validate_still_flags_overvalent_carbon(self):
        """Verifica validate still flags overvalent carbon.

        Returns:
            None.

        """
        graph = MolGraph()
        c = graph.add_atom("C", 0.0, 0.0)
        hs = [graph.add_atom("H", float(k), 1.0) for k in range(5)]
        for h in hs:
            graph.add_bond(c.id, h.id, order=1)
        self.assertIn(c.id, graph.validate())

    def test_coordination_bond_preserves_donor(self):
        """Verifica que el donador se guarda en enlaces coordinativos."""
        graph = MolGraph()
        metal = graph.add_atom("Pd", 0.0, 0.0, is_coordination_center=True)
        donor = graph.add_atom("N", 1.5, 0.0)
        self.assertEqual(metal.sphere_radius, 16.0)
        self.assertTrue(metal.sphere_filled)
        self.assertFalse(metal.sphere_transparent)
        self.assertEqual(metal.sphere_color, "#D9DDE3")
        bond = graph.add_bond(
            donor.id,
            metal.id,
            order=1,
            style=BondStyle.COORDINATION,
            donor_atom_id=donor.id,
        )
        self.assertEqual(bond.style, BondStyle.COORDINATION)
        self.assertEqual(bond.donor_atom_id, donor.id)

        graph.update_bond(bond.id, style=BondStyle.PLAIN)
        self.assertIsNone(graph.get_bond(bond.id).donor_atom_id)

    def test_validate_ignores_coordination_bonds(self):
        """Verifica que validate no cuente enlaces coordinativos para valencia."""
        graph = MolGraph()
        carbon = graph.add_atom("C", 0.0, 0.0)
        for idx in range(6):
            metal = graph.add_atom("Pd", float(idx + 1), 0.0, is_coordination_center=True)
            graph.add_bond(
                carbon.id,
                metal.id,
                order=1,
                style=BondStyle.COORDINATION,
                donor_atom_id=carbon.id,
            )
        self.assertNotIn(carbon.id, graph.validate())

    def test_validate_ignores_interaction_bonds_on_hydrogen(self):
        """Un enlace intermolecular no debe alterar la valencia real de un H."""
        graph = MolGraph()
        nitrogen = graph.add_atom("N", 0.0, 0.0)
        hydrogen = graph.add_atom("H", 1.0, 0.0, is_explicit=True)
        chloride = graph.add_atom("Cl", 2.0, 0.0, charge=-1)

        graph.add_bond(nitrogen.id, hydrogen.id, order=1)
        implicit_before = graph.implicit_h_count(nitrogen.id)
        graph.add_bond(hydrogen.id, chloride.id, style=BondStyle.INTERACTION)

        errors = graph.validate()
        self.assertNotIn(nitrogen.id, errors)
        self.assertNotIn(hydrogen.id, errors)
        self.assertNotIn(chloride.id, errors)
        self.assertEqual(graph.implicit_h_count(nitrogen.id), implicit_before)

    def test_molview_ignores_interaction_bonds_in_connectivity(self):
        """La vista estructural no debe conectar átomos unidos solo por interacción."""
        graph = MolGraph()
        donor = graph.add_atom("O", 0.0, 0.0)
        acceptor = graph.add_atom("H", 1.0, 0.0, is_explicit=True)
        graph.add_bond(donor.id, acceptor.id, style=BondStyle.INTERACTION)

        view = MolView(graph)
        self.assertEqual(view.neighbors(donor.id), [])
        self.assertEqual(view.neighbors(acceptor.id), [])
        self.assertEqual(view.bonds(), [])
        self.assertEqual(view.bond_orders(), [])

    def test_formal_charge_property_and_total_charge(self):
        """La carga formal por átomo y total de molécula deben estar sincronizadas."""
        graph = MolGraph()
        n = graph.add_atom("N", 0.0, 0.0, formal_charge=1)
        cl = graph.add_atom("Cl", 1.0, 0.0, charge=-1)
        self.assertEqual(n.charge, 1)
        self.assertEqual(n.formal_charge, 1)
        self.assertEqual(cl.formal_charge, -1)
        self.assertEqual(graph.total_formal_charge(), 0)

    def test_no_implicit_disables_auto_hydrogens(self):
        """`no_implicit` debe forzar cero H implícitos."""
        graph = MolGraph()
        c1 = graph.add_atom("C", 0.0, 0.0, no_implicit=True)
        c2 = graph.add_atom("C", 1.0, 0.0)
        graph.add_bond(c1.id, c2.id, order=1)

        self.assertEqual(graph.implicit_h_count(c1.id), 0)
        errors = graph.validate()
        self.assertIn(c1.id, errors)


if __name__ == "__main__":
    unittest.main()
