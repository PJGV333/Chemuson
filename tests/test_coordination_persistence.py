"""Pruebas de persistencia para centros y enlaces de coordinación."""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import BondStyle, MolGraph
from chemio.persistence import PersistenceManager


class _FakeCanvas:
    """Canvas mínimo para probar PersistenceManager sin GUI real."""

    def __init__(self) -> None:
        self.model = MolGraph()
        self._loaded_canvas_data = None
        self._rebuilt = False

    def get_persistence_data(self):
        return {"ok": True}

    def load_persistence_data(self, data):
        self._loaded_canvas_data = data

    def _rebuild_items_from_model(self):
        self._rebuilt = True


class CoordinationPersistenceTest(unittest.TestCase):
    """Casos de roundtrip de persistencia para coordinación."""

    def test_roundtrip_coordination_fields(self):
        canvas = _FakeCanvas()
        metal = canvas.model.add_atom(
            "Pd",
            0.0,
            0.0,
            is_coordination_center=True,
            sphere_radius=22.0,
            sphere_color="#00AA55",
            sphere_filled=False,
            sphere_transparent=True,
        )
        donor = canvas.model.add_atom("N", 40.0, 0.0)
        canvas.model.add_bond(
            donor.id,
            metal.id,
            style=BondStyle.COORDINATION,
            donor_atom_id=donor.id,
        )

        data = PersistenceManager.save_to_dict(canvas)
        self.assertEqual(data["model"]["bonds"][0]["style"], BondStyle.COORDINATION.value)
        self.assertEqual(data["model"]["bonds"][0]["type"], BondStyle.COORDINATION.value)
        restored = _FakeCanvas()
        PersistenceManager.load_from_dict(data, restored)

        self.assertTrue(restored._rebuilt)
        self.assertEqual(restored._loaded_canvas_data, {"ok": True})
        self.assertEqual(len(restored.model.atoms), 2)
        self.assertEqual(len(restored.model.bonds), 1)

        restored_metal = restored.model.get_atom(metal.id)
        self.assertTrue(restored_metal.is_coordination_center)
        self.assertEqual(restored_metal.sphere_radius, 22.0)
        self.assertEqual(restored_metal.sphere_color, "#00AA55")
        self.assertFalse(restored_metal.sphere_filled)
        self.assertTrue(restored_metal.sphere_transparent)

        restored_bond = next(iter(restored.model.bonds.values()))
        self.assertEqual(restored_bond.style, BondStyle.COORDINATION)
        self.assertEqual(restored_bond.donor_atom_id, donor.id)

    def test_load_defaults_plain_when_style_missing(self):
        canvas = _FakeCanvas()
        data = {
            "application": "Chemuson",
            "version": PersistenceManager.VERSION,
            "model": {
                "atoms": [
                    {"id": 1, "element": "C", "x": 0.0, "y": 0.0},
                    {"id": 2, "element": "O", "x": 40.0, "y": 0.0},
                ],
                "bonds": [
                    {
                        "id": 1,
                        "a1_id": 1,
                        "a2_id": 2,
                        "order": 1,
                    }
                ],
            },
            "canvas": {"ok": True},
        }
        PersistenceManager.load_from_dict(data, canvas)

        bond = canvas.model.get_bond(1)
        self.assertEqual(bond.style, BondStyle.PLAIN)
        self.assertIsNone(bond.donor_atom_id)

    def test_load_coordination_without_donor_infers_non_center_as_donor(self):
        canvas = _FakeCanvas()
        data = {
            "application": "Chemuson",
            "version": PersistenceManager.VERSION,
            "model": {
                "atoms": [
                    {"id": 1, "element": "N", "x": 40.0, "y": 0.0},
                    {
                        "id": 2,
                        "element": "Pd",
                        "x": 0.0,
                        "y": 0.0,
                        "is_coordination_center": True,
                    },
                ],
                "bonds": [
                    {
                        "id": 1,
                        "a1_id": 1,
                        "a2_id": 2,
                        "order": 1,
                        "type": BondStyle.COORDINATION.value,
                    }
                ],
            },
            "canvas": {"ok": True},
        }
        PersistenceManager.load_from_dict(data, canvas)

        bond = canvas.model.get_bond(1)
        self.assertEqual(bond.style, BondStyle.COORDINATION)
        self.assertEqual(bond.donor_atom_id, 1)

    def test_roundtrip_flexible_bond_style(self):
        canvas = _FakeCanvas()
        a1 = canvas.model.add_atom("O", 0.0, 0.0)
        a2 = canvas.model.add_atom("O", 90.0, 0.0)
        canvas.model.add_bond(
            a1.id,
            a2.id,
            style=BondStyle.FLEX,
            order=1,
            flex_curve_1=0.35,
            flex_curve_2=-0.28,
        )

        data = PersistenceManager.save_to_dict(canvas)
        self.assertEqual(data["model"]["bonds"][0]["style"], BondStyle.FLEX.value)
        self.assertEqual(data["model"]["bonds"][0]["type"], BondStyle.FLEX.value)
        self.assertEqual(data["model"]["bonds"][0]["flex_curve_1"], 0.35)
        self.assertEqual(data["model"]["bonds"][0]["flex_curve_2"], -0.28)

        restored = _FakeCanvas()
        PersistenceManager.load_from_dict(data, restored)
        bond = next(iter(restored.model.bonds.values()))
        self.assertEqual(bond.style, BondStyle.FLEX)
        self.assertIsNone(bond.donor_atom_id)
        self.assertEqual(bond.flex_curve_1, 0.35)
        self.assertEqual(bond.flex_curve_2, -0.28)


if __name__ == "__main__":
    unittest.main()
