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


if __name__ == "__main__":
    unittest.main()
