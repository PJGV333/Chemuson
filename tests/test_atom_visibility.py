"""Pruebas unitarias para test_atom_visibility."""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from PyQt6.QtWidgets import QApplication

from chemuson.core.model import Atom
from chemuson.gui.items import AtomItem


class AtomVisibilityTest(unittest.TestCase):
    """Casos de prueba para AtomVisibilityTest."""
    @classmethod
    def setUpClass(cls) -> None:
        """Función de prueba auxiliar para setUpClass.

        Returns:
            None.

        """
        cls._app = QApplication.instance() or QApplication([])

    def test_carbon_hidden_by_default(self):
        """Verifica carbon hidden by default.

        Returns:
            None.

        """
        atom = Atom(id=1, element="C", x=0, y=0)
        item = AtomItem(atom, show_carbon=False, show_hydrogen=True)
        self.assertFalse(item.label.isVisible())

    def test_nitrogen_always_visible(self):
        """Verifica nitrogen always visible.

        Returns:
            None.

        """
        atom = Atom(id=1, element="N", x=0, y=0)
        item = AtomItem(atom, show_carbon=False, show_hydrogen=False)
        self.assertTrue(item.label.isVisible())

    def test_coordination_center_label_is_centered(self):
        """La etiqueta del centro de coordinación se centra sobre la esfera."""
        atom = Atom(id=1, element="Pd", x=0, y=0, is_coordination_center=True)
        item = AtomItem(atom, show_carbon=False, show_hydrogen=False)
        rect = item.label.boundingRect()
        pos = item.label.pos()
        cx = pos.x() + rect.width() / 2.0
        cy = pos.y() + rect.height() / 2.0
        self.assertAlmostEqual(cx, 0.0, delta=0.2)
        self.assertAlmostEqual(cy, 0.0, delta=0.2)

    def test_coordination_toggle_keeps_label_centered(self):
        """Al activar centro de coordinación, la etiqueta permanece centrada."""
        atom = Atom(id=1, element="Pd", x=0, y=0, is_coordination_center=False)
        item = AtomItem(atom, show_carbon=False, show_hydrogen=False)
        item.set_coordination_center(True)
        rect = item.label.boundingRect()
        pos = item.label.pos()
        cx = pos.x() + rect.width() / 2.0
        cy = pos.y() + rect.height() / 2.0
        self.assertAlmostEqual(cx, 0.0, delta=0.2)
        self.assertAlmostEqual(cy, 0.0, delta=0.2)

    def test_charge_label_uses_superscript_notation(self):
        """La carga se muestra con superíndice Unicode (ej. ²⁺)."""
        atom = Atom(id=1, element="N", x=0, y=0, charge=2)
        item = AtomItem(atom, show_carbon=False, show_hydrogen=False)
        self.assertTrue(item.charge_label.isVisible())
        self.assertEqual(item.charge_label.toPlainText(), "²⁺")
