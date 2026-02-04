import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from PyQt6.QtWidgets import QApplication

from core.model import Atom
from gui.items import AtomItem


class AtomVisibilityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._app = QApplication.instance() or QApplication([])

    def test_carbon_hidden_by_default(self):
        atom = Atom(id=1, element="C", x=0, y=0)
        item = AtomItem(atom, show_carbon=False, show_hydrogen=True)
        self.assertFalse(item.label.isVisible())

    def test_nitrogen_always_visible(self):
        atom = Atom(id=1, element="N", x=0, y=0)
        item = AtomItem(atom, show_carbon=False, show_hydrogen=False)
        self.assertTrue(item.label.isVisible())
