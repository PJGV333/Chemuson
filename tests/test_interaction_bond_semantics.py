"""Regresiones para enlaces de interacción intermolecular."""

import os
import sys

import pytest
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.core.model import BondStyle
from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.commands import AddBondCommand


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def test_interaction_bond_keeps_explicit_carbon_labels() -> None:
    """El enlace intermolecular no debe convertir C explícitos en implícitos."""
    canvas = ChemusonCanvas()
    canvas.state.show_implicit_carbons = False

    left = canvas.model.add_atom("C", 20.0, 20.0, is_explicit=True).id
    right = canvas.model.add_atom("C", 60.0, 20.0, is_explicit=True).id
    canvas._rebuild_items_from_model()

    canvas.undo_stack.push(
        AddBondCommand(
            canvas.model,
            canvas,
            left,
            right,
            style=BondStyle.INTERACTION,
        )
    )

    assert canvas.model.get_atom(left).is_explicit is True
    assert canvas.model.get_atom(right).is_explicit is True
    assert left not in canvas.validate_structure()
    assert right not in canvas.validate_structure()
