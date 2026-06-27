"""Pruebas para ajuste de grosor en flechas con los atajos de grosor."""

import os

import pytest
from PyQt6.QtCore import QPointF, QRectF
from PyQt6.QtWidgets import QApplication


from chemuson.chemio.persistence import PersistenceManager
from chemuson.gui.canvas import ChemusonCanvas


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    return QApplication.instance() or QApplication([])


def test_increase_arrow_thickness_and_undo() -> None:
    canvas = ChemusonCanvas()
    arrow = canvas.add_arrow_item(QPointF(30.0, 40.0), QPointF(120.0, 40.0), kind="forward")
    arrow.setSelected(True)

    default_stroke = float(canvas.drawing_style.stroke_px)
    assert arrow.stroke_px() is None

    canvas.increase_selected_bond_thickness()

    assert arrow.stroke_px() is not None
    assert float(arrow.stroke_px()) > default_stroke

    canvas.undo_stack.undo()

    assert arrow.stroke_px() is None


def test_increase_bracket_thickness_and_undo() -> None:
    canvas = ChemusonCanvas()
    bracket = canvas.add_bracket_item(QRectF(30.0, 40.0, 80.0, 50.0), "[]")
    bracket.setSelected(True)

    default_stroke = float(canvas.drawing_style.stroke_px)
    assert bracket.stroke_px() is None

    canvas.increase_selected_bond_thickness()

    assert bracket.stroke_px() is not None
    assert float(bracket.stroke_px()) > default_stroke

    canvas.undo_stack.undo()

    assert bracket.stroke_px() is None


def test_shortcut_path_adjusts_selected_bonds_arrows_and_brackets() -> None:
    canvas = ChemusonCanvas()
    atom1 = canvas.model.add_atom("C", 20.0, 20.0)
    atom2 = canvas.model.add_atom("C", 80.0, 20.0)
    canvas.add_atom_item(atom1)
    canvas.add_atom_item(atom2)
    bond = canvas.model.add_bond(atom1.id, atom2.id)
    canvas.add_bond_item(bond)
    canvas.state.selected_bonds = {bond.id}

    arrow = canvas.add_arrow_item(QPointF(20.0, 80.0), QPointF(120.0, 80.0), kind="forward")
    arrow.setSelected(True)
    bracket = canvas.add_bracket_item(QRectF(10.0, 100.0, 130.0, 70.0), "[]")
    bracket.setSelected(True)
    # La selección de escena puede limpiar `state.selected_bonds`; lo reponemos.
    canvas.state.selected_bonds = {bond.id}

    canvas.increase_selected_bond_thickness()
    assert canvas.model.get_bond(bond.id).stroke_px is not None
    assert arrow.stroke_px() is not None
    assert bracket.stroke_px() is not None

    canvas.reset_selected_bond_thickness()
    assert canvas.model.get_bond(bond.id).stroke_px is None
    assert arrow.stroke_px() is None
    assert bracket.stroke_px() is None


def test_arrow_thickness_persists_after_save_load(tmp_path) -> None:
    source = ChemusonCanvas()
    arrow = source.add_arrow_item(QPointF(15.0, 15.0), QPointF(110.0, 15.0), kind="forward")
    arrow.setSelected(True)
    source.increase_selected_bond_thickness()
    expected_stroke = arrow.stroke_px()
    assert expected_stroke is not None

    file_path = tmp_path / "arrow_stroke.cmsn"
    PersistenceManager.save_to_file(str(file_path), source)

    restored = ChemusonCanvas()
    PersistenceManager.load_from_file(str(file_path), restored)

    assert len(restored.arrow_items) == 1
    restored_arrow = restored.arrow_items[0]
    assert restored_arrow.stroke_px() == pytest.approx(float(expected_stroke), rel=1e-6)
