"""Pruebas para inversión manual de orientación en dobles enlaces."""

import os
import sys

import pytest
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.chemio.persistence import PersistenceManager
from chemuson.core.model import BondStyle
from chemuson.gui.canvas import ChemusonCanvas


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    return QApplication.instance() or QApplication([])


def _line_mids(item) -> tuple[tuple[float, float], tuple[float, float]]:
    """Devuelve punto medio de línea sigma y línea pi para enlaces dobles."""
    path = item.path()
    assert path.elementCount() >= 4
    e0 = path.elementAt(0)
    e1 = path.elementAt(1)
    e2 = path.elementAt(2)
    e3 = path.elementAt(3)
    sigma_mid = ((e0.x + e1.x) * 0.5, (e0.y + e1.y) * 0.5)
    pi_mid = ((e2.x + e3.x) * 0.5, (e2.y + e3.y) * 0.5)
    return sigma_mid, pi_mid


def _build_square_ring_with_top_double(canvas: ChemusonCanvas) -> int:
    """Construye un anillo cuadrado con doble enlace superior y retorna su ID."""
    a1 = canvas.model.add_atom("C", 100.0, 100.0)
    a2 = canvas.model.add_atom("C", 140.0, 100.0)
    a3 = canvas.model.add_atom("C", 140.0, 140.0)
    a4 = canvas.model.add_atom("C", 100.0, 140.0)

    top = canvas.model.add_bond(a1.id, a2.id, order=2, style=BondStyle.PLAIN, ring_id=1)
    canvas.model.add_bond(a2.id, a3.id, order=1, style=BondStyle.PLAIN, ring_id=1)
    canvas.model.add_bond(a3.id, a4.id, order=1, style=BondStyle.PLAIN, ring_id=1)
    canvas.model.add_bond(a4.id, a1.id, order=1, style=BondStyle.PLAIN, ring_id=1)
    canvas._rebuild_items_from_model()
    return top.id


def test_ring_double_defaults_inward_and_can_flip_outward() -> None:
    canvas = ChemusonCanvas()
    bond_id = _build_square_ring_with_top_double(canvas)

    item = canvas.bond_items[bond_id]
    sigma_mid, pi_mid = _line_mids(item)
    # Centro del anillo está hacia y positiva respecto al enlace superior.
    assert pi_mid[1] > sigma_mid[1]

    item.toggle_double_bond_orientation()
    sigma_mid_2, pi_mid_2 = _line_mids(item)
    assert pi_mid_2[1] < sigma_mid_2[1]
    assert item.manual_pi_offset in {-1, 1}


def test_orientation_persists_after_save_and_load(tmp_path) -> None:
    canvas = ChemusonCanvas()
    bond_id = _build_square_ring_with_top_double(canvas)
    item = canvas.bond_items[bond_id]

    item.setSelected(True)
    canvas._sync_selection_from_scene()
    assert canvas.toggle_selected_double_bond_orientation() is True
    assert canvas.undo_stack.isClean() is False

    sign_before = canvas.model.get_bond(bond_id).pi_offset_sign
    assert sign_before in {-1, 1}

    filepath = tmp_path / "orientation.cmsn"
    PersistenceManager.save_to_file(str(filepath), canvas)

    restored = ChemusonCanvas()
    PersistenceManager.load_from_file(str(filepath), restored)
    restored_bond = restored.model.get_bond(bond_id)
    assert restored_bond.pi_offset_sign == sign_before

    sigma_mid, pi_mid = _line_mids(restored.bond_items[bond_id])
    assert pi_mid[1] < sigma_mid[1]


def test_non_ring_double_keeps_screen_left_convention() -> None:
    canvas = ChemusonCanvas()
    a1 = canvas.model.add_atom("C", 200.0, 220.0)
    a2 = canvas.model.add_atom("C", 260.0, 180.0)
    bond = canvas.model.add_bond(a1.id, a2.id, order=2, style=BondStyle.PLAIN)
    canvas._rebuild_items_from_model()

    sigma_mid, pi_mid = _line_mids(canvas.bond_items[bond.id])
    assert pi_mid[0] < sigma_mid[0]
    assert canvas.model.get_bond(bond.id).pi_offset_sign is None


def test_toggle_shortcut_logic_does_not_apply_to_non_double_bonds() -> None:
    canvas = ChemusonCanvas()
    a1 = canvas.model.add_atom("C", 20.0, 20.0)
    a2 = canvas.model.add_atom("C", 70.0, 20.0)
    triple = canvas.model.add_bond(a1.id, a2.id, order=3, style=BondStyle.PLAIN)
    canvas._rebuild_items_from_model()

    item = canvas.bond_items[triple.id]
    item.setSelected(True)
    canvas._sync_selection_from_scene()

    assert canvas.toggle_selected_double_bond_orientation() is False
    assert canvas.model.get_bond(triple.id).pi_offset_sign is None
