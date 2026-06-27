"""Regresiones para sustituciÃ³n de hidrÃ³genos implÃ­citos dibujados."""


import pytest
from PyQt6.QtWidgets import QApplication


from chemuson.core.model import BondStereo, BondStyle
from chemuson.gui.canvas import ChemusonCanvas


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def test_substitute_drawn_implicit_hydrogen_creates_real_substituent():
    """Clic sobre H implÃ­cito dibujado debe crear sustituyente enlazado."""
    canvas = ChemusonCanvas()
    canvas.state.show_implicit_hydrogens = True

    c1 = canvas.model.add_atom("C", 180.0, 200.0)
    c2 = canvas.model.add_atom("C", 220.0, 200.0)
    bond = canvas.model.add_bond(
        c1.id,
        c2.id,
        order=1,
        style=BondStyle.PLAIN,
        stereo=BondStereo.NONE,
        is_aromatic=False,
    )

    canvas.add_atom_item(c1)
    canvas.add_atom_item(c2)
    canvas.add_bond_item(bond)
    canvas._refresh_implicit_h_overlays([c1.id])

    overlays = canvas._implicit_h_overlays.get(c1.id)
    assert overlays, "Expected implicit H overlays around anchor carbon"

    _line_item, text_item = overlays[0]
    scene_click = text_item.mapToScene(text_item.boundingRect().center())

    anchor_id, angle_value = canvas._pick_implicit_h_overlay(scene_click)
    assert anchor_id == c1.id
    assert canvas._substitute_implicit_hydrogen(anchor_id, "Cl", angle_value)

    chlorine_ids = [atom.id for atom in canvas.model.atoms.values() if atom.element == "Cl"]
    assert len(chlorine_ids) == 1
    cl_id = chlorine_ids[0]

    assert any(
        (bond_item.a1_id == c1.id and bond_item.a2_id == cl_id)
        or (bond_item.a1_id == cl_id and bond_item.a2_id == c1.id)
        for bond_item in canvas.model.bonds.values()
    )
    assert canvas.model.implicit_h_count(c1.id) == 2

