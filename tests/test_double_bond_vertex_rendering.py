"""Regresiones para la geometría de vértices en dobles enlaces."""

import os
import sys

import pytest
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import BondStyle, BondStereo
from gui.canvas import ChemusonCanvas


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def _bond_path_elements(item):
    path = item.path()
    assert path.elementCount() >= 4
    return [path.elementAt(i) for i in range(4)]


def test_cc_double_keeps_trimmed_pi_line_but_cn_is_full():
    """C=C conserva recorte interno; C=N mantiene doble completo."""
    canvas = ChemusonCanvas()

    c1 = canvas.model.add_atom("C", 100.0, 100.0)
    c2 = canvas.model.add_atom("C", 140.0, 100.0)
    b_cc = canvas.model.add_bond(
        c1.id,
        c2.id,
        order=2,
        style=BondStyle.PLAIN,
        stereo=BondStereo.NONE,
        is_aromatic=False,
    )

    c3 = canvas.model.add_atom("C", 100.0, 200.0)
    n1 = canvas.model.add_atom("N", 140.0, 200.0)
    b_cn = canvas.model.add_bond(
        c3.id,
        n1.id,
        order=2,
        style=BondStyle.PLAIN,
        stereo=BondStereo.NONE,
        is_aromatic=False,
    )

    canvas._rebuild_items_from_model()

    cc0, cc1, cc2, cc3 = _bond_path_elements(canvas.bond_items[b_cc.id])
    assert abs(cc0.y - 100.0) < 1e-6
    assert abs(cc1.y - 100.0) < 1e-6
    assert cc2.x > cc0.x + 0.5
    assert cc3.x < cc1.x - 0.5

    cn0, cn1, cn2, cn3 = _bond_path_elements(canvas.bond_items[b_cn.id])
    assert abs(cn0.y - 200.0) < 1e-6
    assert abs(cn1.y - 200.0) < 1e-6
    assert abs(cn2.x - cn0.x) < 1e-6
    assert abs(cn3.x - cn1.x) < 1e-6


def test_double_bond_secondary_line_is_always_on_screen_left():
    """La línea pi del doble enlace debe quedar siempre a la izquierda en pantalla."""
    canvas = ChemusonCanvas()

    # Pendiente ascendente hacia la derecha.
    a1 = canvas.model.add_atom("C", 100.0, 260.0)
    a2 = canvas.model.add_atom("C", 160.0, 180.0)
    b1 = canvas.model.add_bond(
        a1.id,
        a2.id,
        order=2,
        style=BondStyle.PLAIN,
        stereo=BondStereo.NONE,
        is_aromatic=False,
    )

    # Pendiente descendente hacia la derecha (orientación opuesta).
    a3 = canvas.model.add_atom("C", 260.0, 120.0)
    a4 = canvas.model.add_atom("C", 200.0, 200.0)
    b2 = canvas.model.add_bond(
        a3.id,
        a4.id,
        order=2,
        style=BondStyle.PLAIN,
        stereo=BondStereo.NONE,
        is_aromatic=False,
    )

    canvas._rebuild_items_from_model()

    for bond_id in (b1.id, b2.id):
        p0, p1, p2, p3 = _bond_path_elements(canvas.bond_items[bond_id])
        sigma_mid_x = 0.5 * (p0.x + p1.x)
        pi_mid_x = 0.5 * (p2.x + p3.x)
        assert pi_mid_x < sigma_mid_x - 1e-6
