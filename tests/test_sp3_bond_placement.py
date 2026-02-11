"""Regresiones para colocación de enlaces sp3 en carbonos tetravalentes."""

import math
import os
import sys

import pytest
from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import BondStyle, BondStereo
from gui.canvas import ChemusonCanvas
from gui.geom import angle_deg, angle_distance_deg, endpoint_from_angle_len


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def test_sp3_default_fourth_bond_does_not_overlap_existing_bonds():
    """Con 3 enlaces en C(sp3), el 4º enlace por defecto no debe solapar uno existente."""
    canvas = ChemusonCanvas()

    center = canvas.model.add_atom("C", 220.0, 220.0)
    p0 = QPointF(center.x, center.y)
    length = canvas.state.bond_length

    # Geometría trigonal previa (caso reportado): 0°, 120°, 240°.
    for theta in (0.0, 120.0, 240.0):
        p = endpoint_from_angle_len(p0, theta, length)
        neighbor = canvas.model.add_atom("C", p.x(), p.y())
        canvas.model.add_bond(
            center.id,
            neighbor.id,
            order=1,
            style=BondStyle.PLAIN,
            stereo=BondStereo.NONE,
            is_aromatic=False,
        )

    canvas._rebuild_items_from_model()
    canvas.state.fixed_angles = True
    canvas.state.angle_step_deg = 30
    canvas.state.active_bond_order = 1
    canvas.state.active_bond_style = BondStyle.PLAIN
    canvas.state.active_bond_stereo = BondStereo.NONE
    canvas.state.active_bond_aromatic = False

    p1 = canvas._compute_default_bond_endpoint(p0, center.id)

    # No debe coincidir con un enlace ya existente (antes ciclaban orden).
    assert canvas._find_overlapping_bond(p0, p1) is None

    new_angle = angle_deg(p0, p1)
    existing_angles = canvas._get_anchor_bond_angles_deg(center.id)
    assert min(angle_distance_deg(new_angle, existing) for existing in existing_angles) > 5.0

    # Debe mantenerse cercanía a la separación tetraédrica.
    assert math.isclose(canvas._sp3_display_angle_deg(), 109.5, abs_tol=1e-6)


def test_sp3_regular_growth_keeps_grid_alignment():
    """Con 1 enlace en C(sp3) y ángulos fijos, el siguiente enlace debe caer en retícula."""
    canvas = ChemusonCanvas()

    center = canvas.model.add_atom("C", 200.0, 200.0)
    p0 = QPointF(center.x, center.y)
    length = canvas.state.bond_length

    # Un único enlace existente hacia la izquierda (180°).
    p_left = endpoint_from_angle_len(p0, 180.0, length)
    neighbor = canvas.model.add_atom("C", p_left.x(), p_left.y())
    canvas.model.add_bond(
        center.id,
        neighbor.id,
        order=1,
        style=BondStyle.PLAIN,
        stereo=BondStereo.NONE,
        is_aromatic=False,
    )

    canvas._rebuild_items_from_model()
    canvas.state.fixed_angles = True
    canvas.state.angle_step_deg = 30
    canvas.state.active_bond_order = 1
    canvas.state.active_bond_style = BondStyle.PLAIN
    canvas.state.active_bond_stereo = BondStereo.NONE
    canvas.state.active_bond_aromatic = False

    p1 = canvas._compute_default_bond_endpoint(p0, center.id)
    theta = angle_deg(p0, p1)

    # En modo normal (no 4º enlace), debe mantenerse alineado a la retícula.
    snapped = round(theta / 30.0) * 30.0
    assert math.isclose(theta, snapped % 360.0, abs_tol=1e-6)
