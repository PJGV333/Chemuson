"""Regresiones para colocación de enlaces sp3 en carbonos tetravalentes."""

import math
import os
import sys

import pytest
from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.core.model import BondStyle, BondStereo
from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.geom import angle_deg, angle_distance_deg, endpoint_from_angle_len


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
    # Debe salir por una bisectriz abierta (60°, 180° o 300° para este arreglo).
    assert min(angle_distance_deg(new_angle, expected) for expected in (60.0, 180.0, 300.0)) < 1e-6
    # Evita geometrías casi paralelas al enlace existente (regresión: 10.5°).
    assert min(angle_distance_deg(new_angle, existing) for existing in existing_angles) >= 55.0

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


def test_sp3_fourth_bond_prioritizes_open_gap_for_skewed_tripod():
    """Con 3 enlaces tipo trípode sesgado, el 4º debe ir al hueco amplio."""
    canvas = ChemusonCanvas()

    center = canvas.model.add_atom("C", 260.0, 260.0)
    p0 = QPointF(center.x, center.y)
    length = canvas.state.bond_length

    # Caso visual cercano al reportado: 180°, 60°, 300°.
    for theta in (180.0, 60.0, 300.0):
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
    new_angle = angle_deg(p0, p1)
    existing_angles = canvas._get_anchor_bond_angles_deg(center.id)

    # Huecos esperados para ese trípode: 0°, 120°, 240°.
    assert min(angle_distance_deg(new_angle, expected) for expected in (0.0, 120.0, 240.0)) < 1e-6
    assert min(angle_distance_deg(new_angle, existing) for existing in existing_angles) >= 55.0


def test_sp3_congested_fourth_bond_follows_cursor_across_gap_midpoints():
    """Con 3 enlaces sp3 desbalanceados, el 4º debe poder elegir cualquier bisectriz."""
    canvas = ChemusonCanvas()

    center = canvas.model.add_atom("C", 240.0, 240.0)
    p0 = QPointF(center.x, center.y)
    length = canvas.state.bond_length

    # Arreglo no simétrico: huecos de 110°, 140°, 110°.
    for theta in (40.0, 150.0, 290.0):
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
    theta, _ = canvas._pick_bond_direction_deg(
        p0,
        center.id,
        95.0,  # bisectriz del hueco 40°-150°
        1,
        False,
        length,
        apply_collisions=False,
        allow_length_boost=False,
    )

    # Debe seguir la bisectriz cercana al cursor, no forzar solo el hueco máximo.
    assert angle_distance_deg(theta, 95.0) < 1e-6


def test_default_bond_endpoint_ignores_isolated_hidden_carbon_placeholder():
    """Un carbono huérfano invisible no debe alterar el snap ni la geometría del nuevo enlace."""
    canvas = ChemusonCanvas()

    center = canvas.model.add_atom("C", 200.0, 200.0)
    p0 = QPointF(center.x, center.y)
    length = canvas.state.bond_length

    neighbor_pos = endpoint_from_angle_len(p0, 180.0, length)
    neighbor = canvas.model.add_atom("C", neighbor_pos.x(), neighbor_pos.y())
    canvas.model.add_bond(center.id, neighbor.id)
    canvas._rebuild_items_from_model()

    baseline = canvas._compute_default_bond_endpoint(p0, center.id)

    canvas.model.add_atom("C", baseline.x() + 6.0, baseline.y())
    perturbed = canvas._compute_default_bond_endpoint(p0, center.id)

    assert perturbed.x() == pytest.approx(baseline.x())
    assert perturbed.y() == pytest.approx(baseline.y())


def test_deleting_bond_prunes_hidden_carbon_orphans_but_keeps_meaningful_atoms():
    """Al borrar un enlace, los carbonos invisibles aislados deben limpiarse automáticamente."""
    canvas = ChemusonCanvas()

    anchor = canvas.model.add_atom("N", 180.0, 180.0, is_explicit=True)
    orphan = canvas.model.add_atom("C", 220.0, 180.0)
    bond = canvas.model.add_bond(anchor.id, orphan.id)
    canvas._rebuild_items_from_model()

    canvas._delete_selection(set(), {bond.id})

    assert bond.id not in canvas.model.bonds
    assert anchor.id in canvas.model.atoms
    assert orphan.id not in canvas.model.atoms

    canvas.undo_stack.undo()

    assert bond.id in canvas.model.bonds
    assert anchor.id in canvas.model.atoms
    assert orphan.id in canvas.model.atoms
