"""Regresiones de visualización para círculos aromáticos."""

import math
import os
import sys

import pytest
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.core.model import BondStyle, BondStereo
from chemuson.gui.canvas import AROMATIC_CIRCLE_ATOMS_ROLE, ChemusonCanvas


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def _build_hexagon(canvas: ChemusonCanvas, center_x: float = 220.0, center_y: float = 220.0):
    ids = []
    radius = 42.0
    for i in range(6):
        theta = math.radians(60.0 * i)
        atom = canvas.model.add_atom(
            "C",
            center_x + radius * math.cos(theta),
            center_y + radius * math.sin(theta),
        )
        ids.append(atom.id)
    return ids


def _add_aromatic_ring(
    canvas: ChemusonCanvas,
    center_x: float = 220.0,
    center_y: float = 220.0,
) -> list[int]:
    """Construye un anillo aromático simple y devuelve sus IDs de átomos."""
    atom_ids = _build_hexagon(canvas, center_x=center_x, center_y=center_y)
    for i in range(6):
        canvas.model.add_bond(
            atom_ids[i],
            atom_ids[(i + 1) % 6],
            order=1,
            style=BondStyle.PLAIN,
            stereo=BondStereo.NONE,
            is_aromatic=True,
        )
    canvas._kekulize_aromatic_bonds()
    return atom_ids


def _angle_diff_mod_180(a_deg: float, b_deg: float) -> float:
    """Diferencia angular absoluta considerando periodicidad de 180°."""
    delta = (float(a_deg) - float(b_deg) + 90.0) % 180.0 - 90.0
    return abs(delta)


def test_aromatic_circle_survives_single_wedge_styling():
    """Un anillo aromático con una cuña debe conservar el círculo aromático."""
    canvas = ChemusonCanvas()
    atom_ids = _build_hexagon(canvas)

    for i in range(6):
        canvas.model.add_bond(
            atom_ids[i],
            atom_ids[(i + 1) % 6],
            order=1,
            style=BondStyle.PLAIN,
            stereo=BondStereo.NONE,
            is_aromatic=True,
        )
    canvas._kekulize_aromatic_bonds()

    # Simula un caso real donde el enlace estilizado queda no-aromático.
    wedge_bond = canvas.model.get_bond(1)
    canvas.model.update_bond(
        wedge_bond.id,
        style=BondStyle.WEDGE,
        stereo=BondStereo.UP,
        is_aromatic=False,
    )

    canvas._rebuild_items_from_model()
    canvas.state.use_aromatic_circles = True
    canvas.refresh_aromatic_circles()

    assert len(canvas.aromatic_circles) == 1


def test_aromatic_circle_survives_two_stereo_styled_ring_bonds():
    """El círculo aromático se mantiene con dos enlaces de anillo en cuña/hash."""
    canvas = ChemusonCanvas()
    atom_ids = _build_hexagon(canvas, center_x=340.0, center_y=220.0)

    for i in range(6):
        canvas.model.add_bond(
            atom_ids[i],
            atom_ids[(i + 1) % 6],
            order=1,
            style=BondStyle.PLAIN,
            stereo=BondStereo.NONE,
            is_aromatic=True,
        )
    canvas._kekulize_aromatic_bonds()

    bond_1 = canvas.model.get_bond(1)
    bond_2 = canvas.model.get_bond(2)
    canvas.model.update_bond(
        bond_1.id,
        style=BondStyle.WEDGE,
        stereo=BondStereo.UP,
        is_aromatic=False,
    )
    canvas.model.update_bond(
        bond_2.id,
        style=BondStyle.HASHED,
        stereo=BondStereo.DOWN,
        is_aromatic=False,
    )

    canvas._rebuild_items_from_model()
    canvas.state.use_aromatic_circles = True
    canvas.refresh_aromatic_circles()

    assert len(canvas.aromatic_circles) == 1


def test_aromatic_circle_survives_single_bold_styled_ring_bond():
    """El círculo aromático se mantiene con un enlace de anillo en estilo bold."""
    canvas = ChemusonCanvas()
    atom_ids = _build_hexagon(canvas, center_x=280.0, center_y=320.0)

    for i in range(6):
        canvas.model.add_bond(
            atom_ids[i],
            atom_ids[(i + 1) % 6],
            order=1,
            style=BondStyle.PLAIN,
            stereo=BondStereo.NONE,
            is_aromatic=True,
        )
    canvas._kekulize_aromatic_bonds()

    bold_bond = canvas.model.get_bond(1)
    canvas.model.update_bond(
        bold_bond.id,
        style=BondStyle.BOLD,
        stereo=BondStereo.NONE,
        is_aromatic=False,
    )

    canvas._rebuild_items_from_model()
    canvas.state.use_aromatic_circles = True
    canvas.refresh_aromatic_circles()

    assert len(canvas.aromatic_circles) == 1


def test_apply_bold_style_on_aromatic_bond_preserves_aromaticity():
    """Aplicar bold desde la herramienta no debe romper aromaticidad del enlace."""
    canvas = ChemusonCanvas()
    atom_ids = _add_aromatic_ring(canvas, center_x=380.0, center_y=320.0)
    canvas._rebuild_items_from_model()

    target_bond_id = next(iter(canvas.model.bonds.keys()))
    canvas.state.active_bond_order = 1
    canvas.state.active_bond_style = BondStyle.BOLD
    canvas.state.active_bond_stereo = BondStereo.NONE
    canvas.state.active_bond_mode = "set"
    canvas.state.active_bond_aromatic = False
    canvas._apply_bond_style(target_bond_id)

    target_bond = canvas.model.get_bond(target_bond_id)
    assert target_bond.is_aromatic
    assert target_bond.style == BondStyle.BOLD

    canvas.state.use_aromatic_circles = True
    canvas.refresh_aromatic_circles()
    assert len(canvas.aromatic_circles) == 1


def test_aromatic_circle_rotation_clears_stale_trackball_reference():
    """Al rotar en 2D, los círculos deben seguir al anillo aunque exista referencia trackball."""
    canvas = ChemusonCanvas()
    ring_a = _add_aromatic_ring(canvas, center_x=180.0, center_y=220.0)
    ring_b = _add_aromatic_ring(canvas, center_x=300.0, center_y=220.0)
    all_atom_ids = set(ring_a + ring_b)

    canvas._rebuild_items_from_model()
    canvas.state.use_aromatic_circles = True
    canvas.refresh_aromatic_circles()
    assert len(canvas.aromatic_circles) == 2

    # Simula estado residual de trackball previo.
    canvas._rotation_3d_ref_atom_ids = tuple(sorted(all_atom_ids))
    canvas._rotation_3d_ref_positions = {
        atom_id: (canvas.model.get_atom(atom_id).x, canvas.model.get_atom(atom_id).y)
        for atom_id in all_atom_ids
    }
    canvas._rotation_3d_pitch_deg = 0.0
    canvas._rotation_3d_yaw_deg = 0.0

    canvas.state.selected_atoms = set(all_atom_ids)
    canvas.state.selected_bonds = set(canvas.model.bonds.keys())
    canvas.rotate_selection_degrees(0.5, use_start_positions=False)

    assert canvas._rotation_3d_ref_atom_ids == tuple()
    assert canvas._rotation_3d_ref_positions == {}

    for circle in canvas.aromatic_circles:
        ring_atom_ids = tuple(int(atom_id) for atom_id in circle.data(AROMATIC_CIRCLE_ATOMS_ROLE))
        cx = sum(canvas.model.get_atom(atom_id).x for atom_id in ring_atom_ids) / len(ring_atom_ids)
        cy = sum(canvas.model.get_atom(atom_id).y for atom_id in ring_atom_ids) / len(ring_atom_ids)
        center = circle.sceneBoundingRect().center()
        assert math.hypot(center.x() - cx, center.y() - cy) < 0.25


def test_aromatic_circle_keeps_ellipse_on_move_and_rotate():
    """Si el anillo está deformado (elipse), el círculo debe seguir move/rotate."""
    canvas = ChemusonCanvas()
    ring_atom_ids = _add_aromatic_ring(canvas, center_x=260.0, center_y=260.0)
    canvas._rebuild_items_from_model()
    canvas.state.use_aromatic_circles = True
    canvas.refresh_aromatic_circles()
    assert len(canvas.aromatic_circles) == 1

    cx = sum(canvas.model.get_atom(atom_id).x for atom_id in ring_atom_ids) / len(ring_atom_ids)
    cy = sum(canvas.model.get_atom(atom_id).y for atom_id in ring_atom_ids) / len(ring_atom_ids)
    theta = math.radians(32.0)
    cos_t = math.cos(theta)
    sin_t = math.sin(theta)
    sx = 1.35
    sy = 0.72
    m00 = sx * cos_t * cos_t + sy * sin_t * sin_t
    m01 = (sx - sy) * cos_t * sin_t
    m10 = m01
    m11 = sx * sin_t * sin_t + sy * cos_t * cos_t

    for atom_id in ring_atom_ids:
        atom = canvas.model.get_atom(atom_id)
        dx = atom.x - cx
        dy = atom.y - cy
        nx = cx + m00 * dx + m01 * dy
        ny = cy + m10 * dx + m11 * dy
        canvas.model.update_atom_position(atom_id, nx, ny)
        canvas.update_atom_item(atom_id, nx, ny)
    canvas.update_bond_items_for_atoms(set(ring_atom_ids))

    circle = canvas.aromatic_circles[0]
    pos0 = circle.pos()
    rx0 = circle.rect().width() * 0.5
    ry0 = circle.rect().height() * 0.5
    ang0 = float(circle.rotation())
    assert rx0 > ry0 + 0.25

    move_dx = 34.0
    move_dy = -21.0
    for atom_id in ring_atom_ids:
        atom = canvas.model.get_atom(atom_id)
        nx = atom.x + move_dx
        ny = atom.y + move_dy
        canvas.model.update_atom_position(atom_id, nx, ny)
        canvas.update_atom_item(atom_id, nx, ny)
    canvas.update_bond_items_for_atoms(set(ring_atom_ids))

    circle = canvas.aromatic_circles[0]
    pos1 = circle.pos()
    rx1 = circle.rect().width() * 0.5
    ry1 = circle.rect().height() * 0.5
    ang1 = float(circle.rotation())
    assert math.hypot(pos1.x() - (pos0.x() + move_dx), pos1.y() - (pos0.y() + move_dy)) < 0.25
    assert abs(rx1 - rx0) < 0.2
    assert abs(ry1 - ry0) < 0.2
    assert _angle_diff_mod_180(ang1, ang0) < 0.5

    rot_deg = 22.0
    rot_rad = math.radians(rot_deg)
    cos_r = math.cos(rot_rad)
    sin_r = math.sin(rot_rad)
    rot_cx = pos1.x()
    rot_cy = pos1.y()
    for atom_id in ring_atom_ids:
        atom = canvas.model.get_atom(atom_id)
        dx = atom.x - rot_cx
        dy = atom.y - rot_cy
        nx = rot_cx + dx * cos_r - dy * sin_r
        ny = rot_cy + dx * sin_r + dy * cos_r
        canvas.model.update_atom_position(atom_id, nx, ny)
        canvas.update_atom_item(atom_id, nx, ny)
    canvas.update_bond_items_for_atoms(set(ring_atom_ids))

    circle = canvas.aromatic_circles[0]
    pos2 = circle.pos()
    rx2 = circle.rect().width() * 0.5
    ry2 = circle.rect().height() * 0.5
    ang2 = float(circle.rotation())
    assert math.hypot(pos2.x() - rot_cx, pos2.y() - rot_cy) < 0.25
    assert abs(rx2 - rx0) < 0.3
    assert abs(ry2 - ry0) < 0.3
    assert abs(_angle_diff_mod_180(ang2, ang0) - rot_deg) < 1.5


@pytest.mark.parametrize("edit_mode", ["double", "increment"])
def test_partial_aromatic_ring_does_not_hide_pi_lines_with_circle_mode(edit_mode: str):
    """Con círculos activos, un anillo parcialmente aromático no debe ocultar líneas pi."""
    canvas = ChemusonCanvas()
    _add_aromatic_ring(canvas, center_x=260.0, center_y=220.0)
    canvas._rebuild_items_from_model()
    canvas.state.use_aromatic_circles = True
    canvas.refresh_aromatic_circles()
    assert len(canvas.aromatic_circles) == 1

    target_bond_id = next(iter(canvas.model.bonds.keys()))
    if edit_mode == "double":
        canvas.state.active_bond_order = 2
        canvas.state.active_bond_style = BondStyle.PLAIN
        canvas.state.active_bond_stereo = BondStereo.NONE
        canvas.state.active_bond_mode = "set"
        canvas.state.active_bond_aromatic = False
        canvas._apply_bond_style(target_bond_id)
    else:
        canvas.state.active_bond_mode = "increment"
        canvas.state.active_bond_aromatic = False
        canvas._cycle_bond_order(target_bond_id)

    assert len(canvas.aromatic_circles) == 0

    remaining_aromatic = [bond for bond in canvas.model.bonds.values() if bond.is_aromatic]
    assert remaining_aromatic
    probe = remaining_aromatic[0]
    canvas.model.update_bond(
        probe.id,
        order=2,
        style=BondStyle.PLAIN,
        stereo=BondStereo.NONE,
        is_aromatic=True,
        display_order=2,
    )
    canvas.update_bond_item(probe.id)
    canvas.refresh_aromatic_circles()

    assert len(canvas.aromatic_circles) == 0
    assert canvas.bond_items[probe.id].path().elementCount() > 2
