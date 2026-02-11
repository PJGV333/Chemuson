"""Regresiones para la herramienta de cadena zig-zag."""

import os
import sys

import pytest
from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import BondStyle, BondStereo
from gui.canvas import ChemusonCanvas
from gui.geom import angle_deg, angle_distance_deg


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def _expected_first_chain_angle(canvas: ChemusonCanvas, anchor_id: int, scene_pos: QPointF) -> float:
    anchor = canvas.model.get_atom(anchor_id)
    p0 = QPointF(anchor.x, anchor.y)
    raw_angle = angle_deg(p0, scene_pos)
    first_angle, _ = canvas._pick_bond_direction_deg(
        p0,
        anchor_id,
        raw_angle,
        1,
        False,
        canvas.state.bond_length,
        apply_collisions=True,
        allow_length_boost=canvas.state.fixed_lengths,
    )
    return first_angle


def test_chain_from_aromatic_anchor_keeps_zigzag_pattern():
    """Una cadena desde centro aromático no debe degenerar a línea recta."""
    canvas = ChemusonCanvas()
    canvas.state.fixed_angles = True
    canvas.state.fixed_lengths = True

    anchor = canvas.model.add_atom("C", 200.0, 200.0)
    n1 = canvas.model.add_atom("C", 170.0, 182.679492)
    n2 = canvas.model.add_atom("C", 170.0, 217.320508)
    canvas.model.add_bond(
        anchor.id,
        n1.id,
        order=1,
        style=BondStyle.PLAIN,
        stereo=BondStereo.NONE,
        is_aromatic=True,
    )
    canvas.model.add_bond(
        anchor.id,
        n2.id,
        order=1,
        style=BondStyle.PLAIN,
        stereo=BondStereo.NONE,
        is_aromatic=True,
    )

    scene_pos = QPointF(360.0, 200.0)
    canvas._begin_place_chain(anchor.id, scene_pos)
    points = canvas._chain_last_points
    assert points is not None and len(points) >= 3

    first = angle_deg(points[0], points[1])
    second = angle_deg(points[1], points[2])
    assert angle_distance_deg(first, second) > 1.0


def test_chain_first_bond_matches_normal_single_bond_orientation():
    """El primer enlace de la cadena respeta el ángulo de enlace sencillo normal."""
    canvas = ChemusonCanvas()
    canvas.state.fixed_angles = True
    canvas.state.fixed_lengths = True

    anchor = canvas.model.add_atom("C", 200.0, 200.0)
    c1 = canvas.model.add_atom("C", 160.0, 200.0)
    c2 = canvas.model.add_atom("C", 220.0, 234.641016)
    canvas.model.add_bond(
        anchor.id,
        c1.id,
        order=1,
        style=BondStyle.PLAIN,
        stereo=BondStereo.NONE,
        is_aromatic=False,
    )
    canvas.model.add_bond(
        anchor.id,
        c2.id,
        order=1,
        style=BondStyle.PLAIN,
        stereo=BondStereo.NONE,
        is_aromatic=False,
    )

    scene_pos = QPointF(280.0, 140.0)
    expected_first = _expected_first_chain_angle(canvas, anchor.id, scene_pos)
    canvas._begin_place_chain(anchor.id, scene_pos)
    points = canvas._chain_last_points
    assert points is not None and len(points) >= 3

    first = angle_deg(points[0], points[1])
    second = angle_deg(points[1], points[2])
    assert angle_distance_deg(first, expected_first) < 1e-6
    assert angle_distance_deg(first, second) > 1.0
