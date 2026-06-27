"""Regresiones para la herramienta de cadena zig-zag."""


import pytest
from PyQt6.QtCore import QPointF, Qt
from PyQt6.QtWidgets import QApplication


from chemuson.core.model import BondStyle, BondStereo
from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.geom import angle_deg, angle_distance_deg


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


def test_free_chain_zigzag_is_symmetric_around_cursor_axis():
    """Cadena libre: zig-zag simétrico alrededor del eje pedido por el cursor."""
    canvas = ChemusonCanvas()
    canvas.state.fixed_angles = True
    canvas.state.fixed_lengths = True
    canvas.state.angle_step_deg = 30

    start = QPointF(200.0, 220.0)
    end = QPointF(360.0, 220.0)  # eje horizontal
    canvas._begin_place_chain(None, start)
    canvas._update_chain_preview(end, Qt.KeyboardModifier.NoModifier)
    points = canvas._chain_last_points
    assert points is not None and len(points) >= 3

    axis = angle_deg(start, end)
    first = angle_deg(points[0], points[1])
    second = angle_deg(points[1], points[2])

    d1 = ((first - axis + 180.0) % 360.0) - 180.0
    d2 = ((second - axis + 180.0) % 360.0) - 180.0

    assert d1 * d2 < 0.0
    assert abs(abs(d1) - abs(d2)) < 1e-6
