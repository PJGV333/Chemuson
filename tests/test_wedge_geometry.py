"""Pruebas unitarias para test_wedge_geometry."""

import os
import sys
import math

import pytest
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import Atom, Bond, BondStyle, BondStereo
from gui.items import (
    BondItem,
    WEDGE_WIDE_END_MITER_OVERLAP_PX,
    WEDGE_WIDE_END_MITER_OVERLAP_STROKE_MULT,
    WEDGE_WIDE_END_MITER_MAX_DIST_STROKE_MULT,
    WEDGE_WIDE_END_MITER_MAX_DIST_WIDTH_MULT,
    WEDGE_WIDE_END_MITER_BACKTRACK_STROKE_MULT,
)
from gui.wedge_geometry import compute_wedge_points


def _midpoint(p0, p1):
    """Función de prueba auxiliar para  midpoint.

    Args:
        p0: Descripción del parámetro.
        p1: Descripción del parámetro.

    Returns:
        None.

    """
    return ((p0[0] + p1[0]) / 2.0, (p0[1] + p1[1]) / 2.0)


def _normalize(vx, vy):
    length = math.hypot(vx, vy)
    return (vx / length, vy / length)


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def test_wedge_points_ring_no_trim():
    """Verifica wedge points ring no trim.

    Returns:
        None.

    """
    tip, base1, base2 = compute_wedge_points((0.0, 0.0), (10.0, 0.0), 4.0)
    assert tip == (0.0, 0.0)
    mid = _midpoint(base1, base2)
    assert mid[0] == pytest.approx(10.0, abs=1e-6)
    assert mid[1] == pytest.approx(0.0, abs=1e-6)


def test_wedge_points_non_ring_trim():
    """Verifica wedge points non ring trim.

    Returns:
        None.

    """
    tip, base1, base2 = compute_wedge_points(
        (0.0, 0.0),
        (10.0, 0.0),
        4.0,
        trim_start=2.0,
        trim_end=1.0,
    )
    assert tip[0] == pytest.approx(2.0, abs=1e-6)
    assert tip[1] == pytest.approx(0.0, abs=1e-6)
    mid = _midpoint(base1, base2)
    assert mid[0] == pytest.approx(9.0, abs=1e-6)
    assert mid[1] == pytest.approx(0.0, abs=1e-6)


def test_wide_end_miter_chair_case_tracks_neighbor_edges():
    """La base ancha se recorta sobre las edge lines vecinas sin picos."""
    atom_a = Atom(id=1, element="C", x=0.0, y=0.0)
    atom_b = Atom(id=2, element="C", x=40.0, y=0.0)
    stroke_px = 4.0
    wedge_width = 16.0
    bond = Bond(
        id=1,
        a1_id=1,
        a2_id=2,
        style=BondStyle.WEDGE,
        stereo=BondStereo.UP,
        stroke_px=stroke_px,
    )
    item = BondItem(bond, atom_a, atom_b)

    n1x, n1y = _normalize(0.86, 0.50)
    n2x, n2y = _normalize(0.10, -0.99)
    base_cx, base_cy = 40.0, 0.0
    neighbors = [
        (n1x, n1y, stroke_px, base_cx, base_cy),
        (n2x, n2y, stroke_px, base_cx, base_cy),
    ]
    item.set_wedge_join_neighbors([], neighbors)

    tip_pos = (0.0, 0.0)
    tip_neg = (0.0, 0.0)
    base_pos = (40.0, wedge_width * 0.5)
    base_neg = (40.0, -wedge_width * 0.5)
    axis_ux, axis_uy = 1.0, 0.0

    new_pos, new_neg = item._miter_wedge_wide_end_into_neighbors(
        tip_pos,
        tip_neg,
        base_pos,
        base_neg,
        base_cx,
        base_cy,
        axis_ux,
        axis_uy,
        stroke_px,
        wedge_width,
    )

    def edge_line_for_corner(corner, neighbor):
        nux, nuy, nwidth, edge_cx, edge_cy = neighbor
        bnx, bny = -nuy, nux
        half_neighbor = max(nwidth * 0.5, stroke_px * 0.5)
        side = 1.0 if ((corner[0] - base_cx) * bnx + (corner[1] - base_cy) * bny) >= 0.0 else -1.0
        edge_px = edge_cx + bnx * half_neighbor * side
        edge_py = edge_cy + bny * half_neighbor * side
        return (edge_px, edge_py), (nux, nuy), (bnx, bny), half_neighbor

    pos_neighbor = item._pick_wedge_neighbor_for_corner(
        base_pos, base_cx, base_cy, axis_ux, axis_uy, neighbors, stroke_px
    )
    neg_neighbor = item._pick_wedge_neighbor_for_corner(
        base_neg, base_cx, base_cy, axis_ux, axis_uy, neighbors, stroke_px
    )
    assert pos_neighbor is not None
    assert neg_neighbor is not None

    for original_corner, new_corner, neighbor in (
        (base_pos, new_pos, pos_neighbor),
        (base_neg, new_neg, neg_neighbor),
    ):
        (edge_px, edge_py), (nux, nuy), (bnx, bny), half_neighbor = edge_line_for_corner(
            original_corner, neighbor
        )
        perp_dist = abs((new_corner[0] - edge_px) * bnx + (new_corner[1] - edge_py) * bny)
        assert perp_dist < 0.25

        overlap_raw = max(
            WEDGE_WIDE_END_MITER_OVERLAP_PX,
            stroke_px * WEDGE_WIDE_END_MITER_OVERLAP_STROKE_MULT,
        )
        if item._style.cap_style in (Qt.PenCapStyle.RoundCap, Qt.PenCapStyle.SquareCap):
            overlap_clamp = max(overlap_raw, half_neighbor)
        else:
            overlap_clamp = min(overlap_raw, 0.25 * half_neighbor, 0.8 * stroke_px)
        max_dist = max(
            stroke_px * WEDGE_WIDE_END_MITER_MAX_DIST_STROKE_MULT,
            wedge_width * WEDGE_WIDE_END_MITER_MAX_DIST_WIDTH_MULT,
        )
        overlap_clamp = min(overlap_clamp, max_dist * 0.35)
        forward = (new_corner[0] - base_cx) * nux + (new_corner[1] - base_cy) * nuy
        min_back = -max(stroke_px * WEDGE_WIDE_END_MITER_BACKTRACK_STROKE_MULT, 0.2)
        assert forward >= min_back - 1e-6
        assert forward <= max_dist + overlap_clamp + 1e-6


def test_wide_end_miter_uses_neighbor_edge_center_roundcap():
    """Con RoundCap, el miter usa edge_center del vecino, no base_center."""
    atom_a = Atom(id=1, element="C", x=0.0, y=0.0)
    atom_b = Atom(id=2, element="C", x=40.0, y=0.0)
    stroke_px = 4.0
    wedge_width = 16.0
    bond = Bond(
        id=2,
        a1_id=1,
        a2_id=2,
        style=BondStyle.WEDGE,
        stereo=BondStereo.UP,
        stroke_px=stroke_px,
    )
    item = BondItem(bond, atom_a, atom_b)

    nux, nuy = _normalize(0.86, 0.50)
    neighbor = (nux, nuy, stroke_px, 36.5, 0.0)  # shifted edge center
    tip = (0.0, 0.0)
    base_corner = (40.0, 8.0)
    base_center = (40.0, 0.0)
    result = item._miter_wedge_corner_into_neighbor(
        tip,
        base_corner,
        base_center,
        neighbor,
        stroke_px,
        wedge_width,
    )

    bnx, bny = -nuy, nux
    half_neighbor = max(stroke_px * 0.5, stroke_px * 0.5)
    side = 1.0 if ((base_corner[0] - base_center[0]) * bnx + (base_corner[1] - base_center[1]) * bny) >= 0.0 else -1.0
    edge_px = neighbor[3] + bnx * half_neighbor * side
    edge_py = neighbor[4] + bny * half_neighbor * side
    dist_to_edge_center_line = abs((result[0] - edge_px) * bnx + (result[1] - edge_py) * bny)

    old_edge_px = base_center[0] + bnx * half_neighbor * side
    old_edge_py = base_center[1] + bny * half_neighbor * side
    dist_to_base_center_line = abs((result[0] - old_edge_px) * bnx + (result[1] - old_edge_py) * bny)

    assert dist_to_edge_center_line < 0.25
    assert dist_to_base_center_line > 0.75


def test_bond_z_order_stereo_priority():
    """Wedge encima de plain; hashed debajo."""
    atom_a = Atom(id=10, element="C", x=0.0, y=0.0)
    atom_b = Atom(id=11, element="C", x=20.0, y=0.0)
    plain = BondItem(
        Bond(id=10, a1_id=10, a2_id=11, style=BondStyle.PLAIN, stroke_px=2.0),
        atom_a,
        atom_b,
    )
    wedge = BondItem(
        Bond(id=11, a1_id=10, a2_id=11, style=BondStyle.WEDGE, stereo=BondStereo.UP, stroke_px=2.0),
        atom_a,
        atom_b,
    )
    hashed = BondItem(
        Bond(id=12, a1_id=10, a2_id=11, style=BondStyle.HASHED, stereo=BondStereo.DOWN, stroke_px=2.0),
        atom_a,
        atom_b,
    )
    assert wedge.zValue() > plain.zValue()
    assert plain.zValue() > hashed.zValue()
