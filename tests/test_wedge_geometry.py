"""Pruebas unitarias para test_wedge_geometry."""

import os
import sys
import math

import pytest
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.core.model import Atom, Bond, BondStyle, BondStereo
from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.items import (
    BondItem,
    WEDGE_WIDE_END_MITER_OVERLAP_PX,
    WEDGE_WIDE_END_MITER_OVERLAP_STROKE_MULT,
    WEDGE_WIDE_END_MITER_MAX_DIST_STROKE_MULT,
    WEDGE_WIDE_END_MITER_MAX_DIST_WIDTH_MULT,
    WEDGE_WIDE_END_MITER_BACKTRACK_STROKE_MULT,
)
from chemuson.gui.wedge_geometry import compute_wedge_points


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


def _default_wedge_base_corners(item: BondItem, atom_a: Atom, atom_b: Atom):
    """Reconstruye la base ideal previa al miter final."""
    dx = atom_b.x - atom_a.x
    dy = atom_b.y - atom_a.y
    length = math.hypot(dx, dy)
    assert length > 1e-6

    stroke_px = item._stroke_px if item._stroke_px is not None else item._style.stroke_px
    stroke_scale = stroke_px / item._style.stroke_px if item._style.stroke_px > 1e-6 else 1.0
    width = item._style.wedge_width_px * (0.72 + 0.28 * math.sqrt(max(stroke_scale, 1e-6)))
    width = max(width, stroke_px * 2.3)
    width = min(width, length * 0.34)

    trim_start = item._label_shrink_start + item._endpoint_trim_start
    trim_end = item._label_shrink_end + item._endpoint_trim_end
    tip, base1, base2 = compute_wedge_points(
        (atom_a.x, atom_a.y),
        (atom_b.x, atom_b.y),
        width,
        trim_start=trim_start,
        trim_end=trim_end,
    )
    base_cx = (base1[0] + base2[0]) * 0.5
    base_cy = (base1[1] + base2[1]) * 0.5
    wedge_dx = base_cx - tip[0]
    wedge_dy = base_cy - tip[1]
    wedge_len = math.hypot(wedge_dx, wedge_dy)
    assert wedge_len > 1e-6
    w_ux = wedge_dx / wedge_len
    w_uy = wedge_dy / wedge_len
    w_nx = -w_uy
    w_ny = w_ux
    half_w = width * 0.5
    base_pos = (base_cx + w_nx * half_w, base_cy + w_ny * half_w)
    base_neg = (base_cx - w_nx * half_w, base_cy - w_ny * half_w)
    return base_pos, base_neg, (base_cx, base_cy), (w_ux, w_uy), stroke_px


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


def test_aromatic_circle_mode_keeps_wedge_geometry():
    """En modo círculo aromático, una cuña aromática no debe colapsar a línea."""
    atom_a = Atom(id=20, element="C", x=0.0, y=0.0)
    atom_b = Atom(id=21, element="C", x=40.0, y=0.0)
    bond = Bond(
        id=30,
        a1_id=20,
        a2_id=21,
        style=BondStyle.WEDGE,
        stereo=BondStereo.UP,
        is_aromatic=True,
        display_order=2,
        stroke_px=2.0,
    )
    item = BondItem(
        bond,
        atom_a,
        atom_b,
        render_aromatic_as_circle=True,
    )
    rect = item.path().boundingRect()
    assert rect.height() > 0.5


def test_wedge_join_context_includes_bold_neighbor():
    """La integración de cuña en extremo ancho debe considerar enlace bold."""
    canvas = ChemusonCanvas()
    atom_a = canvas.model.add_atom("C", 80.0, 220.0)
    atom_b = canvas.model.add_atom("C", 140.0, 220.0)
    atom_c = canvas.model.add_atom("C", 190.0, 250.0)

    wedge = canvas.model.add_bond(
        atom_a.id,
        atom_b.id,
        style=BondStyle.WEDGE,
        stereo=BondStereo.UP,
    )
    canvas.model.add_bond(
        atom_b.id,
        atom_c.id,
        style=BondStyle.BOLD,
        stereo=BondStereo.NONE,
    )

    canvas._rebuild_items_from_model()
    wedge_item = canvas.bond_items[wedge.id]
    assert wedge_item._wedge_join_end, "No se detectaron vecinos en extremo ancho de la cuña"
    # El vecino bold debe entrar con un ancho mayor que el trazo base.
    max_neighbor_width = max(neighbor[2] for neighbor in wedge_item._wedge_join_end)
    assert max_neighbor_width > canvas.drawing_style.stroke_px * 1.5


def test_ring_wedge_respects_label_clearance_on_labeled_atom():
    """Una cuña en anillo no debe meterse debajo de la etiqueta del heteroátomo."""
    canvas = ChemusonCanvas()

    atom_n = canvas.model.add_atom("N", 200.0, 200.0, formal_charge=1)
    atom_a = canvas.model.add_atom("C", 240.0, 170.0)
    atom_b = canvas.model.add_atom("C", 280.0, 210.0)
    atom_c = canvas.model.add_atom("C", 250.0, 250.0)
    atom_d = canvas.model.add_atom("C", 190.0, 245.0)

    canvas.model.add_bond(atom_n.id, atom_a.id, style=BondStyle.BOLD)
    canvas.model.add_bond(atom_a.id, atom_b.id)
    canvas.model.add_bond(atom_b.id, atom_c.id)
    canvas.model.add_bond(atom_c.id, atom_d.id)
    wedge = canvas.model.add_bond(
        atom_d.id,
        atom_n.id,
        style=BondStyle.WEDGE,
        stereo=BondStereo.UP,
    )

    canvas._rebuild_items_from_model()

    item = canvas.bond_items[wedge.id]
    atom_item = canvas.atom_items[atom_n.id]
    wedge_path = item.mapToScene(item.path())
    label_shape = atom_item.label.mapToScene(atom_item.label.shape())
    charge_shape = atom_item.charge_label.mapToScene(atom_item.charge_label.shape())

    assert item._bond_in_ring
    assert item._label_shrink_end > 0.0
    assert item._label_shrink_end < 19.0
    assert not wedge_path.intersects(label_shape)
    assert not wedge_path.intersects(charge_shape)


def test_wedge_neighbor_selection_keeps_opposite_sides_with_bold_asymmetry():
    """Cada esquina de la base debe elegir el vecino de su lado aunque uno sea bold."""
    canvas = ChemusonCanvas()

    atom_n = canvas.model.add_atom("N", 200.0, 200.0, formal_charge=1)
    atom_tip = canvas.model.add_atom("C", 190.0, 245.0)
    atom_left = canvas.model.add_atom("C", 160.0, 220.0)
    atom_right = canvas.model.add_atom("C", 250.0, 220.0)

    wedge = canvas.model.add_bond(
        atom_tip.id,
        atom_n.id,
        style=BondStyle.WEDGE,
        stereo=BondStereo.UP,
    )
    plain = canvas.model.add_bond(atom_n.id, atom_left.id, style=BondStyle.PLAIN)
    bold = canvas.model.add_bond(atom_n.id, atom_right.id, style=BondStyle.BOLD)

    canvas._rebuild_items_from_model()

    item = canvas.bond_items[wedge.id]
    base_pos, base_neg, (base_cx, base_cy), (axis_ux, axis_uy), stroke_px = _default_wedge_base_corners(
        item,
        canvas.model.get_atom(wedge.a1_id),
        canvas.model.get_atom(wedge.a2_id),
    )

    pos_neighbor = item._pick_wedge_neighbor_for_corner(
        base_pos,
        base_cx,
        base_cy,
        axis_ux,
        axis_uy,
        item._wedge_join_end,
        stroke_px,
    )
    neg_neighbor = item._pick_wedge_neighbor_for_corner(
        base_neg,
        base_cx,
        base_cy,
        axis_ux,
        axis_uy,
        item._wedge_join_end,
        stroke_px,
    )

    assert pos_neighbor is not None
    assert neg_neighbor is not None

    expected_widths = sorted(
        [
            canvas._bond_render_width(canvas.model.get_bond(plain.id)),
            canvas._bond_render_width(canvas.model.get_bond(bold.id)),
        ]
    )
    selected_widths = sorted([pos_neighbor[2], neg_neighbor[2]])
    assert selected_widths == pytest.approx(expected_widths, rel=1e-6)


def test_wedge_tip_on_labeled_nitrogen_is_not_overtrimmed():
    """La punta de la cuña en un heteroátomo etiquetado no debe separarse en exceso."""
    canvas = ChemusonCanvas()

    atom_left = canvas.model.add_atom("C", 430.0, 320.0)
    atom_n = canvas.model.add_atom("N", 490.0, 280.0, formal_charge=1)
    atom_wedge = canvas.model.add_atom("C", 515.0, 338.0)
    atom_plain = canvas.model.add_atom("C", 560.0, 300.0)
    atom_bold = canvas.model.add_atom("C", 545.0, 365.0)

    canvas.model.add_bond(atom_left.id, atom_n.id, style=BondStyle.PLAIN)
    wedge = canvas.model.add_bond(
        atom_n.id,
        atom_wedge.id,
        style=BondStyle.WEDGE,
        stereo=BondStereo.UP,
    )
    canvas.model.add_bond(atom_n.id, atom_plain.id, style=BondStyle.PLAIN)
    canvas.model.add_bond(atom_wedge.id, atom_bold.id, style=BondStyle.BOLD)
    canvas.model.add_bond(atom_plain.id, atom_bold.id, style=BondStyle.PLAIN)

    canvas._rebuild_items_from_model()

    item = canvas.bond_items[wedge.id]
    atom_item = canvas.atom_items[atom_n.id]
    wedge_path = item.mapToScene(item.path())
    label_shape = atom_item.label.mapToScene(atom_item.label.shape())
    charge_shape = atom_item.charge_label.mapToScene(atom_item.charge_label.shape())

    assert item._label_shrink_start < 12.0
    assert not wedge_path.intersects(label_shape)
    assert not wedge_path.intersects(charge_shape)


def test_aromatic_double_bold_uses_thin_secondary_pi_line():
    """En doble aromático bold, el trazo pi debe conservar grosor normal."""
    atom_a = Atom(id=30, element="C", x=0.0, y=0.0)
    atom_b = Atom(id=31, element="C", x=40.0, y=0.0)
    bond = Bond(
        id=40,
        a1_id=30,
        a2_id=31,
        order=2,
        style=BondStyle.BOLD,
        stereo=BondStereo.NONE,
        is_aromatic=True,
        display_order=2,
        stroke_px=2.0,
    )
    item = BondItem(
        bond,
        atom_a,
        atom_b,
        render_aromatic_as_circle=False,
    )
    assert item._has_secondary_path
    assert item._secondary_pen is not None
    assert item.pen().widthF() > item._secondary_pen.widthF()
    assert item._secondary_pen.widthF() == pytest.approx(2.0, abs=1e-6)
