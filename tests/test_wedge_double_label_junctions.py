from __future__ import annotations

import math

import pytest
from PyQt6.QtGui import QColor, QFontMetricsF, QImage, QPainterPath, QPainterPathStroker
from PyQt6.QtWidgets import QApplication

from chemuson.core.model import BondStyle, BondStereo, MolGraph
from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.templates import build_haworth_template
from tools.visual.render_wedge_double_label_gallery import HEIGHT, WIDTH, compose_comparison


@pytest.mark.parametrize("angle_deg,pi_sign", ((0.0, 1), (5.0, -1), (35.0, 1), (90.0, -1)))
def test_hashed_wedge_terminal_mark_integrates_with_double_bond(angle_deg, pi_sign):
    app = QApplication.instance() or QApplication([])
    canvas = ChemusonCanvas()
    tip = canvas.model.add_atom("C", 120.0, 200.0)
    junction = canvas.model.add_atom("C", 180.0, 200.0)
    wedge = canvas.model.add_bond(
        tip.id,
        junction.id,
        style=BondStyle.HASHED,
        stereo=BondStereo.DOWN,
    )
    angle = math.radians(angle_deg)
    double_end = canvas.model.add_atom(
        "C",
        junction.x + 60.0 * math.cos(angle),
        junction.y + 60.0 * math.sin(angle),
    )
    double = canvas.model.add_bond(
        junction.id,
        double_end.id,
        order=2,
        pi_offset_sign=pi_sign,
    )
    atom_snapshot = tuple(
        (atom.id, atom.element, atom.x, atom.y)
        for atom in canvas.model.atoms.values()
    )
    bond_snapshot = tuple(
        (bond.id, bond.a1_id, bond.a2_id, bond.order, bond.style, bond.stereo, bond.pi_offset_sign)
        for bond in canvas.model.bonds.values()
    )

    canvas._rebuild_items_from_model()

    wedge_item = canvas.bond_items[wedge.id]
    double_item = canvas.bond_items[double.id]
    assert wedge_item._hashed_double_join_end
    assert wedge_item.path().elementCount() == 12
    last_start = wedge_item.path().elementAt(10)
    last_end = wedge_item.path().elementAt(11)
    last_x = (last_start.x + last_end.x) * 0.5
    last_y = (last_start.y + last_end.y) * 0.5
    last_mark_gap = junction.x - last_x
    assert last_mark_gap == pytest.approx(2.5, abs=0.2)

    double_shape = double_item.mapToScene(double_item.shape())
    if double_item._has_secondary_path and double_item._secondary_pen is not None:
        stroker = QPainterPathStroker()
        stroker.setWidth(double_item._secondary_pen.widthF())
        stroker.setCapStyle(double_item._secondary_pen.capStyle())
        stroker.setJoinStyle(double_item._secondary_pen.joinStyle())
        double_shape.addPath(
            double_item.mapToScene(stroker.createStroke(double_item._secondary_path))
        )
    assert not wedge_item.mapToScene(wedge_item.shape()).intersects(double_shape)
    assert atom_snapshot == tuple(
        (atom.id, atom.element, atom.x, atom.y)
        for atom in canvas.model.atoms.values()
    )
    assert bond_snapshot == tuple(
        (bond.id, bond.a1_id, bond.a2_id, bond.order, bond.style, bond.stereo, bond.pi_offset_sign)
        for bond in canvas.model.bonds.values()
    )
    canvas.close()
    app.processEvents()


@pytest.mark.parametrize(
    "element,double_angle,label_angle",
    (("O", 0.0, 135.0), ("N", 90.0, 45.0), ("Cl", 45.0, 180.0)),
)
def test_double_bond_neighbor_labels_center_the_visible_anchor(element, double_angle, label_angle):
    app = QApplication.instance() or QApplication([])
    canvas = ChemusonCanvas()
    junction = canvas.model.add_atom("C", 200.0, 200.0)
    double_radians = math.radians(double_angle)
    double_end = canvas.model.add_atom(
        "C",
        junction.x + 60.0 * math.cos(double_radians),
        junction.y + 60.0 * math.sin(double_radians),
    )
    canvas.model.add_bond(junction.id, double_end.id, order=2)
    label_radians = math.radians(label_angle)
    label_atom = canvas.model.add_atom(
        element,
        junction.x + 60.0 * math.cos(label_radians),
        junction.y + 60.0 * math.sin(label_radians),
        is_explicit=True,
    )
    canvas.model.add_bond(junction.id, label_atom.id)

    canvas._rebuild_items_from_model()

    atom_item = canvas.atom_items[label_atom.id]
    anchor = atom_item._label_anchor_override or atom_item._first_label_token(atom_item._display_label)
    assert anchor is not None
    metrics = QFontMetricsF(atom_item.label.font())
    glyph = QPainterPath()
    glyph.addText(0.0, metrics.ascent(), atom_item.label.font(), anchor)
    glyph_center = glyph.boundingRect().center()
    expected_x = label_atom.x + atom_item._label_offset.x()
    expected_y = label_atom.y + atom_item._label_offset.y()
    visible_anchor_x = atom_item.pos().x() + atom_item.label.pos().x() + atom_item._label_anchor_offset + glyph_center.x()
    visible_anchor_y = atom_item.pos().y() + atom_item.label.pos().y() + glyph_center.y()

    assert visible_anchor_x == pytest.approx(expected_x, abs=0.05)
    assert visible_anchor_y == pytest.approx(expected_y, abs=0.1)
    canvas.close()
    app.processEvents()


def _haworth_graph_with_bold_wedge(ring_size: int) -> MolGraph:
    if ring_size == 6:
        graph = build_haworth_template(70.0, anomeric_up=True, bold_front=True)
        ring_atom = graph.get_atom(4)  # C3, incident to two bold front edges.
        for bond in graph.bonds.values():
            atom_a = graph.get_atom(bond.a1_id)
            atom_b = graph.get_atom(bond.a2_id)
            if (
                atom_a.id == ring_atom.id
                and atom_b.element == "O"
                and atom_b.group_h_cap == 1
            ):
                bond.a1_id, bond.a2_id = atom_b.id, atom_a.id
                bond.style = BondStyle.WEDGE
                bond.stereo = BondStereo.UP
                break
        else:  # pragma: no cover - template contract
            raise AssertionError("Haworth hexagon did not expose the expected C3-OH bond")
        return graph

    if ring_size != 5:
        raise ValueError(ring_size)
    graph = MolGraph()
    scale = 78.0
    coordinates = (
        (0.0, -0.72 * scale),
        (0.72 * scale, -0.20 * scale),
        (0.45 * scale, 0.62 * scale),
        (-0.45 * scale, 0.62 * scale),
        (-0.72 * scale, -0.20 * scale),
    )
    elements = ("O", "C", "C", "C", "C")
    ring = [
        graph.add_atom(element, x, y, is_explicit=element != "C")
        for element, (x, y) in zip(elements, coordinates)
    ]
    for index in range(len(ring)):
        graph.add_bond(
            ring[index].id,
            ring[(index + 1) % len(ring)].id,
            style=BondStyle.BOLD if index in {2, 3} else BondStyle.PLAIN,
        )
    anchor = ring[3]
    oh = graph.add_atom("OH", anchor.x - 0.72 * scale, anchor.y - 0.62 * scale, is_explicit=True)
    graph.add_bond(oh.id, anchor.id, style=BondStyle.WEDGE, stereo=BondStereo.UP)
    return graph


def _load_graph(canvas: ChemusonCanvas, graph: MolGraph) -> None:
    for atom in graph.atoms.values():
        canvas.model.add_atom(
            atom.element,
            atom.x,
            atom.y,
            atom_id=atom.id,
            charge=atom.charge,
            explicit_h=atom.explicit_h,
            is_explicit=atom.is_explicit,
            label_scale=atom.label_scale,
        )
    for bond in graph.bonds.values():
        canvas.model.add_bond(
            bond.a1_id,
            bond.a2_id,
            order=bond.order,
            bond_id=bond.id,
            style=bond.style,
            stereo=bond.stereo,
            pi_offset_sign=bond.pi_offset_sign,
        )


@pytest.mark.parametrize("ring_size", (5, 6))
def test_haworth_solid_wedge_joins_bold_ring_edges_without_spikes(ring_size):
    app = QApplication.instance() or QApplication([])
    canvas = ChemusonCanvas()
    graph = _haworth_graph_with_bold_wedge(ring_size)
    atom_snapshot = tuple(
        (atom.id, atom.element, atom.x, atom.y)
        for atom in graph.atoms.values()
    )
    bond_snapshot = tuple(
        (bond.id, bond.a1_id, bond.a2_id, bond.order, bond.style, bond.stereo)
        for bond in graph.bonds.values()
    )
    _load_graph(canvas, graph)

    canvas._rebuild_items_from_model()

    wedge_bond = next(bond for bond in canvas.model.bonds.values() if bond.style == BondStyle.WEDGE)
    wedge_item = canvas.bond_items[wedge_bond.id]
    anchor = canvas.model.get_atom(wedge_bond.a2_id)
    assert wedge_item._wedge_join_end_is_cyclic_bold
    assert len(wedge_item._wedge_join_end) == 2
    polygon = wedge_item.path().toFillPolygons()[0]
    base_distances = [
        math.hypot(point.x() - anchor.x, point.y() - anchor.y)
        for point in polygon
        if math.hypot(point.x() - anchor.x, point.y() - anchor.y) < 20.0
    ]
    assert base_distances
    assert max(base_distances) <= 7.0

    bold_neighbors = [
        bond
        for bond in canvas.model.bonds.values()
        if bond.style == BondStyle.BOLD
        and anchor.id in {bond.a1_id, bond.a2_id}
    ]
    assert len(bold_neighbors) == 2
    wedge_shape = wedge_item.mapToScene(wedge_item.path())
    assert all(
        wedge_shape.intersects(canvas.bond_items[bond.id].mapToScene(canvas.bond_items[bond.id].shape()))
        for bond in bold_neighbors
    )
    tip_item = canvas.atom_items[wedge_bond.a1_id]
    if tip_item.label.isVisible():
        assert not wedge_shape.intersects(tip_item.label.mapToScene(tip_item.label.shape()))
    assert atom_snapshot == tuple(
        (atom.id, atom.element, atom.x, atom.y)
        for atom in canvas.model.atoms.values()
    )
    assert bond_snapshot == tuple(
        (bond.id, bond.a1_id, bond.a2_id, bond.order, bond.style, bond.stereo)
        for bond in canvas.model.bonds.values()
    )
    canvas.close()
    app.processEvents()


def test_focused_comparison_sheet_keeps_before_after_theme_panels_equal_scale(tmp_path):
    app = QApplication.instance() or QApplication([])
    before = tmp_path / "before"
    after = tmp_path / "after"
    before.mkdir()
    after.mkdir()
    for theme in ("light", "dark"):
        for directory, color in ((before, "#FF0000"), (after, "#0000FF")):
            panel = QImage(WIDTH, HEIGHT, QImage.Format.Format_ARGB32_Premultiplied)
            panel.fill(QColor(color))
            assert panel.save(str(directory / f"wedge-double-labels-{theme}.png"))

    output = tmp_path / "comparison.png"
    compose_comparison(before, after, output)
    sheet = QImage(str(output))

    assert sheet.width() == 2 * WIDTH
    assert sheet.height() == 2 * (HEIGHT + 44)
    assert sheet.pixelColor(20, 44).name() == "#ff0000"
    assert sheet.pixelColor(WIDTH + 20, 44).name() == "#0000ff"
    app.processEvents()
