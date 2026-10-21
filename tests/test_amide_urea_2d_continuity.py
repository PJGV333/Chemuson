from __future__ import annotations

import math

import pytest
from PyQt6.QtCore import QPointF, QRectF, Qt
from PyQt6.QtGui import QImage, QPainter, QPainterPath
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QApplication

from chemuson.chemio.persistence import PersistenceManager
from chemuson.core.model import BondStyle, BondStereo
from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.geom import angle_deg, angle_distance_deg, endpoint_from_angle_len
from chemuson.gui.theme import apply_theme


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def _show_canvas(canvas: ChemusonCanvas, zoom: float = 1.0) -> None:
    canvas.resize(1000, 760)
    canvas.show()
    QApplication.processEvents()
    if zoom != 1.0:
        canvas.scale(zoom, zoom)
    canvas.center_on_paper()
    QApplication.processEvents()
    canvas.state.fixed_angles = True
    canvas.state.fixed_lengths = True
    canvas.state.angle_step_deg = 30
    canvas.state.bond_length = 40.0


def _bond_tool(canvas: ChemusonCanvas, element: str, order: int = 1) -> None:
    canvas.set_current_tool("tool_bond")
    canvas.state.default_element = element
    canvas.state.active_bond_order = order
    canvas.state.active_bond_style = BondStyle.PLAIN
    canvas.state.active_bond_stereo = BondStereo.NONE
    canvas.state.active_bond_aromatic = False


def _short_bond(canvas: ChemusonCanvas, anchor_id: int, element: str, order: int = 1) -> int:
    _bond_tool(canvas, element, order)
    anchor = canvas.model.get_atom(anchor_id)
    pos = canvas.mapFromScene(QPointF(anchor.x, anchor.y))
    before = set(canvas.model.atoms)
    QTest.mousePress(canvas.viewport(), Qt.MouseButton.LeftButton, pos=pos)
    QApplication.processEvents()
    preview = canvas._preview_bond_item.path()
    assert preview.elementCount() == 2
    preview_endpoint = QPointF(preview.elementAt(1).x, preview.elementAt(1).y)
    QTest.mouseRelease(canvas.viewport(), Qt.MouseButton.LeftButton, pos=pos)
    QApplication.processEvents()
    created = set(canvas.model.atoms) - before
    assert len(created) == 1
    new_id = created.pop()
    atom = canvas.model.get_atom(new_id)
    assert (atom.x, atom.y) == pytest.approx(
        (preview_endpoint.x(), preview_endpoint.y()), abs=1e-7
    )
    return new_id


def _directed_bond(
    canvas: ChemusonCanvas,
    anchor_id: int,
    element: str,
    angle: float,
    *,
    order: int = 1,
) -> tuple[int, QPointF]:
    _bond_tool(canvas, element, order)
    anchor_atom = canvas.model.get_atom(anchor_id)
    anchor = QPointF(anchor_atom.x, anchor_atom.y)
    end = endpoint_from_angle_len(anchor, angle, canvas.state.bond_length)
    p0 = canvas.mapFromScene(anchor)
    p1 = canvas.mapFromScene(end)
    before = set(canvas.model.atoms)
    QTest.mousePress(canvas.viewport(), Qt.MouseButton.LeftButton, pos=p0)
    QTest.mouseMove(canvas.viewport(), p1, 10)
    QApplication.processEvents()
    preview = canvas._preview_bond_item.path()
    preview_end = QPointF(preview.elementAt(1).x, preview.elementAt(1).y)
    QTest.mouseRelease(canvas.viewport(), Qt.MouseButton.LeftButton, pos=p1)
    QApplication.processEvents()
    created = set(canvas.model.atoms) - before
    assert len(created) == 1
    new_id = created.pop()
    atom = canvas.model.get_atom(new_id)
    assert (atom.x, atom.y) == pytest.approx(
        (preview_end.x(), preview_end.y()), abs=1e-7
    )
    return new_id, preview_end


def _draw_benzene_ring(canvas: ChemusonCanvas) -> int:
    _show_canvas(canvas)
    canvas.state.active_ring_size = 6
    canvas.state.active_ring_aromatic = True
    canvas.set_current_tool("tool_ring")
    center = canvas.mapToScene(canvas.viewport().rect().center())
    ring_end = endpoint_from_angle_len(center, 0.0, canvas.state.bond_length)
    p0 = canvas.mapFromScene(center)
    p1 = canvas.mapFromScene(ring_end)
    QTest.mousePress(canvas.viewport(), Qt.MouseButton.LeftButton, pos=p0)
    QTest.mouseMove(canvas.viewport(), p1, 10)
    QTest.mouseRelease(canvas.viewport(), Qt.MouseButton.LeftButton, pos=p1)
    QApplication.processEvents()
    ring_ids = set(canvas.model.atoms)
    assert len(ring_ids) == 6
    return max(
        (canvas.model.get_atom(atom_id) for atom_id in ring_ids),
        key=lambda atom: atom.x,
    ).id


def _draw_reference_urea(canvas: ChemusonCanvas) -> dict[str, int]:
    ring_anchor_id = _draw_benzene_ring(canvas)
    n_aryl = _short_bond(canvas, ring_anchor_id, "N")
    carbonyl = _short_bond(canvas, n_aryl, "C")
    oxygen = _short_bond(canvas, carbonyl, "O", order=2)
    n_urea = _short_bond(canvas, carbonyl, "N")
    methyl = _short_bond(canvas, n_urea, "C")
    canvas.state.default_element = "C"
    return {
        "aryl_carbon": ring_anchor_id,
        "n_aryl": n_aryl,
        "carbonyl": carbonyl,
        "oxygen": oxygen,
        "n_urea": n_urea,
        "methyl": methyl,
    }


def _bond_between(canvas: ChemusonCanvas, atom_a: int, atom_b: int):
    return next(
        bond
        for bond in canvas.model.bonds.values()
        if {bond.a1_id, bond.a2_id} == {atom_a, atom_b}
    )


def _bond_item(canvas: ChemusonCanvas, atom_a: int, atom_b: int):
    return canvas.bond_items[_bond_between(canvas, atom_a, atom_b).id]


def _bond_path_length(item) -> float:
    lengths = []
    for polygon in item.path().toSubpathPolygons():
        lengths.append(
            sum(
                math.hypot(b.x() - a.x(), b.y() - a.y())
                for a, b in zip(polygon, polygon[1:])
            )
        )
    return max(lengths, default=0.0)


def _graph_snapshot(canvas: ChemusonCanvas) -> tuple:
    return (
        tuple(
            (
                atom.id,
                atom.element,
                atom.x,
                atom.y,
                atom.charge,
                atom.isotope,
                atom.explicit_h,
                atom.group_h_cap,
                atom.no_implicit,
            )
            for atom in sorted(canvas.model.atoms.values(), key=lambda item: item.id)
        ),
        tuple(
            (
                bond.id,
                bond.a1_id,
                bond.a2_id,
                bond.order,
                bond.style,
                bond.stereo,
                bond.is_aromatic,
            )
            for bond in sorted(canvas.model.bonds.values(), key=lambda item: item.id)
        ),
    )


def _make_carbonyl_nitrogen_context(
    canvas: ChemusonCanvas,
    *,
    urea: bool = False,
    charge: int = 0,
) -> tuple[int, int | None]:
    carbon = canvas.model.add_atom("C", 0.0, 0.0, is_explicit=True).id
    oxygen = canvas.model.add_atom("O", 40.0, 0.0, is_explicit=True).id
    nitrogen = canvas.model.add_atom("N", -20.0, 34.0, charge=charge, is_explicit=True).id
    canvas.model.add_bond(carbon, oxygen, order=2)
    canvas.model.add_bond(carbon, nitrogen, order=1)
    second_nitrogen = None
    if urea:
        second_nitrogen = canvas.model.add_atom("N", -20.0, -34.0, is_explicit=True).id
        canvas.model.add_bond(carbon, second_nitrogen, order=1)
    canvas._rebuild_items_from_model()
    return nitrogen, second_nitrogen


def test_nitrogen_context_comes_from_graph_connectivity_and_bond_orders() -> None:
    canvas = ChemusonCanvas()
    amide_n, _ = _make_carbonyl_nitrogen_context(canvas)
    assert canvas._nitrogen_geometry_context(amide_n) == "amide"
    assert canvas._bond_geometry(amide_n, 1, False) == "sp2"
    carbonyl_bond = next(
        bond for bond in canvas.model.bonds.values() if bond.order == 2
    )
    carbonyl_bond.order = 1
    assert canvas._nitrogen_geometry_context(amide_n) == "amine"
    carbonyl_bond.order = 2

    urea = ChemusonCanvas()
    urea_n, second_n = _make_carbonyl_nitrogen_context(urea, urea=True)
    assert second_n is not None
    assert urea._nitrogen_geometry_context(urea_n) == "urea"
    assert urea._nitrogen_geometry_context(second_n) == "urea"
    assert urea._bond_geometry(second_n, 1, False) == "sp2"

    amine = ChemusonCanvas()
    nitrogen = amine.model.add_atom("N", 0.0, 0.0, is_explicit=True).id
    carbon = amine.model.add_atom("C", 40.0, 0.0, is_explicit=True).id
    amine.model.add_bond(nitrogen, carbon, order=1)
    amine._rebuild_items_from_model()
    assert amine._nitrogen_geometry_context(nitrogen) == "amine"
    assert amine._bond_geometry(nitrogen, 1, False) == "sp3"

    aromatic = ChemusonCanvas()
    aromatic_n = aromatic.model.add_atom("N", 0.0, 0.0, is_explicit=True).id
    aromatic_left = aromatic.model.add_atom("C", -20.0, 34.0, is_explicit=True).id
    aromatic_right = aromatic.model.add_atom("C", -20.0, -34.0, is_explicit=True).id
    aromatic.model.add_bond(aromatic_n, aromatic_left, is_aromatic=True)
    aromatic.model.add_bond(aromatic_n, aromatic_right, is_aromatic=True)
    aromatic._rebuild_items_from_model()
    assert aromatic._nitrogen_geometry_context(aromatic_n) == "aromatic"
    assert aromatic._bond_geometry(aromatic_n, 1, False) == "sp2"

    coordinated = ChemusonCanvas()
    coordination_n = coordinated.model.add_atom(
        "N", 0.0, 0.0, is_coordination_center=True, is_explicit=True
    ).id
    assert coordinated._nitrogen_geometry_context(coordination_n) == "coordination"
    assert coordinated._bond_geometry(coordination_n, 1, False) == "sp3"

    charged = ChemusonCanvas()
    charged_n, _ = _make_carbonyl_nitrogen_context(charged, charge=1)
    assert charged._nitrogen_geometry_context(charged_n) == "other"
    assert charged._bond_geometry(charged_n, 1, False) == "sp3"

    overcoordinated = ChemusonCanvas()
    over_n = overcoordinated.model.add_atom("N", 0.0, 0.0, is_explicit=True).id
    for index in range(4):
        angle = 90.0 * index
        point = endpoint_from_angle_len(QPointF(0.0, 0.0), angle, 40.0)
        carbon = overcoordinated.model.add_atom(
            "C", point.x(), point.y(), is_explicit=True
        ).id
        overcoordinated.model.add_bond(over_n, carbon, order=1)
    assert overcoordinated._nitrogen_geometry_context(over_n) == "other"


def test_qtest_urea_geometry_keeps_reference_lengths_and_planar_angles() -> None:
    canvas = ChemusonCanvas()
    ids = _draw_reference_urea(canvas)

    assert len(canvas.model.atoms) == 11
    assert len(canvas.model.bonds) == 11
    assert canvas.model.get_atom(ids["n_aryl"]).element == "N"
    assert canvas._nitrogen_geometry_context(ids["n_aryl"]) == "urea"
    assert canvas._nitrogen_geometry_context(ids["n_urea"]) == "urea"
    assert canvas.model.implicit_h_count(ids["n_aryl"]) == 1
    assert canvas.model.implicit_h_count(ids["n_urea"]) == 1
    assert canvas.model.get_atom(ids["oxygen"]).element == "O"

    distances = []
    for bond in canvas.model.bonds.values():
        atom1 = canvas.model.get_atom(bond.a1_id)
        atom2 = canvas.model.get_atom(bond.a2_id)
        distances.append(math.hypot(atom2.x - atom1.x, atom2.y - atom1.y))
    assert distances == pytest.approx([40.0] * 11, abs=1e-7)

    for atom_id in (ids["carbonyl"], ids["n_aryl"], ids["n_urea"]):
        atom = canvas.model.get_atom(atom_id)
        directions = []
        for bond in canvas.model.bonds.values():
            if atom_id not in (bond.a1_id, bond.a2_id):
                continue
            other_id = bond.a2_id if bond.a1_id == atom_id else bond.a1_id
            other = canvas.model.get_atom(other_id)
            directions.append(angle_deg(QPointF(atom.x, atom.y), QPointF(other.x, other.y)))
        if len(directions) >= 2:
            assert angle_distance_deg(*directions[:2]) == pytest.approx(120.0, abs=1e-6)


def test_amide_nh_trimming_uses_glyph_intersection_not_projection() -> None:
    canvas = ChemusonCanvas()
    ids = _draw_reference_urea(canvas)
    snapshot = _graph_snapshot(canvas)

    n_carbonyl = _bond_item(canvas, ids["n_aryl"], ids["carbonyl"])
    carbonyl_n = _bond_item(canvas, ids["carbonyl"], ids["n_urea"])
    aryl_n = _bond_item(canvas, ids["aryl_carbon"], ids["n_aryl"])
    n_methyl = _bond_item(canvas, ids["n_urea"], ids["methyl"])

    amide_trims = [
        n_carbonyl._label_shrink_start + n_carbonyl._label_shrink_end,
        carbonyl_n._label_shrink_start + carbonyl_n._label_shrink_end,
    ]
    assert amide_trims == pytest.approx([8.5645, 8.5645], abs=0.12)
    assert _bond_path_length(n_carbonyl) > 31.0
    assert _bond_path_length(carbonyl_n) > 31.0
    for atom_id, bond_item in (
        (ids["n_aryl"], n_carbonyl),
        (ids["n_urea"], carbonyl_n),
    ):
        label_item = canvas.atom_items[atom_id].label
        label_path = label_item.mapToScene(label_item.shape())
        bond_path = bond_item.mapToScene(bond_item.shape())
        assert not bond_path.intersects(label_path)
    assert aryl_n._label_shrink_start + aryl_n._label_shrink_end == pytest.approx(
        6.4875, abs=0.08
    )
    assert n_methyl._label_shrink_start + n_methyl._label_shrink_end == pytest.approx(
        6.4875, abs=0.12
    )
    assert _graph_snapshot(canvas) == snapshot


def test_explicit_drag_direction_is_kept_for_amide_n_substituent() -> None:
    canvas = ChemusonCanvas()
    ids = _draw_reference_urea(canvas)
    canvas.state.fixed_angles = False
    anchor = canvas.model.get_atom(ids["n_urea"])
    drag_angle = 37.0
    new_id, preview = _directed_bond(canvas, ids["n_urea"], "C", drag_angle)
    created = canvas.model.get_atom(new_id)
    actual = angle_deg(QPointF(anchor.x, anchor.y), QPointF(created.x, created.y))
    assert actual == pytest.approx(drag_angle, abs=0.6)
    assert (created.x, created.y) == pytest.approx((preview.x(), preview.y()), abs=1e-7)
    assert math.hypot(created.x - anchor.x, created.y - anchor.y) == pytest.approx(40.0, abs=1e-7)


def test_qtest_urea_undo_redo_and_cmsn_roundtrip_preserve_graph(tmp_path) -> None:
    canvas = ChemusonCanvas()
    ids = _draw_reference_urea(canvas)
    expected = _graph_snapshot(canvas)
    expected_h_positions = {
        atom_id: canvas.atom_items[atom_id].implicit_hydrogen_position
        for atom_id in (ids["n_aryl"], ids["n_urea"])
    }

    canvas.undo_stack.undo()
    assert ids["methyl"] not in canvas.model.atoms
    canvas.undo_stack.redo()
    assert _graph_snapshot(canvas) == expected
    assert {
        atom_id: canvas.atom_items[atom_id].implicit_hydrogen_position
        for atom_id in expected_h_positions
    } == expected_h_positions

    path = tmp_path / "amide-urea.cmsn"
    PersistenceManager.save_to_file(str(path), canvas)
    restored = ChemusonCanvas()
    PersistenceManager.load_from_file(str(path), restored)
    assert _graph_snapshot(restored) == expected
    assert {
        atom_id: restored.atom_items[atom_id].implicit_hydrogen_position
        for atom_id in expected_h_positions
    } == expected_h_positions
    restored_n_carbonyl = _bond_item(restored, ids["n_aryl"], ids["carbonyl"])
    assert restored_n_carbonyl._label_shrink_start + restored_n_carbonyl._label_shrink_end == pytest.approx(
        8.5645, abs=0.12
    )


@pytest.mark.parametrize(
    ("zoom", "dpr", "theme"),
    [(1.0, 1.0, "light"), (1.25, 1.25, "dark"), (1.5, 2.0, "light")],
)
def test_amide_label_clearance_is_stable_across_zoom_dpr_and_theme(
    tmp_path, zoom: float, dpr: float, theme: str
) -> None:
    canvas = ChemusonCanvas()
    ids = _draw_reference_urea(canvas)
    assert apply_theme(canvas, theme) == theme
    canvas.resetTransform()
    if zoom != 1.0:
        canvas.scale(zoom, zoom)
    QApplication.processEvents()

    item = _bond_item(canvas, ids["n_aryl"], ids["carbonyl"])
    shrink = item._label_shrink_start + item._label_shrink_end
    assert shrink == pytest.approx(8.5645, abs=0.12)

    image = QImage(500, 360, QImage.Format.Format_ARGB32_Premultiplied)
    image.setDevicePixelRatio(dpr)
    image.fill(Qt.GlobalColor.white)
    painter = QPainter(image)
    center = canvas.mapToScene(canvas.viewport().rect().center())
    canvas.scene.render(
        painter,
        QRectF(0.0, 0.0, 500.0, 360.0),
        QRectF(center.x() - 100.0, center.y() - 80.0, 200.0, 160.0),
    )
    painter.end()
    assert image.save(str(tmp_path / f"amide-urea-{theme}-{zoom:g}x-dpr-{dpr:g}.png"))


def test_nh2_and_nh_labels_follow_implicit_hydrogen_count() -> None:
    canvas = ChemusonCanvas()
    carbonyl = canvas.model.add_atom("C", 0.0, 0.0, is_explicit=True).id
    oxygen = canvas.model.add_atom("O", 40.0, 0.0, is_explicit=True).id
    nitrogen = canvas.model.add_atom("N", -20.0, 34.0, is_explicit=True).id
    canvas.model.add_bond(carbonyl, oxygen, order=2)
    canvas.model.add_bond(carbonyl, nitrogen, order=1)
    canvas._rebuild_items_from_model()
    label, _anchor, _offset = canvas._build_atom_label(canvas.model.get_atom(nitrogen))
    assert label in {"NH2", "H2N"}
    assert canvas.model.implicit_h_count(nitrogen) == 2

    extra = canvas.model.add_atom("C", -60.0, 34.0, is_explicit=True).id
    canvas.model.add_bond(nitrogen, extra, order=1)
    canvas._rebuild_items_from_model()
    label, _anchor, _offset = canvas._build_atom_label(canvas.model.get_atom(nitrogen))
    assert label in {"NH", "HN"}
    assert canvas.model.implicit_h_count(nitrogen) == 1


@pytest.mark.parametrize(
    ("degree", "expected_label"),
    [(1, "NH2"), (2, "NH"), (3, "N")],
    ids=("primary-amine", "secondary-amine", "tertiary-amine"),
)
def test_primary_secondary_and_tertiary_amine_contexts_stay_distinct_from_amides(
    degree: int, expected_label: str
) -> None:
    canvas = ChemusonCanvas()
    nitrogen = canvas.model.add_atom("N", 0.0, 0.0, is_explicit=True).id
    for index in range(degree):
        angle = 360.0 * index / max(1, degree)
        point = endpoint_from_angle_len(QPointF(0.0, 0.0), angle, 40.0)
        carbon = canvas.model.add_atom(
            "C", point.x(), point.y(), is_explicit=True
        ).id
        canvas.model.add_bond(nitrogen, carbon, order=1)
    canvas._rebuild_items_from_model()

    label, _anchor, _offset = canvas._build_atom_label(
        canvas.model.get_atom(nitrogen)
    )
    alternate_label = {"NH2": "H2N", "NH": "HN", "N": "N"}[expected_label]
    assert label in {expected_label, alternate_label}
    assert canvas._nitrogen_geometry_context(nitrogen) == "amine"
    assert canvas._bond_geometry(nitrogen, 1, False) == "sp3"


def test_qtest_acetanilide_amide_keeps_trigonal_carbonyl_and_nitrogen() -> None:
    canvas = ChemusonCanvas()
    aryl_carbon = _draw_benzene_ring(canvas)
    nitrogen = _short_bond(canvas, aryl_carbon, "N")
    carbonyl = _short_bond(canvas, nitrogen, "C")
    oxygen = _short_bond(canvas, carbonyl, "O", order=2)
    methyl = _short_bond(canvas, carbonyl, "C")

    assert canvas._nitrogen_geometry_context(nitrogen) == "amide"
    assert canvas.model.implicit_h_count(nitrogen) == 1
    assert canvas.model.get_atom(oxygen).element == "O"
    assert canvas.model.get_atom(methyl).element == "C"
    assert math.hypot(
        canvas.model.get_atom(oxygen).x - canvas.model.get_atom(carbonyl).x,
        canvas.model.get_atom(oxygen).y - canvas.model.get_atom(carbonyl).y,
    ) == pytest.approx(40.0, abs=1e-7)

    carbonyl_angles = []
    center = canvas.model.get_atom(carbonyl)
    for bond in canvas.model.bonds.values():
        if carbonyl not in (bond.a1_id, bond.a2_id):
            continue
        neighbor_id = bond.a2_id if bond.a1_id == carbonyl else bond.a1_id
        neighbor = canvas.model.get_atom(neighbor_id)
        carbonyl_angles.append(
            angle_deg(QPointF(center.x, center.y), QPointF(neighbor.x, neighbor.y))
        )
    assert [
        angle_distance_deg(carbonyl_angles[0], angle)
        for angle in carbonyl_angles[1:]
    ] == pytest.approx([120.0, 120.0], abs=1e-6)

    _assert_nh_glyph_geometry(canvas, nitrogen, "below")
    n_c = _bond_item(canvas, nitrogen, carbonyl)
    assert n_c._label_shrink_start + n_c._label_shrink_end < 10.0
    assert _bond_path_length(n_c) > 30.0


def test_qtest_aniline_and_aromatic_ring_are_unchanged_controls() -> None:
    canvas = ChemusonCanvas()
    aryl_carbon = _draw_benzene_ring(canvas)
    ring_bonds = list(canvas.model.bonds.values())
    ring_coordinates = {
        atom_id: (canvas.model.get_atom(atom_id).x, canvas.model.get_atom(atom_id).y)
        for atom_id in canvas.model.atoms
    }
    nitrogen = _short_bond(canvas, aryl_carbon, "N")
    methyl = _short_bond(canvas, nitrogen, "C")

    assert canvas._nitrogen_geometry_context(nitrogen) == "amine"
    assert canvas.model.implicit_h_count(nitrogen) == 1
    assert canvas.model.get_atom(methyl).element == "C"
    assert angle_distance_deg(
        angle_deg(
            QPointF(canvas.model.get_atom(nitrogen).x, canvas.model.get_atom(nitrogen).y),
            QPointF(canvas.model.get_atom(aryl_carbon).x, canvas.model.get_atom(aryl_carbon).y),
        ),
        angle_deg(
            QPointF(canvas.model.get_atom(nitrogen).x, canvas.model.get_atom(nitrogen).y),
            QPointF(canvas.model.get_atom(methyl).x, canvas.model.get_atom(methyl).y),
        ),
    ) == pytest.approx(120.0, abs=1e-6)

    for bond in ring_bonds:
        assert bond.is_aromatic is True
        item = canvas.bond_items[bond.id]
        assert item._label_shrink_start == pytest.approx(0.0)
        assert item._label_shrink_end == pytest.approx(0.0)
        first = canvas.model.get_atom(bond.a1_id)
        second = canvas.model.get_atom(bond.a2_id)
        assert math.hypot(second.x - first.x, second.y - first.y) == pytest.approx(
            40.0, abs=1e-7
        )
    assert {
        atom_id: (canvas.model.get_atom(atom_id).x, canvas.model.get_atom(atom_id).y)
        for atom_id in ring_coordinates
    } == ring_coordinates


def _first_glyph_scene_rect(text_item) -> QRectF:
    block = text_item.document().begin()
    line = block.layout().lineForTextPosition(0)
    assert line.isValid()
    x = line.cursorToX(0)
    if isinstance(x, tuple):
        x = x[0]
    cursor = text_item.textCursor()
    cursor.setPosition(0)
    cursor.setPosition(1, cursor.MoveMode.KeepAnchor)
    path = QPainterPath()
    block_rect = text_item.document().documentLayout().blockBoundingRect(block)
    path.addText(
        block_rect.x() + line.position().x() + float(x),
        block_rect.y() + line.position().y()
        + line.ascent(),
        cursor.charFormat().font(),
        block.text()[0],
    )
    return text_item.mapToScene(path).boundingRect()


def _assert_nh_glyph_geometry(canvas: ChemusonCanvas, atom_id: int, position: str) -> None:
    atom = canvas.model.get_atom(atom_id)
    atom_item = canvas.atom_items[atom_id]
    n_rect = _first_glyph_scene_rect(atom_item.label)
    h_rect = _first_glyph_scene_rect(atom_item.implicit_hydrogen_label)
    assert atom_item.implicit_hydrogen_position == position
    assert atom_item.label.isVisible()
    assert atom_item.implicit_hydrogen_label.isVisible()
    assert atom_item.implicit_hydrogen_label.rotation() == pytest.approx(0.0)
    assert n_rect.width() >= 5.0 and n_rect.height() >= 8.0
    assert h_rect.width() >= 5.0 and h_rect.height() >= 8.0
    n_cursor = atom_item.label.textCursor()
    n_cursor.setPosition(0)
    n_cursor.setPosition(1, n_cursor.MoveMode.KeepAnchor)
    h_cursor = atom_item.implicit_hydrogen_label.textCursor()
    h_cursor.setPosition(0)
    h_cursor.setPosition(1, h_cursor.MoveMode.KeepAnchor)
    assert n_cursor.charFormat().foreground().color().name() == "#3050f8"
    assert h_cursor.charFormat().foreground().color().name() == "#000000"
    assert n_rect.center().x() == pytest.approx(atom.x, abs=0.1)
    assert n_rect.center().y() == pytest.approx(atom.y, abs=0.1)
    assert not n_rect.intersects(h_rect)
    union_rect = atom_item.label.mapToScene(atom_item.label.shape()).boundingRect()
    assert union_rect.contains(n_rect)
    assert union_rect.contains(h_rect)
    if position == "below":
        assert h_rect.top() - n_rect.bottom() >= 1.0
    elif position == "above":
        assert n_rect.top() - h_rect.bottom() >= 1.0
    elif position == "right":
        assert h_rect.left() - n_rect.right() >= 1.0
    else:
        assert n_rect.left() - h_rect.right() >= 1.0


@pytest.mark.parametrize("theme", ["light", "dark"])
def test_qtest_urea_places_hydrogens_outside_bond_sectors_and_keeps_clearance(
    theme: str,
) -> None:
    canvas = ChemusonCanvas()
    ids = _draw_reference_urea(canvas)
    assert apply_theme(canvas, theme) == theme
    _assert_nh_glyph_geometry(canvas, ids["n_aryl"], "below")
    _assert_nh_glyph_geometry(canvas, ids["n_urea"], "above")

    for nitrogen_id in (ids["n_aryl"], ids["n_urea"]):
        label_path = canvas.atom_items[nitrogen_id].label.mapToScene(
            canvas.atom_items[nitrogen_id].label.shape()
        )
        for bond in canvas.model.bonds.values():
            if nitrogen_id not in (bond.a1_id, bond.a2_id):
                continue
            bond_item = canvas.bond_items[bond.id]
            assert not label_path.intersects(bond_item.mapToScene(bond_item.shape()))
        for other_id, other_item in canvas.atom_items.items():
            if other_id == nitrogen_id or not other_item.label.isVisible():
                continue
            other_path = other_item.label.mapToScene(other_item.label.shape())
            assert not label_path.intersects(other_path)

    distances = [
        math.hypot(
            canvas.model.get_atom(bond.a2_id).x - canvas.model.get_atom(bond.a1_id).x,
            canvas.model.get_atom(bond.a2_id).y - canvas.model.get_atom(bond.a1_id).y,
        )
        for bond in canvas.model.bonds.values()
    ]
    assert distances == pytest.approx([40.0] * len(distances), abs=1e-7)


@pytest.mark.parametrize(
    ("degree", "expected_hydrogens", "expected_position"),
    [(1, 2, "left"), (2, 1, "above")],
    ids=("primary-amine", "secondary-amine"),
)
def test_qtest_amine_nh_glyphs_reposition_without_moving_the_nitrogen(
    degree: int, expected_hydrogens: int, expected_position: str
) -> None:
    canvas = ChemusonCanvas()
    _show_canvas(canvas)
    canvas.set_current_tool("tool_atom")
    canvas.state.default_element = "N"
    center = canvas.mapToScene(canvas.viewport().rect().center())
    pos = canvas.mapFromScene(center)
    before = set(canvas.model.atoms)
    QTest.mousePress(canvas.viewport(), Qt.MouseButton.LeftButton, pos=pos)
    QTest.mouseRelease(canvas.viewport(), Qt.MouseButton.LeftButton, pos=pos)
    QApplication.processEvents()
    nitrogen_ids = set(canvas.model.atoms) - before
    assert len(nitrogen_ids) == 1
    nitrogen_id = nitrogen_ids.pop()

    for _ in range(degree):
        _short_bond(canvas, nitrogen_id, "C")
    nitrogen = canvas.model.get_atom(nitrogen_id)
    assert canvas.model.implicit_h_count(nitrogen_id) == expected_hydrogens
    expected_label = "H2N" if degree == 1 else "HN"
    assert canvas._build_atom_label(nitrogen)[0] == expected_label
    item = canvas.atom_items[nitrogen_id]
    assert item.implicit_hydrogen_position == expected_position
    _assert_nh_glyph_geometry(canvas, nitrogen_id, expected_position)
    assert (nitrogen.x, nitrogen.y) == pytest.approx((center.x(), center.y()), abs=1e-7)


def test_qtest_primary_amine_suffix_h2_label_places_hydrogen_on_free_side() -> None:
    canvas = ChemusonCanvas()
    _show_canvas(canvas)
    canvas.set_current_tool("tool_atom")
    canvas.state.default_element = "N"
    center = canvas.mapToScene(canvas.viewport().rect().center())
    pos = canvas.mapFromScene(center)
    QTest.mousePress(canvas.viewport(), Qt.MouseButton.LeftButton, pos=pos)
    QTest.mouseRelease(canvas.viewport(), Qt.MouseButton.LeftButton, pos=pos)
    QApplication.processEvents()
    nitrogen_id = next(iter(canvas.model.atoms))

    _directed_bond(canvas, nitrogen_id, "C", 180.0)

    assert canvas._build_atom_label(canvas.model.get_atom(nitrogen_id))[0] == "NH2"
    _assert_nh_glyph_geometry(canvas, nitrogen_id, "right")
    nitrogen = canvas.model.get_atom(nitrogen_id)
    assert (nitrogen.x, nitrogen.y) == pytest.approx((center.x(), center.y()), abs=1e-7)
