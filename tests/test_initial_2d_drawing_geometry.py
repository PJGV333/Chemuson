from __future__ import annotations

import math

import pytest
from PyQt6.QtCore import QPointF, QRectF, Qt
from PyQt6.QtGui import QImage, QPainter
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QApplication

from chemuson.chemio.persistence import PersistenceManager
from chemuson.chemio.rdkit_io import molgraph_to_smiles_isolated_or_error
from chemuson.core.model import BondStyle, BondStereo
from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.geom import (
    SP3_BOND_ANGLE_DEG,
    SP3_SKELETAL_ANGLE_2D_DEG,
    angle_deg,
    angle_distance_deg,
    candidate_directions_deg,
    endpoint_from_angle_len,
)


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def _show_canvas(canvas: ChemusonCanvas, *, view_scale: float = 1.0) -> None:
    canvas.resize(900, 650)
    canvas.show()
    QApplication.processEvents()
    if view_scale != 1.0:
        canvas.scale(view_scale, view_scale)
    canvas.center_on_paper()
    QApplication.processEvents()


def _gesture(
    canvas: ChemusonCanvas,
    start: QPointF,
    end: QPointF,
) -> None:
    start_pos = canvas.mapFromScene(start)
    end_pos = canvas.mapFromScene(end)
    assert canvas.viewport().rect().contains(start_pos)
    assert canvas.viewport().rect().contains(end_pos)
    QTest.mousePress(canvas.viewport(), Qt.MouseButton.LeftButton, pos=start_pos)
    QTest.mouseMove(canvas.viewport(), end_pos, 10)
    QTest.mouseRelease(canvas.viewport(), Qt.MouseButton.LeftButton, pos=end_pos)
    QApplication.processEvents()


def _draw_free_ring(
    canvas: ChemusonCanvas,
    *,
    size: int = 6,
    angle_deg_value: float = 0.0,
    view_scale: float = 1.0,
) -> tuple[int, ...]:
    _show_canvas(canvas, view_scale=view_scale)
    canvas.state.active_ring_size = size
    canvas.state.active_ring_aromatic = size == 6
    canvas.set_current_tool("tool_ring")
    center = canvas.mapToScene(canvas.viewport().rect().center())
    cursor = endpoint_from_angle_len(
        center,
        angle_deg_value,
        canvas.state.bond_length,
    )
    _gesture(canvas, center, cursor)
    return tuple(sorted(canvas.model.atoms))


def _interior_angles(points: list[QPointF]) -> list[float]:
    angles = []
    for index, point in enumerate(points):
        previous = points[(index - 1) % len(points)]
        following = points[(index + 1) % len(points)]
        a = angle_deg(point, previous)
        b = angle_deg(point, following)
        angles.append(angle_distance_deg(a, b))
    return angles


def _cross_product(p: QPointF, q: QPointF, r: QPointF) -> float:
    return (q.x() - p.x()) * (r.y() - p.y()) - (q.y() - p.y()) * (r.x() - p.x())


def _non_adjacent_crossing_count(canvas: ChemusonCanvas) -> int:
    crossings = 0
    bonds = list(canvas.model.bonds.values())
    for index, first in enumerate(bonds):
        p0 = canvas.model.get_atom(first.a1_id)
        p1 = canvas.model.get_atom(first.a2_id)
        a = QPointF(p0.x, p0.y)
        b = QPointF(p1.x, p1.y)
        for second in bonds[index + 1 :]:
            if {first.a1_id, first.a2_id} & {second.a1_id, second.a2_id}:
                continue
            q0 = canvas.model.get_atom(second.a1_id)
            q1 = canvas.model.get_atom(second.a2_id)
            c = QPointF(q0.x, q0.y)
            d = QPointF(q1.x, q1.y)
            if (
                _cross_product(a, b, c) * _cross_product(a, b, d) < -1e-8
                and _cross_product(c, d, a) * _cross_product(c, d, b) < -1e-8
            ):
                crossings += 1
    return crossings


def _regular_cycle_points(canvas: ChemusonCanvas, ring_ids: set[int]) -> list[QPointF]:
    adjacency = {atom_id: [] for atom_id in ring_ids}
    for bond in canvas.model.bonds.values():
        if bond.a1_id in ring_ids and bond.a2_id in ring_ids:
            adjacency[bond.a1_id].append(bond.a2_id)
            adjacency[bond.a2_id].append(bond.a1_id)
    assert all(len(neighbors) == 2 for neighbors in adjacency.values())
    start = min(ring_ids)
    ordered = [start]
    previous = None
    current = start
    while len(ordered) < len(ring_ids):
        following = next(atom_id for atom_id in adjacency[current] if atom_id != previous)
        if following == start:
            break
        ordered.append(following)
        previous, current = current, following
    assert len(ordered) == len(ring_ids)
    return [
        QPointF(canvas.model.get_atom(atom_id).x, canvas.model.get_atom(atom_id).y)
        for atom_id in ordered
    ]


def _ring_snapshot(canvas: ChemusonCanvas) -> tuple:
    return (
        tuple(
            (atom_id, canvas.model.get_atom(atom_id).x, canvas.model.get_atom(atom_id).y)
            for atom_id in sorted(canvas.model.atoms)
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
                bond.ring_id,
            )
            for bond in sorted(canvas.model.bonds.values(), key=lambda item: item.id)
        ),
    )


def _short_bond(canvas: ChemusonCanvas, anchor_id: int, *, element: str = "C", order: int = 1) -> int:
    canvas.set_current_tool("tool_bond")
    canvas.state.default_element = element
    canvas.state.active_bond_order = order
    canvas.state.active_bond_style = BondStyle.PLAIN
    canvas.state.active_bond_stereo = BondStereo.NONE
    canvas.state.active_bond_aromatic = False
    anchor = canvas.model.get_atom(anchor_id)
    pos = canvas.mapFromScene(QPointF(anchor.x, anchor.y))
    before = set(canvas.model.atoms)
    QTest.mousePress(canvas.viewport(), Qt.MouseButton.LeftButton, pos=pos)
    QApplication.processEvents()
    preview_path = canvas._preview_bond_item.path()
    assert preview_path.elementCount() == 2
    preview_endpoint = QPointF(
        preview_path.elementAt(1).x,
        preview_path.elementAt(1).y,
    )
    QTest.mouseRelease(canvas.viewport(), Qt.MouseButton.LeftButton, pos=pos)
    QApplication.processEvents()
    created = set(canvas.model.atoms) - before
    assert len(created) == 1, f"short click at atom {anchor_id} created {created!r}"
    new_atom_id = created.pop()
    new_atom = canvas.model.get_atom(new_atom_id)
    assert (new_atom.x, new_atom.y) == pytest.approx(
        (preview_endpoint.x(), preview_endpoint.y()),
        abs=1e-7,
    )
    return new_atom_id


def test_regular_polygon_geometry_is_size_and_rigid_transform_invariant() -> None:
    canvas = ChemusonCanvas()
    side = 40.0
    rotation = 23.0
    shift = QPointF(123.5, -42.25)

    for size in (3, 4, 5, 6):
        start = QPointF(-11.0, 7.5) + shift
        end = endpoint_from_angle_len(start, rotation, side)
        vertices = canvas._regular_polygon_from_edge(start, end, size, direction=1)
        assert len(vertices) == size

        lengths = [
            math.hypot(
                vertices[(index + 1) % size].x() - vertex.x(),
                vertices[(index + 1) % size].y() - vertex.y(),
            )
            for index, vertex in enumerate(vertices)
        ]
        expected_angle = 180.0 - 360.0 / size
        assert lengths == pytest.approx([side] * size, abs=1e-9)
        assert _interior_angles(vertices) == pytest.approx(
            [expected_angle] * size,
            abs=1e-8,
        )


def test_ring_tool_coordinates_are_independent_of_view_scale_and_qimage_dpr(tmp_path) -> None:
    app = QApplication.instance() or QApplication([])
    for view_scale, dpr in ((1.0, 1.0), (1.25, 1.25)):
        canvas = ChemusonCanvas()
        ring_ids = _draw_free_ring(
            canvas,
            angle_deg_value=17.0,
            view_scale=view_scale,
        )
        assert len(ring_ids) == 6
        snapshot = _ring_snapshot(canvas)
        points = [
            QPointF(canvas.model.get_atom(atom_id).x, canvas.model.get_atom(atom_id).y)
            for atom_id in ring_ids
        ]
        lengths = [
            math.hypot(
                points[(index + 1) % 6].x() - point.x(),
                points[(index + 1) % 6].y() - point.y(),
            )
            for index, point in enumerate(points)
        ]
        assert lengths == pytest.approx([40.0] * 6, abs=1e-8)
        assert _interior_angles(points) == pytest.approx([120.0] * 6, abs=1e-8)

        image = QImage(300, 300, QImage.Format.Format_ARGB32_Premultiplied)
        image.setDevicePixelRatio(dpr)
        image.fill(Qt.GlobalColor.white)
        painter = QPainter(image)
        center = canvas.mapToScene(canvas.viewport().rect().center())
        source = QRectF(center.x() - 70.0, center.y() - 70.0, 140.0, 140.0)
        canvas.scene.render(painter, QRectF(0.0, 0.0, 300.0, 300.0), source)
        painter.end()
        assert image.save(str(tmp_path / f"benzene-{view_scale:g}x-dpr-{dpr:g}.png"))
        assert _ring_snapshot(canvas) == snapshot
        canvas.close()
        app.processEvents()


def test_ring_anchored_to_existing_bond_reuses_edge_and_keeps_heteroatom_geometry() -> None:
    canvas = ChemusonCanvas()
    canvas.state.active_ring_size = 6
    canvas.state.active_ring_aromatic = False
    canvas.state.default_element = "C"
    carbon = canvas.model.add_atom("C", 300.0, 300.0)
    oxygen = canvas.model.add_atom("O", 340.0, 300.0)
    shared_edge = canvas.model.add_bond(carbon.id, oxygen.id, order=1)
    canvas._rebuild_items_from_model()
    original = {
        atom_id: (canvas.model.get_atom(atom_id).x, canvas.model.get_atom(atom_id).y)
        for atom_id in (carbon.id, oxygen.id)
    }

    canvas._create_ring_from_bond_hotkey(shared_edge.id, 6)

    ring_ids = {
        atom_id
        for bond in canvas.model.bonds.values()
        if bond.ring_id is not None
        for atom_id in (bond.a1_id, bond.a2_id)
    }
    assert len(ring_ids) == 6
    assert len(canvas.model.atoms) == 6
    assert len(canvas.model.bonds) == 6
    assert sum(bond.id == shared_edge.id for bond in canvas.model.bonds.values()) == 1
    assert canvas.model.get_atom(oxygen.id).element == "O"
    assert {
        atom_id: (canvas.model.get_atom(atom_id).x, canvas.model.get_atom(atom_id).y)
        for atom_id in (carbon.id, oxygen.id)
    } == original

    points = _regular_cycle_points(canvas, ring_ids)
    lengths = [
        math.hypot(
            points[(index + 1) % 6].x() - point.x(),
            points[(index + 1) % 6].y() - point.y(),
        )
        for index, point in enumerate(points)
    ]
    assert lengths == pytest.approx([40.0] * 6, abs=1e-8)
    assert _interior_angles(points) == pytest.approx([120.0] * 6, abs=1e-8)


def test_ring_anchored_to_atom_preserves_existing_substituent() -> None:
    canvas = ChemusonCanvas()
    canvas.state.active_ring_size = 6
    anchor = canvas.model.add_atom("C", 300.0, 300.0)
    substituent_pos = endpoint_from_angle_len(QPointF(anchor.x, anchor.y), 180.0, 40.0)
    substituent = canvas.model.add_atom("C", substituent_pos.x(), substituent_pos.y())
    canvas.model.add_bond(anchor.id, substituent.id)
    canvas._rebuild_items_from_model()

    canvas._create_ring_hotkey(anchor.id, 6)

    assert (
        canvas.model.get_atom(anchor.id).x,
        canvas.model.get_atom(anchor.id).y,
        canvas.model.get_atom(substituent.id).x,
        canvas.model.get_atom(substituent.id).y,
    ) == pytest.approx((300.0, 300.0, substituent_pos.x(), substituent_pos.y()))
    assert canvas.model.find_bond_between(anchor.id, substituent.id) is not None
    ring_ids = {
        atom_id
        for bond in canvas.model.bonds.values()
        if bond.ring_id is not None
        for atom_id in (bond.a1_id, bond.a2_id)
    }
    points = _regular_cycle_points(canvas, ring_ids)
    lengths = [
        math.hypot(
            points[(index + 1) % 6].x() - point.x(),
            points[(index + 1) % 6].y() - point.y(),
        )
        for index, point in enumerate(points)
    ]
    assert lengths == pytest.approx([40.0] * 6, abs=1e-8)
    assert _interior_angles(points) == pytest.approx([120.0] * 6, abs=1e-8)


def test_tetrahedral_legacy_default_and_planar_candidate_override_are_distinct() -> None:
    assert SP3_BOND_ANGLE_DEG == 109.5
    assert SP3_SKELETAL_ANGLE_2D_DEG == 120.0

    legacy = candidate_directions_deg("sp3", [0.0], 0.0)
    planar = candidate_directions_deg("sp3", [0.0], 0.0, sp3_angle_deg=120.0)
    assert legacy == pytest.approx([109.5, 250.5])
    assert planar == pytest.approx([120.0, 240.0])


def test_default_sp3_suggestion_uses_editorial_angle_not_tetrahedral_angle() -> None:
    canvas = ChemusonCanvas()
    canvas.state.fixed_angles = False
    canvas.state.fixed_lengths = True
    center = canvas.model.add_atom("C", 200.0, 200.0)
    p0 = QPointF(center.x, center.y)
    neighbor_pos = endpoint_from_angle_len(p0, 180.0, canvas.state.bond_length)
    neighbor = canvas.model.add_atom("C", neighbor_pos.x(), neighbor_pos.y())
    canvas.model.add_bond(center.id, neighbor.id)
    canvas._rebuild_items_from_model()

    endpoint = canvas._compute_default_bond_endpoint(p0, center.id)
    local_angle = angle_distance_deg(
        angle_deg(p0, endpoint),
        angle_deg(p0, neighbor_pos),
    )
    assert local_angle == pytest.approx(120.0, abs=1e-8)


def test_free_chain_zigzag_has_editorial_120_degree_internal_angles() -> None:
    canvas = ChemusonCanvas()
    canvas.state.fixed_angles = True
    canvas.state.fixed_lengths = True
    canvas.state.angle_step_deg = 30
    start = QPointF(200.0, 220.0)
    canvas._begin_place_chain(None, start)
    canvas._update_chain_preview(
        QPointF(360.0, 220.0),
        Qt.KeyboardModifier.NoModifier,
    )
    points = canvas._chain_last_points
    assert points is not None and len(points) == 5

    lengths = [
        math.hypot(b.x() - a.x(), b.y() - a.y())
        for a, b in zip(points, points[1:])
    ]
    assert lengths == pytest.approx([canvas.state.bond_length] * 4, abs=1e-8)
    chain_angles = [
        angle_distance_deg(
            angle_deg(points[index], points[index - 1]),
            angle_deg(points[index], points[index + 1]),
        )
        for index in range(1, len(points) - 1)
    ]
    assert chain_angles == pytest.approx([120.0] * 3, abs=1e-8)


def test_fixed_angle_step_and_alt_shift_bond_modifiers_keep_existing_semantics() -> None:
    canvas = ChemusonCanvas()
    canvas.state.fixed_angles = True
    canvas.state.fixed_lengths = True
    canvas.state.angle_step_deg = 30
    anchor = canvas.model.add_atom("C", 300.0, 300.0, is_explicit=True)
    canvas._rebuild_items_from_model()
    p0 = QPointF(anchor.x, anchor.y)
    canvas._begin_place_bond(anchor.id, p0)

    requested = endpoint_from_angle_len(p0, 47.0, 40.0)
    snapped = canvas._compute_bond_endpoint(p0, requested, Qt.KeyboardModifier.NoModifier)
    snapped_angle = angle_deg(p0, snapped)
    assert snapped_angle / 30.0 == pytest.approx(round(snapped_angle / 30.0), abs=1e-8)
    assert math.hypot(snapped.x() - p0.x(), snapped.y() - p0.y()) == pytest.approx(40.0)

    alt_endpoint = canvas._compute_bond_endpoint(
        p0,
        requested,
        Qt.KeyboardModifier.AltModifier,
    )
    assert angle_distance_deg(angle_deg(p0, alt_endpoint), 47.0) < 1e-8
    assert math.hypot(alt_endpoint.x() - p0.x(), alt_endpoint.y() - p0.y()) == pytest.approx(40.0)

    longer_cursor = endpoint_from_angle_len(p0, 47.0, 53.0)
    shift_endpoint = canvas._compute_bond_endpoint(
        p0,
        longer_cursor,
        Qt.KeyboardModifier.ShiftModifier,
    )
    assert math.hypot(shift_endpoint.x() - p0.x(), shift_endpoint.y() - p0.y()) == pytest.approx(53.0)
    canvas.close()
    QApplication.processEvents()


def test_short_click_on_rotated_benzene_vertex_adds_radial_bond_without_mutating_ring() -> None:
    app = QApplication.instance() or QApplication([])
    canvas = ChemusonCanvas()
    ring_ids = _draw_free_ring(canvas, angle_deg_value=17.0)
    ring_snapshot = _ring_snapshot(canvas)
    ring_anchor = max(
        (canvas.model.get_atom(atom_id) for atom_id in ring_ids),
        key=lambda atom: atom.x,
    )
    center = QPointF(
        sum(canvas.model.get_atom(atom_id).x for atom_id in ring_ids) / 6.0,
        sum(canvas.model.get_atom(atom_id).y for atom_id in ring_ids) / 6.0,
    )
    anchor_pos = QPointF(ring_anchor.x, ring_anchor.y)
    expected_endpoint = canvas._compute_default_bond_endpoint(anchor_pos, ring_anchor.id)
    before_atom_ids = set(canvas.model.atoms)

    new_atom_id = _short_bond(canvas, ring_anchor.id)

    assert new_atom_id not in before_atom_ids
    assert tuple(
        (
            atom_id,
            canvas.model.get_atom(atom_id).x,
            canvas.model.get_atom(atom_id).y,
        )
        for atom_id in ring_ids
    ) == ring_snapshot[0]
    after_bonds = tuple(
        (
            bond.id,
            bond.a1_id,
            bond.a2_id,
            bond.order,
            bond.style,
            bond.stereo,
            bond.is_aromatic,
            bond.ring_id,
        )
        for bond in sorted(canvas.model.bonds.values(), key=lambda item: item.id)
        if bond.id in {saved[0] for saved in ring_snapshot[1]}
    )
    assert after_bonds == ring_snapshot[1]
    assert len(canvas.model.bonds) == 7
    new_bond = next(
        bond
        for bond in canvas.model.bonds.values()
        if new_atom_id in (bond.a1_id, bond.a2_id)
    )
    assert new_bond.order == 1
    assert not new_bond.is_aromatic

    new_atom = canvas.model.get_atom(new_atom_id)
    new_pos = QPointF(new_atom.x, new_atom.y)
    assert new_pos.x() == pytest.approx(expected_endpoint.x(), abs=1e-7)
    assert new_pos.y() == pytest.approx(expected_endpoint.y(), abs=1e-7)
    radial_angle = angle_deg(center, anchor_pos)
    substituent_angle = angle_deg(anchor_pos, new_pos)
    assert angle_distance_deg(radial_angle, substituent_angle) < 1e-6

    ring_neighbor_angles = []
    for bond in canvas.model.bonds.values():
        if bond.id == new_bond.id or ring_anchor.id not in (bond.a1_id, bond.a2_id):
            continue
        other_id = bond.a2_id if bond.a1_id == ring_anchor.id else bond.a1_id
        other = canvas.model.get_atom(other_id)
        ring_neighbor_angles.append(
            angle_deg(anchor_pos, QPointF(other.x, other.y))
        )
    assert len(ring_neighbor_angles) == 2
    assert [
        angle_distance_deg(substituent_angle, angle)
        for angle in ring_neighbor_angles
    ] == pytest.approx([120.0, 120.0], abs=1e-6)

    grown_snapshot = _ring_snapshot(canvas)
    canvas.undo_stack.undo()
    app.processEvents()
    assert _ring_snapshot(canvas) == ring_snapshot
    canvas.undo_stack.redo()
    app.processEvents()
    assert _ring_snapshot(canvas) == grown_snapshot
    canvas.close()
    app.processEvents()


def test_explicit_drag_preserves_cursor_direction_and_preview_equals_commit() -> None:
    app = QApplication.instance() or QApplication([])
    canvas = ChemusonCanvas()
    _show_canvas(canvas)
    canvas.state.fixed_angles = True
    canvas.state.fixed_lengths = True
    canvas.state.angle_step_deg = 30
    anchor = canvas.model.add_atom("C", 400.0, 400.0, is_explicit=True)
    canvas._rebuild_items_from_model()
    canvas.set_current_tool("tool_bond")
    start = QPointF(anchor.x, anchor.y)
    cursor = endpoint_from_angle_len(start, 47.0, canvas.state.bond_length)
    start_pos = canvas.mapFromScene(start)
    cursor_pos = canvas.mapFromScene(cursor)
    assert canvas._pick_hover_target(canvas.mapToScene(start_pos))[0] == anchor.id
    QTest.mousePress(canvas.viewport(), Qt.MouseButton.LeftButton, pos=start_pos)
    QTest.mouseMove(canvas.viewport(), cursor_pos, 10)
    app.processEvents()

    preview_path = canvas._preview_bond_item.path()
    assert preview_path.elementCount() == 2
    preview_endpoint = QPointF(
        preview_path.elementAt(1).x,
        preview_path.elementAt(1).y,
    )
    preview_angle = angle_deg(start, preview_endpoint)
    cursor_scene = canvas.mapToScene(cursor_pos)
    cursor_angle = angle_deg(start, cursor_scene)
    assert angle_distance_deg(preview_angle, cursor_angle) < 1.5

    QTest.mouseRelease(canvas.viewport(), Qt.MouseButton.LeftButton, pos=cursor_pos)
    app.processEvents()
    assert len(canvas.model.atoms) == 2
    created = next(atom for atom_id, atom in canvas.model.atoms.items() if atom_id != anchor.id)
    assert (created.x, created.y) == pytest.approx(
        (preview_endpoint.x(), preview_endpoint.y()),
        abs=1e-7,
    )
    canvas.close()
    app.processEvents()


def test_1_phenylbutan_2_one_from_short_clicks_has_expected_graph_and_geometry(tmp_path) -> None:
    app = QApplication.instance() or QApplication([])
    canvas = ChemusonCanvas()
    ring_ids = _draw_free_ring(canvas)
    ring_anchor = max(
        (canvas.model.get_atom(atom_id) for atom_id in ring_ids),
        key=lambda atom: atom.x,
    )

    phenyl_methylene = _short_bond(canvas, ring_anchor.id)
    carbonyl_carbon = _short_bond(canvas, phenyl_methylene)
    oxygen = _short_bond(canvas, carbonyl_carbon, element="O", order=2)
    alpha_methylene = _short_bond(canvas, carbonyl_carbon)
    terminal_methyl = _short_bond(canvas, alpha_methylene)

    assert len(canvas.model.atoms) == 11
    assert len(canvas.model.bonds) == 11
    assert sum(atom.element == "C" for atom in canvas.model.atoms.values()) == 10
    assert sum(atom.element == "O" for atom in canvas.model.atoms.values()) == 1
    carbonyl_bond = canvas.model.find_bond_between(carbonyl_carbon, oxygen)
    assert carbonyl_bond is not None and carbonyl_bond.order == 2
    assert len([bond for bond in canvas.model.bonds.values() if bond.is_aromatic]) == 6
    assert canvas.model.find_bond_between(ring_anchor.id, phenyl_methylene) is not None
    assert canvas.model.find_bond_between(phenyl_methylene, carbonyl_carbon) is not None
    assert canvas.model.find_bond_between(carbonyl_carbon, alpha_methylene) is not None
    assert canvas.model.find_bond_between(alpha_methylene, terminal_methyl) is not None
    assert molgraph_to_smiles_isolated_or_error(
        canvas.model,
        timeout_s=2.0,
    ) == "CCC(=O)Cc1ccccc1"

    lengths = [
        math.hypot(
            canvas.model.get_atom(bond.a1_id).x - canvas.model.get_atom(bond.a2_id).x,
            canvas.model.get_atom(bond.a1_id).y - canvas.model.get_atom(bond.a2_id).y,
        )
        for bond in canvas.model.bonds.values()
    ]
    assert lengths == pytest.approx([canvas.state.bond_length] * 11, abs=1e-7)

    center_ids = (ring_anchor.id, phenyl_methylene, carbonyl_carbon, alpha_methylene)
    local_angles = []
    for atom_id in center_ids:
        atom = canvas.model.get_atom(atom_id)
        vectors = []
        for bond in canvas.model.bonds.values():
            if atom_id not in (bond.a1_id, bond.a2_id):
                continue
            other_id = bond.a2_id if bond.a1_id == atom_id else bond.a1_id
            other = canvas.model.get_atom(other_id)
            vectors.append(angle_deg(QPointF(atom.x, atom.y), QPointF(other.x, other.y)))
        local_angles.extend(
            angle_distance_deg(a, b)
            for index, a in enumerate(vectors)
            for b in vectors[index + 1 :]
        )
    assert local_angles == pytest.approx([120.0] * len(local_angles), abs=1e-6)
    assert _non_adjacent_crossing_count(canvas) == 0

    snapshot = _ring_snapshot(canvas)
    canvas.undo_stack.undo()
    app.processEvents()
    assert len(canvas.model.atoms) == 10
    canvas.undo_stack.redo()
    app.processEvents()
    assert _ring_snapshot(canvas) == snapshot

    document_path = tmp_path / "1-phenylbutan-2-one.cmsn"
    PersistenceManager.save_to_file(str(document_path), canvas)
    restored = ChemusonCanvas()
    PersistenceManager.load_from_file(str(document_path), restored)
    assert _ring_snapshot(restored) == snapshot
    restored.close()
    canvas.close()
    app.processEvents()
