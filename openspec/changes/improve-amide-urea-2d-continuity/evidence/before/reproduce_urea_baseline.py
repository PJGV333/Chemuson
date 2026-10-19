from __future__ import annotations

import json
import math
import os
from pathlib import Path

from PyQt6.QtCore import QPointF, QRectF, Qt
from PyQt6.QtGui import QBrush, QColor, QPainter, QPainterPath, QPainterPathStroker, QTransform
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QApplication

from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.geom import angle_deg, angle_distance_deg, endpoint_from_angle_len
from chemuson.gui.theme import apply_theme

OUT = Path(os.environ.get("UREA_OUT", Path(__file__).parent))
OUT.mkdir(parents=True, exist_ok=True)
app = QApplication.instance() or QApplication([])
canvas = ChemusonCanvas()
canvas.resize(1000, 760)
canvas.show()
app.processEvents()
canvas.center_on_paper()
app.processEvents()
canvas.state.fixed_angles = True
canvas.state.fixed_lengths = True
canvas.state.angle_step_deg = 30
canvas.state.bond_length = 40.0


def scene_pos(atom_id: int) -> QPointF:
    atom = canvas.model.get_atom(atom_id)
    return QPointF(atom.x, atom.y)


def drag(start: QPointF, end: QPointF, *, element: str = "C", order: int = 1) -> int:
    canvas.set_current_tool("tool_bond")
    canvas.state.default_element = element
    canvas.state.active_bond_order = order
    canvas.state.active_bond_aromatic = False
    p0 = canvas.mapFromScene(start)
    p1 = canvas.mapFromScene(end)
    before = set(canvas.model.atoms)
    QTest.mousePress(canvas.viewport(), Qt.MouseButton.LeftButton, pos=p0)
    QTest.mouseMove(canvas.viewport(), p1, 10)
    app.processEvents()
    preview = canvas._preview_bond_item.path()
    preview_points = [
        [preview.elementAt(i).x, preview.elementAt(i).y]
        for i in range(preview.elementCount())
    ]
    QTest.mouseRelease(canvas.viewport(), Qt.MouseButton.LeftButton, pos=p1)
    app.processEvents()
    created = set(canvas.model.atoms) - before
    if len(created) != 1:
        raise RuntimeError(f"drag expected one atom, got {created}")
    return created.pop()


def short_click(anchor_id: int, *, element: str = "C", order: int = 1) -> tuple[int, dict]:
    canvas.set_current_tool("tool_bond")
    canvas.state.default_element = element
    canvas.state.active_bond_order = order
    canvas.state.active_bond_aromatic = False
    pos = canvas.mapFromScene(scene_pos(anchor_id))
    before = set(canvas.model.atoms)
    QTest.mousePress(canvas.viewport(), Qt.MouseButton.LeftButton, pos=pos)
    app.processEvents()
    preview = canvas._preview_bond_item.path()
    preview_end = [preview.elementAt(1).x, preview.elementAt(1).y] if preview.elementCount() > 1 else None
    QTest.mouseRelease(canvas.viewport(), Qt.MouseButton.LeftButton, pos=pos)
    app.processEvents()
    created = set(canvas.model.atoms) - before
    if len(created) != 1:
        raise RuntimeError(f"short click expected one atom at {anchor_id}, got {created}")
    new_id = created.pop()
    return new_id, {"preview_end": preview_end, "committed": [canvas.model.get_atom(new_id).x, canvas.model.get_atom(new_id).y]}


def glyph_positions(atom_id: int) -> dict:
    atom = canvas.model.get_atom(atom_id)
    item = canvas.atom_items[atom_id]
    text_item = item.label
    block = text_item.document().begin()
    text = block.text()
    layout = block.layout()
    doc_layout = text_item.document().documentLayout()
    block_rect = doc_layout.blockBoundingRect(block)
    result = {
        "text": text_item.toPlainText(),
        "atom_center": [atom.x, atom.y],
        "label_item_origin_scene": [text_item.mapToScene(QPointF(0.0, 0.0)).x(), text_item.mapToScene(QPointF(0.0, 0.0)).y()],
        "label_bounding_rect_local": [text_item.boundingRect().x(), text_item.boundingRect().y(), text_item.boundingRect().width(), text_item.boundingRect().height()],
        "glyphs": {},
    }
    cursor = text_item.textCursor()
    for index, char in enumerate(text):
        cursor.setPosition(index)
        cursor.setPosition(index + 1, cursor.MoveMode.KeepAnchor)
        fmt = cursor.charFormat()
        font = fmt.font()
        line = layout.lineForTextPosition(index)
        if not line.isValid():
            continue
        x1 = line.cursorToX(index)
        x2 = line.cursorToX(index + 1)
        if isinstance(x1, tuple):
            x1 = x1[0]
        if isinstance(x2, tuple):
            x2 = x2[0]
        baseline = block_rect.y() + line.position().y() + line.ascent()
        glyph = QPainterPath()
        glyph.addText(block_rect.x() + line.position().x() + float(x1), baseline, font, char)
        glyph_scene = text_item.mapToScene(glyph)
        bounds = glyph_scene.boundingRect()
        result["glyphs"].setdefault(char, []).append({
            "layout_cursor_x": [float(x1), float(x2)],
            "font_point_size": font.pointSizeF(),
            "scene_bbox": [bounds.x(), bounds.y(), bounds.width(), bounds.height()],
            "scene_center": [bounds.center().x(), bounds.center().y()],
        })
    return result


# Aromatic ring via a real ring-tool press/move/release gesture.
canvas.set_current_tool("tool_ring")
canvas.state.active_ring_size = 6
canvas.state.active_ring_aromatic = True
center = canvas.mapToScene(canvas.viewport().rect().center())
ring_end = endpoint_from_angle_len(center, 0.0, 40.0)
ring_p0 = canvas.mapFromScene(center)
ring_p1 = canvas.mapFromScene(ring_end)
QTest.mousePress(canvas.viewport(), Qt.MouseButton.LeftButton, pos=ring_p0)
QTest.mouseMove(canvas.viewport(), ring_p1, 10)
QTest.mouseRelease(canvas.viewport(), Qt.MouseButton.LeftButton, pos=ring_p1)
app.processEvents()
ring_ids = sorted(canvas.model.atoms)
ring_anchor = max((canvas.model.get_atom(i) for i in ring_ids), key=lambda atom: atom.x)

# Ar–NH–C(=O)–NH–CH3: short clicks for local-assist, explicit drags for O and methyl.
n1, gesture_n1 = short_click(ring_anchor.id, element="N")
carbonyl, gesture_c = short_click(n1, element="C")
co_angle = angle_deg(scene_pos(carbonyl), scene_pos(n1)) + 120.0
oxygen = drag(scene_pos(carbonyl), endpoint_from_angle_len(scene_pos(carbonyl), co_angle, 40.0), element="O", order=2)
n2, gesture_n2 = short_click(carbonyl, element="N")
n2_to_c_angle = angle_deg(scene_pos(n2), scene_pos(carbonyl)) - 120.0
methyl = drag(scene_pos(n2), endpoint_from_angle_len(scene_pos(n2), n2_to_c_angle, 40.0), element="C")
canvas.state.default_element = "C"
app.processEvents()

atoms = list(canvas.model.atoms.values())
bonds = list(canvas.model.bonds.values())
path_ids = [ring_anchor.id, n1, carbonyl, oxygen, n2, methyl]
path_pairs = list(zip(path_ids, path_ids[1:]))


def label_ray_hit(atom_id: int, ux: float, uy: float, width: float) -> dict:
    atom = canvas.model.get_atom(atom_id)
    item = canvas.atom_items.get(atom_id)
    if item is None or not item.label.isVisible():
        return {"intersects": False, "max_forward_extent": None}
    text_path = item.label.mapToScene(item.label.shape())
    ray = QPainterPath()
    ray.moveTo(atom.x, atom.y)
    ray.lineTo(atom.x + ux * canvas.state.bond_length, atom.y + uy * canvas.state.bond_length)
    stroker = QPainterPathStroker()
    stroker.setWidth(width)
    corridor = stroker.createStroke(ray)
    hit = text_path.intersected(corridor)
    projections = [
        (point.x() - atom.x) * ux + (point.y() - atom.y) * uy
        for polygon in hit.toSubpathPolygons()
        for point in polygon
    ]
    return {"intersects": not hit.isEmpty(), "max_forward_extent": max(projections) if projections else None}


def expanded_label_ray_hit(atom_id: int, ux: float, uy: float, radius: float) -> dict:
    atom = canvas.model.get_atom(atom_id)
    item = canvas.atom_items.get(atom_id)
    if item is None or not item.label.isVisible():
        return {"intersects": False, "max_centerline_extent": None}
    text_path = item.label.mapToScene(item.label.shape())
    expander = QPainterPathStroker()
    expander.setWidth(max(0.0, 2.0 * radius))
    expander.setCapStyle(Qt.PenCapStyle.RoundCap)
    expander.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    expanded = text_path.united(expander.createStroke(text_path))
    ray = QPainterPath()
    ray.moveTo(atom.x, atom.y)
    ray.lineTo(atom.x + ux * canvas.state.bond_length, atom.y + uy * canvas.state.bond_length)
    ray_stroker = QPainterPathStroker()
    ray_stroker.setWidth(0.01)
    hit = expanded.intersected(ray_stroker.createStroke(ray))
    projections = [
        (point.x() - atom.x) * ux + (point.y() - atom.y) * uy
        for polygon in hit.toSubpathPolygons()
        for point in polygon
    ]
    return {"intersects": not hit.isEmpty(), "max_centerline_extent": max(projections) if projections else None}


def bond_geometry(bond) -> dict:
    a1 = canvas.model.get_atom(bond.a1_id)
    a2 = canvas.model.get_atom(bond.a2_id)
    item = canvas.bond_items[bond.id]
    dx, dy = a2.x - a1.x, a2.y - a1.y
    stroke_width = float(item.pen().widthF())
    clearance = max(1.0, stroke_width * 0.5 + 1.0)
    center_length = math.hypot(dx, dy)
    ux, uy = (dx / center_length, dy / center_length) if center_length > 1e-9 else (0.0, 0.0)
    subpaths = []
    for polygon in item.path().toSubpathPolygons():
        points = [item.mapToScene(point) for point in polygon]
        visible_length = sum(math.hypot(b.x() - a.x(), b.y() - a.y()) for a, b in zip(points, points[1:]))
        subpaths.append({
            "points": [[point.x(), point.y()] for point in points],
            "visible_path_length": visible_length,
        })
    return {
        "id": bond.id,
        "a1": bond.a1_id,
        "a2": bond.a2_id,
        "elements": [a1.element, a2.element],
        "order": bond.order,
        "aromatic": bond.is_aromatic,
        "center_distance": center_length,
        "reference_length": canvas.state.bond_length,
        "direction_deg": angle_deg(QPointF(a1.x, a1.y), QPointF(a2.x, a2.y)),
        "render_stroke_width": stroke_width,
        "nominal_clearance": clearance,
        "trim_start_label": item._label_shrink_start,
        "trim_end_label": item._label_shrink_end,
        "trim_start_endpoint": item._endpoint_trim_start,
        "trim_end_endpoint": item._endpoint_trim_end,
        "trim_start_total": item._label_shrink_start + item._endpoint_trim_start,
        "trim_end_total": item._label_shrink_end + item._endpoint_trim_end,
        "main_path_subpaths": subpaths,
        "secondary_path_subpaths": [
            {
                "points": [[p.x(), p.y()] for p in [item.mapToScene(poly[0]), *[item.mapToScene(poly[i]) for i in range(1, len(poly))]]],
                "visible_path_length": sum(math.hypot(poly[i].x() - poly[i-1].x(), poly[i].y() - poly[i-1].y()) for i in range(1, len(poly))),
            }
            for poly in item._secondary_path.toSubpathPolygons()
        ],
        "trim_projections": {
            "start_label_extent": canvas._label_projection_extent_for_atom(bond.a1_id, ux, uy),
            "end_label_extent": canvas._label_projection_extent_for_atom(bond.a2_id, -ux, -uy),
            "start_ray_intersection_stroke_1px": label_ray_hit(bond.a1_id, ux, uy, 1.0),
            "end_ray_intersection_stroke_1px": label_ray_hit(bond.a2_id, -ux, -uy, 1.0),
            "start_ray_intersection_actual_stroke": label_ray_hit(bond.a1_id, ux, uy, stroke_width),
            "end_ray_intersection_actual_stroke": label_ray_hit(bond.a2_id, -ux, -uy, stroke_width),
            "start_expanded_label_centerline_extent": expanded_label_ray_hit(bond.a1_id, ux, uy, clearance),
            "end_expanded_label_centerline_extent": expanded_label_ray_hit(bond.a2_id, -ux, -uy, clearance),
        },
    }

pair_angles = {}
for atom_id in path_ids:
    atom = canvas.model.get_atom(atom_id)
    neighbor_angles = []
    for bond in bonds:
        if atom_id not in (bond.a1_id, bond.a2_id):
            continue
        other_id = bond.a2_id if bond.a1_id == atom_id else bond.a1_id
        other = canvas.model.get_atom(other_id)
        neighbor_angles.append((other_id, angle_deg(QPointF(atom.x, atom.y), QPointF(other.x, other.y))))
    pair_angles[str(atom_id)] = [
        {"neighbors": [a[0], b[0]], "angle_deg": angle_distance_deg(a[1], b[1])}
        for index, a in enumerate(neighbor_angles)
        for b in neighbor_angles[index + 1 :]
    ]

report = {
    "structure": "Ar-NH-C(=O)-NH-CH3",
    "path_atom_ids": path_ids,
    "gesture_events": {
        "ring": "QTest press/move/release ring tool",
        "Ar-N": {"type": "short click", **gesture_n1},
        "N-C": {"type": "short click", **gesture_c},
        "C-O": {"type": "directed drag", "atom_id": oxygen},
        "C-N": {"type": "short click", **gesture_n2},
        "N-CH3": {"type": "directed drag", "atom_id": methyl},
    },
    "atoms": [
        {
            "id": atom.id,
            "element": atom.element,
            "x": atom.x,
            "y": atom.y,
            "implicit_h": canvas.model.implicit_h_count(atom.id),
            "degree": canvas.model.atom_degree(atom.id),
            "label": glyph_positions(atom.id) if atom.id in {n1, n2, oxygen} else None,
        }
        for atom in atoms
    ],
    "bonds": [bond_geometry(bond) for bond in bonds],
    "reference_path_bonds": [
        bond_geometry(canvas.model.find_bond_between(a, b))
        for a, b in path_pairs
        if canvas.model.find_bond_between(a, b) is not None
    ],
    "angles": pair_angles,
}
# Capture the current production canvas in both app themes and zoom levels.
canvas._update_hover(QPointF(-100000.0, -100000.0))
for theme, background in (("light", "#F1F5F9"), ("dark", "#111827")):
    apply_theme(canvas, theme)
    canvas.scene.setBackgroundBrush(QBrush(QColor(background)))
    for scale in (1.0, 1.25):
        canvas.setTransform(QTransform().scale(scale, scale))
        canvas.center_on_paper()
        app.processEvents()
        image = canvas.viewport().grab().toImage()
        if not image.save(str(OUT / f"urea-baseline-{theme}-{scale:g}x.png")):
            raise RuntimeError(f"failed to save {theme} {scale}x capture")
(OUT / "urea-baseline-metrics.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
print(json.dumps(report, indent=2))
canvas.close()
app.processEvents()
