from __future__ import annotations

import hashlib

from PyQt6.QtCore import QPointF
from PyQt6.QtGui import QImage, QPainterPath, QPolygonF
from PyQt6.QtWidgets import QApplication

from chemuson.gui.projections import (
    HaworthParameters,
    HaworthSubstituent,
    build_haworth_geometry,
)
from tools.visual.depiction_lab_app.renderer import ProductionCanvasAdapter


def _digest(image: QImage) -> str:
    image = image.convertToFormat(QImage.Format.Format_ARGB32)
    bits = image.bits()
    bits.setsize(image.sizeInBytes())
    return hashlib.sha256(bytes(bits)).hexdigest()


def _check_rendered_label_clearance(adapter: ProductionCanvasAdapter) -> None:
    canvas = adapter.canvas
    for atom_id, atom_item in canvas.atom_items.items():
        glyphs = [
            glyph
            for glyph in (atom_item.label, atom_item.charge_label)
            if glyph.isVisible() and glyph.toPlainText().strip()
        ]
        connected = [
            canvas.bond_items[bond.id]
            for bond in canvas.model.bonds.values()
            if atom_id in (bond.a1_id, bond.a2_id) and bond.id in canvas.bond_items
        ]
        for glyph in glyphs:
            glyph_shape = glyph.mapToScene(glyph.shape())
            assert not any(
                bond.mapToScene(bond.shape()).intersects(glyph_shape)
                for bond in connected
            ), f"bond intersects {glyph.toPlainText()} label at atom {atom_id}"


def _check_no_nonbonded_bond_crossings(adapter: ProductionCanvasAdapter) -> None:
    canvas = adapter.canvas
    bonds = list(canvas.model.bonds.values())
    shapes = {
        bond.id: canvas.bond_items[bond.id].mapToScene(canvas.bond_items[bond.id].shape())
        for bond in bonds
    }
    for index, bond in enumerate(bonds):
        for other in bonds[index + 1 :]:
            if {bond.a1_id, bond.a2_id}.isdisjoint({other.a1_id, other.a2_id}):
                assert not shapes[bond.id].intersects(shapes[other.id]), (
                    "nonbonded Haworth bonds cross",
                    bond.id,
                    other.id,
                )


def test_haworth_beta_d_glucose_uses_real_renderer_without_c3_label_collisions() -> None:
    app = QApplication.instance() or QApplication([])
    geometry = build_haworth_geometry(HaworthParameters())
    instance = geometry.to_projection_instance()
    adapter = ProductionCanvasAdapter.create(instance.render_graph)
    adapter.canvas.resize(720, 520)
    adapter.canvas.show()
    for _ in range(4):
        app.processEvents()

    try:
        semantics = geometry.semantics
        assert semantics is not None
        assert [(item.ring_position, item.orientation) for item in semantics.substituents] == [
            (1, "up"), (2, "down"), (3, "up"), (4, "down"), (5, "up")
        ]
        ring_path = QPainterPath()
        ring_polygon = QPolygonF(
            [QPointF(geometry.atoms[index].x, geometry.atoms[index].y) for index in semantics.ring_atom_indices]
        )
        ring_path.addPolygon(ring_polygon)
        c3 = next(item for item in semantics.substituents if item.ring_position == 3)
        c3_label = adapter.canvas.atom_items[c3.group_atom_index + 1].label
        label_center = c3_label.mapToScene(c3_label.boundingRect().center())
        assert not ring_path.contains(label_center)

        digests = set()
        for zoom, angle in ((0.75, 0.0), (1.0, 0.0), (1.5, 30.0), (2.0, 60.0)):
            adapter.set_view(zoom, angle)
            for _ in range(3):
                app.processEvents()
            _check_rendered_label_clearance(adapter)
            _check_no_nonbonded_bond_crossings(adapter)
            target = adapter.canvas.mapFromScene(adapter.scene_rect.center())
            viewport_center = adapter.canvas.viewport().rect().center()
            assert abs(target.x() - viewport_center.x()) <= 10
            assert abs(target.y() - viewport_center.y()) <= 10
            image = adapter.canvas.viewport().grab().toImage()
            assert image.width() >= 380 and image.height() >= 300
            dark_pixels = 0
            for y in range(image.height()):
                for x in range(image.width()):
                    color = image.pixelColor(x, y)
                    dark_pixels += max(color.red(), color.green(), color.blue()) < 170
            assert dark_pixels > 20
            digests.add(_digest(image))
        assert len(digests) == 4
    finally:
        adapter.close()
        app.processEvents()


def test_furanose_vertex_policy_renders_labels_clear_of_the_real_ring() -> None:
    app = QApplication.instance() or QApplication([])
    geometry = build_haworth_geometry(
        HaworthParameters(
            ring_size=5,
            sugar_family="generic",
            substituents=(
                HaworthSubstituent(1, "OH", "up"),
                HaworthSubstituent(2, "OH", "down"),
                HaworthSubstituent(3, "OH", "up"),
                HaworthSubstituent(4, "OH", "down"),
            ),
        )
    )
    adapter = ProductionCanvasAdapter.create(geometry.to_projection_instance().render_graph)
    adapter.canvas.resize(680, 500)
    adapter.canvas.show()
    for _ in range(4):
        app.processEvents()

    try:
        for zoom, angle in ((0.75, 0.0), (1.0, 25.0), (1.5, 0.0)):
            adapter.set_view(zoom, angle)
            app.processEvents()
            _check_rendered_label_clearance(adapter)
            _check_no_nonbonded_bond_crossings(adapter)
            image = adapter.canvas.viewport().grab().toImage()
            assert image.width() >= 380 and image.height() >= 300
            target = adapter.canvas.mapFromScene(adapter.scene_rect.center())
            center = adapter.canvas.viewport().rect().center()
            assert abs(target.x() - center.x()) <= 10
            assert abs(target.y() - center.y()) <= 10
    finally:
        adapter.close()
        app.processEvents()
