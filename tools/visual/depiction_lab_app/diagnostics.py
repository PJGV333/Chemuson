"""Optional scene-only geometry overlays for visual inspection."""
from __future__ import annotations

from PyQt6.QtCore import QPointF, QRectF, Qt
from PyQt6.QtGui import QBrush, QColor, QFont, QPen
from PyQt6.QtWidgets import (
    QGraphicsEllipseItem,
    QGraphicsItem,
    QGraphicsLineItem,
    QGraphicsPathItem,
    QGraphicsRectItem,
    QGraphicsSimpleTextItem,
)

from chemuson.core.model import BondStyle

from .renderer import ProductionCanvasAdapter

COLORS = {
    "atom": QColor("#2563EB"),
    "bond": QColor("#7C3AED"),
    "label": QColor("#DB2777"),
    "endpoint": QColor("#D97706"),
    "overlap": QColor("#DC2626"),
}


class GeometryDiagnostics:
    """Install/remove diagnostic marks without changing the chemical graph."""

    def __init__(self, adapter: ProductionCanvasAdapter) -> None:
        self.adapter = adapter
        self.items: list[QGraphicsItem] = []

    def clear(self) -> None:
        scene = self.adapter.canvas.scene
        for item in self.items:
            scene.removeItem(item)
        self.items.clear()
        self.adapter.diagnostic_items = []

    def _add(self, item: QGraphicsItem, z: float = 1000.0) -> None:
        item.setZValue(z)
        item.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
        item.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable, False)
        self.adapter.canvas.scene.addItem(item)
        self.items.append(item)

    def _mark_point(self, point: QPointF, color: QColor, radius: float = 2.6) -> None:
        marker = QGraphicsEllipseItem(
            point.x() - radius,
            point.y() - radius,
            radius * 2.0,
            radius * 2.0,
        )
        marker.setPen(QPen(color, 1.3))
        marker.setBrush(QBrush(Qt.BrushStyle.NoBrush))
        self._add(marker, z=1002.0)

    def _mark_rendered_ends(self, bond_item) -> None:
        path = bond_item.path()
        points = [
            QPointF(path.elementAt(index).x, path.elementAt(index).y)
            for index in range(path.elementCount())
        ]
        if len(points) < 2:
            return
        if bond_item.style == BondStyle.HASHED:
            first = (points[0] + points[1]) * 0.5
            last = (points[-2] + points[-1]) * 0.5
        elif bond_item.style == BondStyle.WEDGE:
            polygons = path.toFillPolygons()
            if not polygons or polygons[0].size() < 4:
                return
            polygon = polygons[0]
            first_tip = QPointF(polygon[0])
            last_tip = QPointF(polygon[polygon.size() - 2])
            tip_midpoint = (first_tip + last_tip) * 0.5
            atom1 = bond_item._atom1_ref
            atom2 = bond_item._atom2_ref
            ux = atom2.x - atom1.x
            uy = atom2.y - atom1.y
            length = max(1e-9, (ux * ux + uy * uy) ** 0.5)
            ux, uy = ux / length, uy / length
            candidates = [QPointF(polygon[index]) for index in range(1, polygon.size() - 1)]
            candidates.sort(
                key=lambda point: (point.x() - tip_midpoint.x()) * ux
                + (point.y() - tip_midpoint.y()) * uy,
                reverse=True,
            )
            base_midpoint = (
                (candidates[0] + candidates[1]) * 0.5
                if len(candidates) > 1
                else tip_midpoint
            )
            first, last = tip_midpoint, base_midpoint
        else:
            first, last = points[0], points[1]
        self._mark_point(bond_item.mapToScene(first), QColor("#059669"))
        self._mark_point(bond_item.mapToScene(last), QColor("#059669"))

    def show(self) -> None:
        self.clear()
        canvas = self.adapter.canvas
        for atom_id, atom in canvas.model.atoms.items():
            center = (float(atom.x), float(atom.y))
            center_mark = QGraphicsEllipseItem(center[0] - 3, center[1] - 3, 6, 6)
            center_mark.setPen(QPen(COLORS["atom"], 1.4))
            center_mark.setBrush(QBrush(Qt.BrushStyle.NoBrush))
            self._add(center_mark)
            atom_item = canvas.atom_items.get(atom_id)
            if atom_item is None:
                continue
            for text_item in (atom_item.label, atom_item.charge_label):
                if not text_item.isVisible() or not text_item.toPlainText().strip():
                    continue
                rect = text_item.mapRectToScene(text_item.boundingRect())
                label_box = QGraphicsRectItem(rect)
                label_box.setPen(QPen(COLORS["label"], 1.0, Qt.PenStyle.DashLine))
                label_box.setBrush(QBrush(Qt.BrushStyle.NoBrush))
                self._add(label_box)

        for bond_item in canvas.bond_items.values():
            self._mark_rendered_ends(bond_item)
            path = bond_item.mapToScene(bond_item.path())
            outline = QGraphicsPathItem(path)
            outline.setPen(QPen(COLORS["bond"], 1.0, Qt.PenStyle.DashLine))
            outline.setBrush(QBrush(Qt.BrushStyle.NoBrush))
            self._add(outline)
            bounds = path.boundingRect()
            box = QGraphicsRectItem(bounds)
            box.setPen(QPen(COLORS["endpoint"], 0.8, Qt.PenStyle.DotLine))
            box.setBrush(QBrush(Qt.BrushStyle.NoBrush))
            self._add(box, z=999.0)
            for atom_id in (bond_item.a1_id, bond_item.a2_id):
                atom = canvas.model.get_atom(atom_id)
                endpoint = QGraphicsEllipseItem(atom.x - 2.2, atom.y - 2.2, 4.4, 4.4)
                endpoint.setPen(QPen(COLORS["endpoint"], 1.1))
                endpoint.setBrush(QBrush(Qt.BrushStyle.NoBrush))
                self._add(endpoint, z=1001.0)

        note = QGraphicsSimpleTextItem(
            "Blue: atom centers  ·  orange: endpoints/bounds  ·  purple: renderer paths  ·  pink: label boxes  ·  red: geometric overlap"
        )
        note.setFont(QFont("DejaVu Sans", 9))
        note.setBrush(QBrush(QColor("#111827")))
        note.setPos(self.adapter.scene_rect.left() + 8, self.adapter.scene_rect.top() + 8)
        self._add(note, z=1003.0)
        ref_x = self.adapter.scene_rect.left() + 12.0
        ref_y = self.adapter.scene_rect.top() + 34.0
        reference = QGraphicsLineItem(
            ref_x,
            ref_y,
            ref_x + 40.0,
            ref_y,
        )
        reference.setPen(QPen(COLORS["atom"], 1.8))
        self._add(reference, z=1003.0)
        reference_label = QGraphicsSimpleTextItem("Lref = 40 scene units")
        reference_label.setFont(QFont("DejaVu Sans", 8))
        reference_label.setBrush(QBrush(COLORS["atom"]))
        reference_label.setPos(ref_x + 44.0, ref_y - 8.0)
        self._add(reference_label, z=1003.0)

        for atom_id, atom_item in canvas.atom_items.items():
            glyphs = [glyph for glyph in (atom_item.label, atom_item.charge_label) if glyph.isVisible()]
            connected = [
                canvas.bond_items[bond.id]
                for bond in canvas.model.bonds.values()
                if atom_id in {bond.a1_id, bond.a2_id} and bond.id in canvas.bond_items
            ]
            has_intersection = any(
                bond.mapToScene(bond.shape()).intersects(
                    glyph.mapToScene(glyph.shape())
                )
                for bond in connected
                for glyph in glyphs
            )
            if has_intersection:
                rect = atom_item.label.mapRectToScene(atom_item.label.boundingRect())
                marker = QGraphicsRectItem(QRectF(rect))
                marker.setPen(QPen(COLORS["overlap"], 2.0))
                marker.setBrush(QBrush(Qt.BrushStyle.NoBrush))
                self._add(marker, z=1002.0)

        self.adapter.diagnostic_items = list(self.items)

    def set_enabled(self, enabled: bool) -> None:
        if enabled:
            self.show()
        else:
            self.clear()
