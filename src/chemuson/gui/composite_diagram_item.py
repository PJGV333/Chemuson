"""Item compuesto para diagramas electrónicos semánticos."""

from __future__ import annotations

from PyQt6.QtCore import QPointF, QRectF, QSizeF, Qt
from PyQt6.QtGui import QFont, QPainter, QPainterPath
from PyQt6.QtWidgets import QGraphicsItem, QGraphicsSceneMouseEvent, QGraphicsTextItem, QStyleOptionGraphicsItem, QWidget

from chemuson.gui.diagram_layout import (
    SEMANTIC_TEXT_ID_ROLE,
    SEMANTIC_TEXT_KIND_ROLE,
    build_items_from_semantic_diagram,
)
from chemuson.gui.diagram_models import SemanticDiagram
from chemuson.gui.energy_diagrams import refresh_semantic_diagram_metadata
from chemuson.gui.items import BaseItem, EnergyDiagramItem, SemanticConnectorItem


class CompositeDiagramItem(BaseItem):
    """Agrupa un diagrama semántico completo como una sola entidad seleccionable."""

    def __init__(self, semantic_diagram: SemanticDiagram) -> None:
        super().__init__()
        self.semantic_diagram = semantic_diagram
        self._children: list[BaseItem] = []
        self._level_items: dict[str, EnergyDiagramItem] = {}
        self._display_size = QSizeF(1.0, 1.0)

        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsFocusable)
        self.setAcceptedMouseButtons(Qt.MouseButton.LeftButton | Qt.MouseButton.RightButton)
        self.setZValue(44)
        self._build_children()
        self.setTransformOriginPoint(self.boundingRect().center())

    def _build_children(self) -> None:
        raw_items = build_items_from_semantic_diagram(self.semantic_diagram)
        if not raw_items:
            self._display_size = QSizeF(1.0, 1.0)
            return

        bbox: QRectF | None = None
        for item in raw_items:
            candidate = item.sceneBoundingRect()
            if not candidate.isValid() or candidate.isNull():
                continue
            bbox = candidate if bbox is None else bbox.united(candidate)
        if bbox is None:
            self._display_size = QSizeF(1.0, 1.0)
            return

        for item in raw_items:
            self._shift_item_to_local_origin(item, bbox.topLeft())
            item.setParentItem(self)
            item.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable, False)
            item.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsFocusable, isinstance(item, EnergyDiagramItem))
            if hasattr(item, "setAcceptedMouseButtons"):
                item.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
            if isinstance(item, EnergyDiagramItem):
                level_id = str(item.metadata().get("semantic_level_id", "") or "")
                if level_id:
                    self._level_items[level_id] = item
            self._children.append(item)

        self._display_size = QSizeF(max(1.0, bbox.width()), max(1.0, bbox.height()))

    @staticmethod
    def _shift_item_to_local_origin(item: BaseItem, offset: QPointF) -> None:
        if isinstance(item, EnergyDiagramItem):
            rect = item.display_rect()
            item.set_display_rect(
                QRectF(
                    rect.x() - offset.x(),
                    rect.y() - offset.y(),
                    rect.width(),
                    rect.height(),
                )
            )
            return
        if isinstance(item, SemanticConnectorItem):
            start, end = item.points()
            item.set_points(
                QPointF(start.x() - offset.x(), start.y() - offset.y()),
                QPointF(end.x() - offset.x(), end.y() - offset.y()),
            )
            return
        item.setPos(item.pos() - offset)

    def boundingRect(self) -> QRectF:
        return QRectF(0.0, 0.0, self._display_size.width(), self._display_size.height())

    def shape(self) -> QPainterPath:
        path = QPainterPath()
        path.addRect(self.boundingRect())
        return path

    def contains(self, point: QPointF) -> bool:
        return self.boundingRect().contains(point)

    def paint(
        self,
        painter: QPainter,
        option: QStyleOptionGraphicsItem,
        widget: QWidget | None = None,
    ) -> None:
        del painter, option, widget

    def child_items(self) -> list[BaseItem]:
        return list(self._children)

    def level_items(self) -> list[EnergyDiagramItem]:
        return list(self._level_items.values())

    def _editor_view(self):
        scene = self.scene()
        if scene is None:
            return None
        for view in scene.views():
            if hasattr(view, "_prompt_semantic_diagram_title"):
                return view
        return None

    def _text_item_at(self, point: QPointF) -> QGraphicsTextItem | None:
        for item in reversed(self._children):
            if not isinstance(item, QGraphicsTextItem):
                continue
            local_point = item.mapFromParent(point)
            try:
                if item.contains(local_point) or item.boundingRect().contains(local_point):
                    return item
            except RuntimeError:
                continue
        return None

    def edit_target_at(self, point: QPointF) -> dict[str, object] | None:
        text_item = self._text_item_at(point)
        if text_item is not None:
            kind = str(text_item.data(SEMANTIC_TEXT_KIND_ROLE) or "").strip()
            if kind == "diagram_title":
                return {"kind": "diagram_title"}
            if kind == "lane_title":
                return {
                    "kind": "lane_title",
                    "lane_id": str(text_item.data(SEMANTIC_TEXT_ID_ROLE) or ""),
                }

        level_item = self._level_item_at(point)
        if level_item is None:
            return None
        local_point = level_item.mapFromParent(point)
        slot_index = level_item._slot_index_at(local_point)
        level_id = str(level_item.metadata().get("semantic_level_id", "") or "")
        if slot_index is not None:
            return {
                "kind": "occupancy",
                "level_id": level_id,
                "slot_index": int(slot_index),
            }
        if level_id:
            return {"kind": "level_label", "level_id": level_id}
        return None

    def _level_item_at(self, point: QPointF) -> EnergyDiagramItem | None:
        for item in reversed(self._children):
            if not isinstance(item, EnergyDiagramItem):
                continue
            local_point = item.mapFromParent(point)
            try:
                if item.shape().contains(local_point) or item.boundingRect().contains(local_point):
                    return item
            except RuntimeError:
                continue
        return None

    def mouseDoubleClickEvent(self, event: QGraphicsSceneMouseEvent) -> None:
        target = self.edit_target_at(event.pos())
        if target is None:
            super().mouseDoubleClickEvent(event)
            return
        self.setSelected(True)
        view = self._editor_view()
        kind = str(target.get("kind", "") or "")
        if kind == "occupancy":
            level_item = self._level_items.get(str(target.get("level_id", "") or ""))
            if level_item is None:
                super().mouseDoubleClickEvent(event)
                return
            level_item.begin_direct_edit(int(target.get("slot_index", 0) or 0))
            event.accept()
            return
        if view is None:
            super().mouseDoubleClickEvent(event)
            return
        if kind == "diagram_title":
            view._prompt_semantic_diagram_title(self)
            event.accept()
            return
        if kind == "lane_title":
            view._prompt_semantic_diagram_lane_title(
                self,
                str(target.get("lane_id", "") or ""),
            )
            event.accept()
            return
        if kind == "level_label":
            view._prompt_semantic_diagram_level_label(
                self,
                str(target.get("level_id", "") or ""),
            )
            event.accept()
            return
        super().mouseDoubleClickEvent(event)

    def _sync_semantic_from_children(self) -> None:
        occupancy_map = {
            "empty": 0,
            "up": 1,
            "down": 1,
            "pair": 2,
            "upup": 2,
            "downdown": 2,
        }
        for level in self.semantic_diagram.levels:
            item = self._level_items.get(level.id)
            if item is None:
                continue
            level.label = item.label()
            level.occupancies = [
                occupancy_map.get(value, 0) for value in item.occupancies()
            ]
        refresh_semantic_diagram_metadata(self.semantic_diagram)

    def _clear_children(self) -> None:
        for item in list(self._children):
            try:
                if item.parentItem() is self:
                    item.setParentItem(None)
                scene = item.scene()
                if scene is not None:
                    scene.removeItem(item)
            except RuntimeError:
                continue
        self._children.clear()
        self._level_items.clear()

    def rebuild_from_semantic_diagram(self, diagram: SemanticDiagram) -> None:
        """Reconstruye el item desde un nuevo modelo preservando su presentación."""
        rect = self.display_rect()
        rotation = self.rotation()
        z_value = self.zValue()
        opacity = self.opacity()
        selected = self.isSelected()
        self.semantic_diagram = SemanticDiagram.from_json_dict(diagram.to_json_dict())
        refresh_semantic_diagram_metadata(self.semantic_diagram)
        self.prepareGeometryChange()
        self._clear_children()
        self._build_children()
        self.set_display_rect(rect)
        self.setRotation(rotation)
        self.setZValue(z_value)
        self.setOpacity(opacity)
        self.setTransformOriginPoint(self.boundingRect().center())
        if selected:
            self.setSelected(True)
        self.update()

    def set_diagram_title(self, text: str) -> bool:
        normalized = str(text or "")
        if normalized == self.semantic_diagram.title:
            return False
        self._sync_semantic_from_children()
        self.semantic_diagram.title = normalized
        self.rebuild_from_semantic_diagram(self.semantic_diagram)
        return True

    def set_lane_title(self, lane_id: str, text: str) -> bool:
        target_lane_id = str(lane_id or "")
        normalized = str(text or "")
        for lane in self.semantic_diagram.lanes:
            if lane.id != target_lane_id:
                continue
            if normalized == lane.title:
                return False
            self._sync_semantic_from_children()
            lane.title = normalized
            self.rebuild_from_semantic_diagram(self.semantic_diagram)
            return True
        return False

    def set_level_label(self, level_id: str, text: str) -> bool:
        target_level_id = str(level_id or "")
        normalized = str(text or "")
        for level in self.semantic_diagram.levels:
            if level.id != target_level_id:
                continue
            if normalized == level.label:
                return False
            level.label = normalized
            level_item = self._level_items.get(target_level_id)
            if level_item is not None:
                level_item.set_label(normalized)
            return True
        return False

    def display_rect(self) -> QRectF:
        return QRectF(self.pos().x(), self.pos().y(), self._display_size.width(), self._display_size.height())

    def set_display_rect(self, rect: QRectF) -> None:
        normalized = QRectF(rect).normalized()
        old_width = max(1.0, self._display_size.width())
        old_height = max(1.0, self._display_size.height())
        scale_x = max(0.01, normalized.width() / old_width)
        scale_y = max(0.01, normalized.height() / old_height)

        for item in self._children:
            if isinstance(item, EnergyDiagramItem):
                child_rect = item.display_rect()
                item.set_display_rect(
                    QRectF(
                        child_rect.x() * scale_x,
                        child_rect.y() * scale_y,
                        child_rect.width() * scale_x,
                        child_rect.height() * scale_y,
                    )
                )
            elif isinstance(item, SemanticConnectorItem):
                start, end = item.points()
                item.set_points(
                    QPointF(start.x() * scale_x, start.y() * scale_y),
                    QPointF(end.x() * scale_x, end.y() * scale_y),
                )
            elif isinstance(item, QGraphicsTextItem):
                item.setPos(item.pos().x() * scale_x, item.pos().y() * scale_y)
                font = QFont(item.font())
                size = font.pointSizeF()
                if size <= 0.0:
                    size = float(font.pointSize()) if font.pointSize() > 0 else 10.0
                font.setPointSizeF(max(1.0, size * min(scale_x, scale_y)))
                item.setFont(font)
            else:
                item.setPos(item.pos().x() * scale_x, item.pos().y() * scale_y)

        self.prepareGeometryChange()
        self.setPos(normalized.topLeft())
        self._display_size = QSizeF(max(1.0, normalized.width()), max(1.0, normalized.height()))
        self.setTransformOriginPoint(self.boundingRect().center())
        self.update()

    def to_json(self) -> dict:
        self._sync_semantic_from_children()
        rect = self.display_rect()
        return {
            "semantic_diagram": self.semantic_diagram.to_json_dict(),
            "x": rect.x(),
            "y": rect.y(),
            "width": rect.width(),
            "height": rect.height(),
            "rotation": self.rotation(),
            "z": self.zValue(),
            "opacity": self.opacity(),
        }

    @classmethod
    def from_json(cls, payload: dict) -> "CompositeDiagramItem":
        item = cls(SemanticDiagram.from_json_dict(payload.get("semantic_diagram", {})))
        item.set_display_rect(
            QRectF(
                float(payload.get("x", 0.0)),
                float(payload.get("y", 0.0)),
                float(payload.get("width", 1.0)),
                float(payload.get("height", 1.0)),
            )
        )
        item.setRotation(float(payload.get("rotation", 0.0)))
        item.setZValue(float(payload.get("z", 44.0)))
        raw_opacity = payload.get("opacity", 1.0)
        try:
            item.setOpacity(1.0 if raw_opacity is None else float(raw_opacity))
        except Exception:
            item.setOpacity(1.0)
        return item
