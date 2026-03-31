"""Layout engine para diagramas electrónicos semánticos."""

from __future__ import annotations

from PyQt6.QtCore import QPointF, QRectF, Qt
from PyQt6.QtGui import QColor, QFont
from PyQt6.QtWidgets import QGraphicsTextItem

from chemuson.gui.diagram_models import DiagramLevel, SemanticDiagram
from chemuson.gui.items import BaseItem, EnergyDiagramItem, SemanticConnectorItem


ENERGY_SCALE = 48.0
LEVEL_HEIGHT = 30.0
LANE_TITLE_GAP = 18.0
TITLE_GAP = 10.0
DEFAULT_LEVEL_WIDTH = 96.0
SEMANTIC_TEXT_KIND_ROLE = 0x5F01
SEMANTIC_TEXT_ID_ROLE = 0x5F02


def _effective_level_width(level: DiagramLevel) -> float:
    """Compone un ancho visual compacto cuando el modelo usa el ancho default."""
    try:
        raw_width = float(level.width)
    except Exception:
        raw_width = DEFAULT_LEVEL_WIDTH
    if raw_width > 0.0 and abs(raw_width - DEFAULT_LEVEL_WIDTH) > 0.01:
        return max(28.0, raw_width)

    degeneracy = max(1, int(level.degeneracy))
    label_text = str(level.label or "").strip()
    slot_width = 18.0 if level.representation == "boxes" else 16.0
    slot_gap = 4.0
    slots_band = degeneracy * slot_width + max(0, degeneracy - 1) * slot_gap
    label_band = 0.0
    if label_text:
        label_band = min(74.0, max(18.0, len(label_text) * 6.2 + 8.0))
    return max(42.0, slots_band + label_band + 14.0)


def _label_side_for_x(x: float) -> str:
    if x > 80.0:
        return "right"
    return "left"


def _center_text_item(item: QGraphicsTextItem, center_x: float, top_y: float) -> None:
    rect = item.boundingRect()
    item.setPos(center_x - rect.width() * 0.5, top_y)


def _connector_points(source: EnergyDiagramItem, target: EnergyDiagramItem) -> tuple[QPointF, QPointF]:
    source_points = source.level_anchor_points()
    target_points = target.level_anchor_points()
    source_rect = source.sceneBoundingRect()
    target_rect = target.sceneBoundingRect()
    dx = target_rect.center().x() - source_rect.center().x()
    dy = target_rect.center().y() - source_rect.center().y()
    if abs(dx) >= abs(dy):
        if dx >= 0.0:
            return source_points["right_center"], target_points["left_center"]
        return source_points["left_center"], target_points["right_center"]
    if dy >= 0.0:
        return source_points["bottom_center"], target_points["top_center"]
    return source_points["top_center"], target_points["bottom_center"]


def _level_display_rect(level_item: EnergyDiagramItem) -> QRectF:
    return level_item.display_rect()


def build_items_from_semantic_diagram(diagram: SemanticDiagram) -> list[BaseItem]:
    """Convierte un SemanticDiagram a items Qt reutilizando EnergyDiagramItem."""
    lane_map = {lane.id: lane for lane in diagram.lanes}
    items: list[BaseItem] = []
    level_items: dict[str, EnergyDiagramItem] = {}

    for level in diagram.levels:
        lane = lane_map.get(level.lane_id)
        if lane is None:
            continue
        level_width = _effective_level_width(level)
        style_payload: dict[str, object] = {}
        if level.representation == "bar":
            style_payload = {
                "fill_visible": False,
                "box_stroke_visible": False,
                "box_base_visible": True,
            }
        item = EnergyDiagramItem(
            "custom_level",
            label=level.label,
            label_side=_label_side_for_x(float(lane.x)),
            occupancies=list(level.occupancies),
            slot_count=int(level.degeneracy),
            style_payload=style_payload,
            metadata={
                "semantic_level_id": level.id,
                "semantic_lane_id": level.lane_id,
                "semantic_representation": level.representation,
                **dict(level.metadata),
            },
            width=level_width,
            height=LEVEL_HEIGHT,
        )
        y_scene = -float(level.energy) * ENERGY_SCALE
        item.set_display_rect(
            QRectF(
                float(lane.x) - level_width * 0.5,
                y_scene - LEVEL_HEIGHT * 0.5,
                level_width,
                LEVEL_HEIGHT,
            )
        )
        level_items[level.id] = item
        items.append(item)

    connector_items: list[BaseItem] = []
    for connector in diagram.connectors:
        source_item = level_items.get(connector.source_level_id)
        target_item = level_items.get(connector.target_level_id)
        if source_item is None or target_item is None:
            continue
        start, end = _connector_points(source_item, target_item)
        connector_item = SemanticConnectorItem(start, end, style=connector.style)
        connector_items.append(connector_item)

    title_items: list[BaseItem] = []
    if diagram.title or diagram.lanes:
        content_top = min(
            (_level_display_rect(item).top() for item in level_items.values()),
            default=0.0,
        )
        font = QFont("Arial", 10)
        font.setBold(True)
        lane_font = QFont("Arial", 8)
        lane_font.setBold(True)
        show_lane_titles = bool(diagram.metadata.get("show_lane_titles", True))
        lane_title_top = content_top

        lane_items: list[tuple[QGraphicsTextItem, float]] = []
        if show_lane_titles:
            for lane in diagram.lanes:
                lane_title = str(lane.title or "").strip()
                if not lane_title:
                    continue
                lane_item = QGraphicsTextItem(lane_title)
                lane_item.setFont(lane_font)
                lane_item.setDefaultTextColor(QColor("#4A4A4A"))
                lane_item.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
                lane_item.setData(SEMANTIC_TEXT_KIND_ROLE, "lane_title")
                lane_item.setData(SEMANTIC_TEXT_ID_ROLE, str(lane.id))
                lane_items.append((lane_item, float(lane.x)))
            if lane_items:
                lane_title_top = (
                    content_top
                    - max(item.boundingRect().height() for item, _x in lane_items)
                    - LANE_TITLE_GAP
                )
                for lane_item, lane_x in lane_items:
                    _center_text_item(lane_item, lane_x, lane_title_top)
                    title_items.append(lane_item)

        if str(diagram.title or "").strip():
            title_item = QGraphicsTextItem(str(diagram.title))
            title_item.setFont(font)
            title_item.setDefaultTextColor(QColor("#202020"))
            title_item.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
            title_item.setData(SEMANTIC_TEXT_KIND_ROLE, "diagram_title")
            lane_centers = [float(lane.x) for lane in diagram.lanes] or [0.0]
            center_x = (min(lane_centers) + max(lane_centers)) * 0.5
            title_top = content_top - title_item.boundingRect().height() - TITLE_GAP
            if lane_items:
                title_top = lane_title_top - title_item.boundingRect().height() - TITLE_GAP
            _center_text_item(title_item, center_x, title_top)
            title_items.append(title_item)

    return [*connector_items, *items, *title_items]
