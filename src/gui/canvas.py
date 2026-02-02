"""
Chemuson Canvas
Page-based canvas using QGraphicsView/QGraphicsScene for document-style editing.
"""
from __future__ import annotations

import json
import math
from typing import Dict, Iterable, List, Optional, Tuple

from PyQt6.QtWidgets import (
    QApplication,
    QMenu,
    QGraphicsView,
    QGraphicsScene,
    QGraphicsRectItem,
    QGraphicsPathItem,
    QGraphicsDropShadowEffect,
    QGraphicsPixmapItem,
    QGraphicsEllipseItem,
    QGraphicsLineItem,
    QGraphicsTextItem,
    QGraphicsItem,
    QInputDialog,
)
from PyQt6.QtGui import (
    QPainter,
    QPen,
    QColor,
    QBrush,
    QWheelEvent,
    QUndoStack,
    QImage,
    QPixmap,
    QPainterPath,
    QFont,
    QTextCharFormat,
)
from PyQt6.QtCore import Qt, QPointF, QRectF, QRect, QSize, QBuffer, QMimeData, pyqtSignal

try:
    from PyQt6.QtSvg import QSvgGenerator
except Exception:  # Optional Qt module at runtime
    QSvgGenerator = None

from core.model import Bond, BondStyle, BondStereo, ChemState, MolGraph, VALENCE_MAP
from gui.items import (
    AtomItem,
    BondItem,
    AromaticCircleItem,
    ArrowItem,
    PreviewArrowItem,
    BracketItem,
    HoverAtomIndicatorItem,
    HoverBondIndicatorItem,
    OptimizeZoneItem,
    PreviewBondItem,
    PreviewRingItem,
    PreviewChainItem,
    PreviewChainLabelItem,
    PreviewChainItem,
    PreviewChainLabelItem,
    TextAnnotationItem,
    ABBREVIATION_LABELS,
)
from gui.style import CHEMDOODLE_LIKE, DrawingStyle
from gui.geom import (
    angle_deg,
    snap_angle_deg,
    endpoint_from_angle_len,
    choose_optimal_direction,
    geometry_for_bond,
    SP3_BOND_ANGLE_DEG,
    candidate_directions_deg,
    filter_occupied_angles_deg,
    pick_closest_direction_deg,
    angle_distance_deg,
    segments_intersect,
    segment_min_distance,
    segments_nearly_equal,
    closest_atom,
    closest_bond,
    bond_side,
)
from gui.commands import (
    AddAtomCommand,
    AddBondCommand,
    AddRingCommand,
    AddChainCommand,
    AddArrowCommand,
    AddBracketCommand,
    ChangeAtomCommand,
    ChangeBondCommand,
    ChangeBondLengthCommand,
    ChangeChargeCommand,
    DeleteSelectionCommand,
    MoveAtomsCommand,
    MoveTextItemsCommand,
)
from gui.dialogs import AtomLabelDialog
from chemio.rdkit_io import (
    molgraph_to_molfile,
    molgraph_to_smiles,
    molfile_to_molgraph,
    smiles_to_molgraph,
    kekulize_display_orders,
)


# Paper dimensions (default: A4-ish, in px)
DEFAULT_PAPER_WIDTH = 800
DEFAULT_PAPER_HEIGHT = 1000
PAPER_MARGIN = 40  # Margins where atoms shouldn't be placed
ATOM_HIT_RADIUS = 16.0
DEFAULT_BOND_LENGTH = 40.0
HOVER_ATOM_RADIUS = 16.0
HOVER_BOND_DISTANCE = 10.0
OPTIMIZE_ZONE_SCALE = 1.2
CHAIN_MAX_BONDS = 12
ANGLE_OCCUPIED_TOLERANCE_DEG = 20.0
MIN_ATOM_DIST_SCALE = 0.65
MIN_BOND_DIST_SCALE = 0.2
COLLISION_LENGTH_BOOST = 1.2
BOND_OVERLAP_TOLERANCE_PX = 5.0
BOND_DRAG_THRESHOLD_PX = 6.0
BOND_LAST_ANGLE_TOLERANCE_DEG = 20.0
CLIPBOARD_RENDER_SCALE = 5.0
FUNCTIONAL_GROUP_LABELS = [
    "NH2",
    "NO2",
    "OH",
    "COOH",
    "CO2H",
    "CHO",
    "CN",
    "SO3H",
    "SH",
    "OMe",
    "OEt",
    "Me",
    "Et",
    "iPr",
    "tBu",
    "TBS",
    "Si",
    "Ph",
    "R",
]
FUNCTIONAL_GROUP_ALIASES = {label.lower(): label for label in FUNCTIONAL_GROUP_LABELS}
ELEMENT_SYMBOLS = {"C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "H", "Si", "B", "Se", "Li", "Na", "K", "Mg"}
IMPLICIT_H_ELEMENTS = {"C", "N", "O", "S", "P"}
LABEL_OFFSET_SCALE = 0.6
LABEL_OFFSET_MIN_PX = 4.0
RULER_THICKNESS_PX = 24
RULER_MAJOR_STEP_PX = 100
RULER_MINOR_STEP_PX = 20
GRID_MINOR_STEP_PX = 20
GRID_MAJOR_STEP_PX = 100
SELECTION_BOX_PADDING_PX = 6.0
SELECTION_HANDLE_RADIUS_PX = 6.0
SELECTION_HANDLE_OFFSET_PX = 18.0
GRID_MINOR_STEP_PX = 20
GRID_MAJOR_STEP_PX = 100
SELECTION_BOX_PADDING_PX = 6.0
SELECTION_HANDLE_RADIUS_PX = 6.0
SELECTION_HANDLE_OFFSET_PX = 18.0


class ChemusonCanvas(QGraphicsView):
    """
    Page-based canvas for drawing molecules.
    Uses QGraphicsView/QGraphicsScene with a centered paper sheet.
    """
    
    selection_changed = pyqtSignal(int, int, int, dict)

    def __init__(self, parent=None) -> None:
        super().__init__(parent)

        self.scene = QGraphicsScene()
        self.setScene(self.scene)

        self.paper_width = DEFAULT_PAPER_WIDTH
        self.paper_height = DEFAULT_PAPER_HEIGHT
        self.model = MolGraph()
        self.state = ChemState()
        self.undo_stack = QUndoStack(self)
        self.drawing_style: DrawingStyle = CHEMDOODLE_LIKE
        self._ring_centers: dict[int, QPointF] = {}
        self._next_ring_id = 1

        self.atom_items: dict[int, AtomItem] = {}
        self.bond_items: dict[int, BondItem] = {}
        self.aromatic_circles: list[AromaticCircleItem] = []

        self.current_tool: str = self.state.active_tool
        self.bond_anchor_id: Optional[int] = None
        self.hovered_atom_id: Optional[int] = None
        self.hovered_bond_id: Optional[int] = None
        self._last_scene_pos = QPointF(0, 0)
        self._bond_last_angle: Optional[float] = None
        self._bond_zigzag_sign = 1
        self._drag_free_orientation = False
        self._bond_drag_start: Optional[QPointF] = None
        self._arrow_start_pos: Optional[QPointF] = None
        self.arrow_items: list[ArrowItem] = []
        self._bracket_drag_start: Optional[QPointF] = None
        self._bracket_preview: Optional[QGraphicsRectItem] = None
        self.bracket_items: list[BracketItem] = []
        
        # Text tool state
        self._current_text_settings = {
            "family": "Arial",
            "size": 12,
            "bold": False,
            "italic": False,
            "underline": False,
            "sub": False,
            "sup": False,
            "color": QColor("black")
        }

        self._dragging_selection = False
        self._drag_start_pos: Optional[QPointF] = None
        self._drag_start_positions: Dict[int, Tuple[float, float]] = {}
        self._drag_has_moved = False
        self._select_drag_mode: Optional[str] = None
        self._select_start_pos: Optional[QPointF] = None
        self._select_path: Optional[QPainterPath] = None
        self._select_preview_path: Optional[QGraphicsPathItem] = None
        self._select_preview_rect: Optional[QGraphicsRectItem] = None
        self._select_additive = False

        self._drag_mode = "none"
        self._drag_anchor: Optional[dict] = None
        self._ring_last_vertices: Optional[List[QPointF]] = None
        self._chain_last_points: Optional[List[QPointF]] = None
        self._overlays_ready = False
        self.show_rulers = False
        self.show_grid = False
        self._grid_minor_item: Optional[QGraphicsPathItem] = None
        self._grid_major_item: Optional[QGraphicsPathItem] = None
        self._selection_box: Optional[QGraphicsRectItem] = None
        self._selection_handle: Optional[QGraphicsEllipseItem] = None
        self._selection_move_handle: Optional[QGraphicsEllipseItem] = None
        self._selection_scale_handle: Optional[QGraphicsRectItem] = None
        self._rotation_dragging = False
        self._rotation_center: Optional[QPointF] = None
        self._rotation_start_angle: Optional[float] = None
        self._rotation_start_positions: Dict[int, Tuple[float, float]] = {}
        self._scale_dragging = False
        self._scale_anchor: Optional[QPointF] = None
        self._scale_start_handle: Optional[QPointF] = None
        self._scale_start_length = 0.0
        self._scale_start_positions: Dict[int, Tuple[float, float]] = {}
        self._scale_start_text_positions: Dict[TextAnnotationItem, Tuple[QPointF, float]] = {}
        self._scale_has_moved = False
        self._implicit_h_overlays: Dict[int, list[tuple[QGraphicsLineItem, QGraphicsTextItem]]] = {}
        self._group_anchor_overrides: Dict[int, str] = {}

        self._setup_view()
        self._create_paper()
        self._create_overlays()

        self.scene.selectionChanged.connect(self._sync_selection_from_scene)
        self.undo_stack.indexChanged.connect(self._on_undo_stack_changed)
        self.setMouseTracking(True)
        self.viewport().setMouseTracking(True)

    @property
    def graph(self) -> MolGraph:
        """Alias for model for compatibility."""
        return self.model

    def _setup_view(self) -> None:
        self.setBackgroundBrush(QBrush(QColor("#E0E0E0")))
        self.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.setRenderHint(QPainter.RenderHint.TextAntialiasing)

        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)

        self._update_scene_rect()

        self.setDragMode(QGraphicsView.DragMode.NoDrag)

        self._zoom_factor = 1.0
        self._min_zoom = 0.25
        self._max_zoom = 4.0

    def set_show_rulers(self, enabled: bool) -> None:
        self.show_rulers = enabled
        self.viewport().update()

    def set_show_grid(self, enabled: bool) -> None:
        self.show_grid = enabled
        if self._grid_minor_item is not None:
            self._grid_minor_item.setVisible(enabled)
        if self._grid_major_item is not None:
            self._grid_major_item.setVisible(enabled)
        self.viewport().update()

    def paintEvent(self, event) -> None:
        super().paintEvent(event)
        if not self.show_rulers:
            return
        painter = QPainter(self.viewport())
        painter.setRenderHint(QPainter.RenderHint.TextAntialiasing)
        if self.show_rulers:
            self._paint_rulers(painter)
        painter.end()


    def _paint_rulers(self, painter: QPainter) -> None:
        rect = self.viewport().rect()
        thickness = RULER_THICKNESS_PX
        bg = QColor("#F2F2F2")
        border = QColor("#B0B0B0")
        text_color = QColor("#4A4A4A")
        painter.fillRect(0, 0, rect.width(), thickness, bg)
        painter.fillRect(0, 0, thickness, rect.height(), bg)
        painter.setPen(QPen(border, 1))
        painter.drawLine(thickness, thickness, rect.width(), thickness)
        painter.drawLine(thickness, thickness, thickness, rect.height())

        font = painter.font()
        font.setPointSize(8)
        painter.setFont(font)
        painter.setPen(QPen(text_color, 1))

        top_left_scene = self.mapToScene(0, 0)
        bottom_right_scene = self.mapToScene(rect.width(), rect.height())
        x_start = math.floor(top_left_scene.x() / RULER_MAJOR_STEP_PX) * RULER_MAJOR_STEP_PX
        x_end = math.ceil(bottom_right_scene.x() / RULER_MAJOR_STEP_PX) * RULER_MAJOR_STEP_PX
        y_start = math.floor(top_left_scene.y() / RULER_MAJOR_STEP_PX) * RULER_MAJOR_STEP_PX
        y_end = math.ceil(bottom_right_scene.y() / RULER_MAJOR_STEP_PX) * RULER_MAJOR_STEP_PX

        minor_step_view = abs(
            self.mapFromScene(QPointF(RULER_MINOR_STEP_PX, 0)).x()
            - self.mapFromScene(QPointF(0, 0)).x()
        )
        draw_minor = minor_step_view >= 4

        x = x_start
        while x <= x_end:
            x_view = self.mapFromScene(QPointF(x, 0)).x()
            if x_view >= thickness:
                painter.drawLine(
                    x_view, thickness - 1, x_view, thickness - 9
                )
                painter.drawText(x_view + 2, thickness - 11, str(int(x)))
                if draw_minor:
                    minor = x + RULER_MINOR_STEP_PX
                    while minor < x + RULER_MAJOR_STEP_PX and minor <= x_end:
                        minor_view = self.mapFromScene(QPointF(minor, 0)).x()
                        if minor_view >= thickness:
                            painter.drawLine(
                                minor_view, thickness - 1, minor_view, thickness - 5
                            )
                        minor += RULER_MINOR_STEP_PX
            x += RULER_MAJOR_STEP_PX

        y = y_start
        while y <= y_end:
            y_view = self.mapFromScene(QPointF(0, y)).y()
            if y_view >= thickness:
                painter.drawLine(
                    thickness - 1, y_view, thickness - 9, y_view
                )
                painter.drawText(2, y_view - 2, str(int(y)))
                if draw_minor:
                    minor = y + RULER_MINOR_STEP_PX
                    while minor < y + RULER_MAJOR_STEP_PX and minor <= y_end:
                        minor_view = self.mapFromScene(QPointF(0, minor)).y()
                        if minor_view >= thickness:
                            painter.drawLine(
                                thickness - 1, minor_view, thickness - 5, minor_view
                            )
                        minor += RULER_MINOR_STEP_PX
            y += RULER_MAJOR_STEP_PX

    def leaveEvent(self, event) -> None:
        super().leaveEvent(event)

    def _create_paper(self) -> None:
        self.paper = QGraphicsRectItem(0, 0, self.paper_width, self.paper_height)
        self.paper.setBrush(QBrush(Qt.GlobalColor.white))
        self.paper.setPen(QPen(QColor("#CCCCCC"), 1))
        self.paper.setZValue(-10)

        shadow = QGraphicsDropShadowEffect()
        shadow.setBlurRadius(20)
        shadow.setColor(QColor(0, 0, 0, 80))
        shadow.setOffset(5, 5)
        self.paper.setGraphicsEffect(shadow)

        self.scene.addItem(self.paper)
        self._create_grid()
        self.center_on_paper()

    def _create_grid(self) -> None:
        if self._grid_minor_item is not None:
            self.scene.removeItem(self._grid_minor_item)
            self._grid_minor_item = None
        if self._grid_major_item is not None:
            self.scene.removeItem(self._grid_major_item)
            self._grid_major_item = None
        minor_path = QPainterPath()
        for x in range(0, self.paper_width + 1, GRID_MINOR_STEP_PX):
            minor_path.moveTo(x, 0)
            minor_path.lineTo(x, self.paper_height)
        for y in range(0, self.paper_height + 1, GRID_MINOR_STEP_PX):
            minor_path.moveTo(0, y)
            minor_path.lineTo(self.paper_width, y)

        major_path = QPainterPath()
        for x in range(0, self.paper_width + 1, GRID_MAJOR_STEP_PX):
            major_path.moveTo(x, 0)
            major_path.lineTo(x, self.paper_height)
        for y in range(0, self.paper_height + 1, GRID_MAJOR_STEP_PX):
            major_path.moveTo(0, y)
            major_path.lineTo(self.paper_width, y)

        self._grid_minor_item = QGraphicsPathItem(minor_path)
        self._grid_minor_item.setPen(QPen(QColor("#D8D8D8"), 0))
        self._grid_minor_item.setZValue(-9)
        self._grid_minor_item.setVisible(self.show_grid)
        self.scene.addItem(self._grid_minor_item)

        self._grid_major_item = QGraphicsPathItem(major_path)
        self._grid_major_item.setPen(QPen(QColor("#C4C4C4"), 0))
        self._grid_major_item.setZValue(-9)
        self._grid_major_item.setVisible(self.show_grid)
        self.scene.addItem(self._grid_major_item)

    def _create_overlays(self) -> None:
        self._hover_atom_indicator = HoverAtomIndicatorItem()
        self._hover_bond_indicator = HoverBondIndicatorItem()
        self._optimize_zone = OptimizeZoneItem()
        self._preview_bond_item = PreviewBondItem()
        self._preview_ring_item = PreviewRingItem()
        self._preview_chain_item = PreviewChainItem()
        self._preview_chain_label = PreviewChainLabelItem()
        self._preview_arrow_item = PreviewArrowItem()

        self.scene.addItem(self._hover_atom_indicator)
        self.scene.addItem(self._hover_bond_indicator)
        self.scene.addItem(self._optimize_zone)
        self.scene.addItem(self._preview_bond_item)
        self.scene.addItem(self._preview_ring_item)
        self.scene.addItem(self._preview_chain_item)
        self.scene.addItem(self._preview_chain_label)
        self.scene.addItem(self._preview_arrow_item)
        self._overlays_ready = True

    def set_current_tool(self, tool_id: str) -> None:
        if tool_id.startswith("atom_"):
            self.state.default_element = tool_id.split("_", 1)[1]
            tool_id = "tool_atom"
        elif tool_id.startswith("bond_"):
            tool_id = "tool_bond"
        elif tool_id.startswith("tool_brackets_"):
            bracket_key = tool_id.split("tool_brackets_", 1)[1]
            mapping = {
                "round": "()",
                "square": "[]",
                "curly": "{}",
            }
            self.state.active_bracket_type = mapping.get(bracket_key, "[]")
            tool_id = "tool_brackets"

        self.current_tool = tool_id
        self.state.active_tool = tool_id
        self._clear_bond_anchor()
        self._cancel_drag()
        self._arrow_start_pos = None
        self._clear_bracket_preview()

        if tool_id in {"tool_select", "tool_select_lasso"}:
            self.setCursor(Qt.CursorShape.ArrowCursor)
            self.setDragMode(QGraphicsView.DragMode.NoDrag)
        else:
            self.setCursor(Qt.CursorShape.CrossCursor)
            self.setDragMode(QGraphicsView.DragMode.NoDrag)

    def set_active_bond(self, bond_spec: dict) -> None:
        self.state.active_bond_order = bond_spec.get("order", 1)
        self.state.active_bond_style = bond_spec.get("style", BondStyle.PLAIN)
        self.state.active_bond_stereo = bond_spec.get("stereo", BondStereo.NONE)
        self.state.active_bond_mode = bond_spec.get("mode", "increment")
        self.state.active_bond_aromatic = bond_spec.get("aromatic", False)
        if self.state.active_bond_aromatic:
            self.state.active_bond_order = 1
            self.state.active_bond_style = BondStyle.PLAIN
            self.state.active_bond_stereo = BondStereo.NONE
        self.set_current_tool("tool_bond")

    def set_active_ring(self, ring_spec: dict) -> None:
        self.state.active_ring_size = ring_spec.get("size", 6)
        self.state.active_ring_aromatic = ring_spec.get("aromatic", False)
        self.set_current_tool("tool_ring")

    def set_active_element(self, element: str) -> None:
        self.state.default_element = element
        self.set_current_tool("tool_atom")

    def center_on_paper(self) -> None:
        self.centerOn(self.paper_width / 2, self.paper_height / 2)

    def set_paper_size(self, width: int, height: int) -> None:
        width = max(200, int(width))
        height = max(200, int(height))
        if width == self.paper_width and height == self.paper_height:
            return
        self.paper_width = width
        self.paper_height = height
        if getattr(self, "paper", None) is not None:
            self.scene.removeItem(self.paper)
            self.paper = None
        self._update_scene_rect()
        self._create_paper()
        self.viewport().update()

    def mousePressEvent(self, event) -> None:
        scene_pos = self.mapToScene(event.pos())
        self._last_scene_pos = scene_pos

        if not self._is_on_paper(scene_pos.x(), scene_pos.y()) and self.current_tool not in {"tool_select", "tool_select_lasso"}:
            super().mousePressEvent(event)
            return

        if event.button() == Qt.MouseButton.RightButton:
            clicked_item = self._get_item_at(scene_pos)
            if isinstance(
                clicked_item,
                (AtomItem, BondItem, ArrowItem, BracketItem, TextAnnotationItem),
            ):
                if not (event.modifiers() & Qt.KeyboardModifier.ShiftModifier):
                    self.scene.clearSelection()
                clicked_item.setSelected(True)
                self._sync_selection_from_scene()
            self._show_context_menu(
                scene_pos,
                event.globalPosition().toPoint(),
                clicked_item=clicked_item,
            )
            return

        if event.button() != Qt.MouseButton.LeftButton:
            super().mousePressEvent(event)
            return

        if self.current_tool in {"tool_select", "tool_select_lasso"}:
            clicked_item = self._get_item_at(scene_pos)
            if (
                (event.modifiers() & Qt.KeyboardModifier.AltModifier)
                and isinstance(clicked_item, AtomItem)
                and self._cycle_anchor_override(clicked_item.atom_id)
            ):
                return
            if event.button() == Qt.MouseButton.LeftButton and self._hit_selection_scale_handle(scene_pos):
                self._begin_scale_drag(scene_pos)
                return
            if event.button() == Qt.MouseButton.LeftButton and self._hit_selection_handle(scene_pos):
                self._begin_rotation_drag(scene_pos)
                return
            if event.button() == Qt.MouseButton.LeftButton and self._hit_selection_move_handle(scene_pos):
                if self._selection_move_handle is not None:
                    self._selection_move_handle.setCursor(Qt.CursorShape.ClosedHandCursor)
                self._begin_drag(scene_pos)
                return
            if clicked_item is None:
                additive = bool(event.modifiers() & Qt.KeyboardModifier.ShiftModifier)
                free_select = self.current_tool == "tool_select_lasso" or bool(
                    event.modifiers() & Qt.KeyboardModifier.ControlModifier
                )
                self._begin_selection_drag(scene_pos, free_select, additive)
                return
            if event.modifiers() & Qt.KeyboardModifier.ShiftModifier:
                clicked_item.setSelected(not clicked_item.isSelected())
            else:
                if not clicked_item.isSelected():
                    self.scene.clearSelection()
                    clicked_item.setSelected(True)
            self._sync_selection_from_scene()
            if isinstance(clicked_item, AtomItem) and clicked_item.isSelected():
                self._begin_drag(scene_pos)
            return

        clicked_atom_id, clicked_bond_id = self._pick_hover_target(scene_pos)

        arrow_tools = {
            "tool_arrow_forward": "forward",
            "tool_arrow_forward_open": "forward_open",
            "tool_arrow_forward_dashed": "forward_dashed",
            "tool_arrow_retro": "retro",
            "tool_arrow_retro_open": "retro_open",
            "tool_arrow_retro_dashed": "retro_dashed",
            "tool_arrow_both": "both",
            "tool_arrow_both_open": "both_open",
            "tool_arrow_both_dashed": "both_dashed",
            "tool_arrow_equilibrium": "equilibrium",
            "tool_arrow_equilibrium_dashed": "equilibrium_dashed",
            "tool_arrow_retrosynthetic": "retrosynthetic",
            "tool_arrow_curved": "curved",
            "tool_arrow_curved_fishhook": "curved_fishhook",
        }
        if self.current_tool in arrow_tools:
            if self._arrow_start_pos is None:
                self._arrow_start_pos = scene_pos
                return
            if (scene_pos - self._arrow_start_pos).manhattanLength() < 2.0:
                self._arrow_start_pos = None
                self._preview_arrow_item.hide_preview()
                return
            arrow_kind = arrow_tools[self.current_tool]
            cmd = AddArrowCommand(self, self._arrow_start_pos, scene_pos, arrow_kind)
            self.undo_stack.push(cmd)
            self._arrow_start_pos = None
            self._preview_arrow_item.hide_preview()
            return

        if self.current_tool == "tool_brackets":
            self._bracket_drag_start = scene_pos
            preview = self._ensure_bracket_preview()
            preview.setRect(QRectF(scene_pos, scene_pos))
            return

        if self.current_tool in {"tool_charge", "tool_charge_plus", "tool_charge_minus"}:
            if clicked_atom_id is None:
                item = self._get_item_at(scene_pos)
                if isinstance(item, AtomItem):
                    clicked_atom_id = item.atom_id
            if clicked_atom_id is None:
                return
            atom = self.model.get_atom(clicked_atom_id)
            if self.current_tool == "tool_charge_plus":
                new_charge = 1
            elif self.current_tool == "tool_charge_minus":
                new_charge = -1
            else:
                if atom.charge == 0:
                    new_charge = 1
                elif atom.charge > 0:
                    new_charge = -1
                else:
                    new_charge = 0
            cmd = ChangeChargeCommand(self.model, self, clicked_atom_id, new_charge)
            self.undo_stack.push(cmd)
            return

        if self.current_tool == "tool_atom":
            element = self.state.default_element
            if clicked_atom_id is None:
                cmd = AddAtomCommand(self.model, self, element, scene_pos.x(), scene_pos.y())
                self.undo_stack.push(cmd)
            else:
                cmd = ChangeAtomCommand(self.model, self, clicked_atom_id, element)
                self.undo_stack.push(cmd)
            return

        if self.current_tool == "tool_bond":
            if clicked_bond_id is not None and clicked_atom_id is None:
                if self.state.active_bond_mode == "increment" and not self.state.active_bond_aromatic:
                    self._cycle_bond_order(clicked_bond_id)
                else:
                    self._apply_bond_style(clicked_bond_id)
                return

            if clicked_atom_id is not None:
                self._begin_place_bond(clicked_atom_id, scene_pos)
                return
            self._begin_place_bond(None, scene_pos)
            return

        if self.current_tool == "tool_ring":
            if clicked_bond_id is not None and clicked_atom_id is None:
                self._begin_place_ring("bond", clicked_bond_id, scene_pos)
                return
            if clicked_atom_id is not None:
                self._begin_place_ring("atom", clicked_atom_id, scene_pos)
                return
            self._begin_place_ring("free", None, scene_pos)
            return

        if self.current_tool == "tool_chain":
            if clicked_atom_id is not None:
                self._begin_place_chain(clicked_atom_id, scene_pos)
                return
            if clicked_bond_id is None:
                self._begin_place_chain(None, scene_pos)
                return

        if self.current_tool == "tool_erase":
            if clicked_atom_id is not None:
                self._delete_selection({clicked_atom_id}, set())
            elif clicked_bond_id is not None:
                self._delete_selection(set(), {clicked_bond_id})
            return

        elif self.current_tool == "tool_text":
            # Create text annotation at position
            item = TextAnnotationItem(
                "",
                scene_pos.x(),
                scene_pos.y()
            )
            self._apply_text_settings(item)
            # Enter edit mode immediately
            item.setTextInteractionFlags(Qt.TextInteractionFlag.TextEditorInteraction)
            self.scene.addItem(item)
            item.setFocus()
            # Position cursor at end to ensure it's visible
            cursor = item.textCursor()
            cursor.movePosition(cursor.MoveOperation.End)
            item.setTextCursor(cursor)
            
            self.scene.clearSelection()
            item.setSelected(True)
            return

        super().mousePressEvent(event)

    def mouseMoveEvent(self, event) -> None:
        scene_pos = self.mapToScene(event.pos())
        self._last_scene_pos = scene_pos

        if self._scale_dragging:
            self._update_scale_drag(scene_pos)
            return
        if self._rotation_dragging:
            self._update_rotation_drag(scene_pos)
            return

        if self._select_drag_mode is not None:
            self._update_selection_drag(scene_pos)
            return

        arrow_tools = {
            "tool_arrow_forward": "forward",
            "tool_arrow_forward_open": "forward_open",
            "tool_arrow_forward_dashed": "forward_dashed",
            "tool_arrow_retro": "retro",
            "tool_arrow_retro_open": "retro_open",
            "tool_arrow_retro_dashed": "retro_dashed",
            "tool_arrow_both": "both",
            "tool_arrow_both_open": "both_open",
            "tool_arrow_both_dashed": "both_dashed",
            "tool_arrow_equilibrium": "equilibrium",
            "tool_arrow_equilibrium_dashed": "equilibrium_dashed",
            "tool_arrow_retrosynthetic": "retrosynthetic",
            "tool_arrow_curved": "curved",
            "tool_arrow_curved_fishhook": "curved_fishhook",
        }
        if self.current_tool in arrow_tools and self._arrow_start_pos is not None:
            arrow_kind = arrow_tools[self.current_tool]
            self._preview_arrow_item.update_preview(self._arrow_start_pos, scene_pos, arrow_kind)
            return

        if self._drag_mode == "place_bond":
            self._update_bond_drag_state(scene_pos)
            self._update_bond_preview(scene_pos, event.modifiers())
            return
        if self._drag_mode == "place_ring":
            self._update_ring_preview(scene_pos, event.modifiers())
            return
        if self._drag_mode == "place_chain":
            self._update_chain_preview(scene_pos, event.modifiers())
            return

        if self._bracket_drag_start is not None and self.current_tool == "tool_brackets":
            preview = self._ensure_bracket_preview()
            preview.setRect(QRectF(self._bracket_drag_start, scene_pos).normalized())
            return

        if self._dragging_selection and self._drag_start_pos is not None:
            delta = scene_pos - self._drag_start_pos
            if delta.manhattanLength() > 0:
                self._drag_has_moved = True
                
                # Move atoms
                for atom_id, (x, y) in self._drag_start_positions.items():
                    nx = x + delta.x()
                    ny = y + delta.y()
                    self.model.update_atom_position(atom_id, nx, ny)
                    self.update_atom_item(atom_id, nx, ny)
                self.update_bond_items_for_atoms(set(self._drag_start_positions.keys()))
                
                # Move text items
                if hasattr(self, "_drag_start_text_positions"):
                    for item, (pos, rot) in self._drag_start_text_positions.items():
                        new_pos = pos + delta
                        item.setPos(new_pos)
                
                self._update_selection_overlay()
            return

        self._update_hover(scene_pos)

        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event) -> None:
        if self._scale_dragging:
            self._finalize_scale_drag()
            return
        if self._rotation_dragging:
            self._finalize_rotation_drag()
            return
        if self._select_drag_mode is not None:
            self._finalize_selection_drag()
            return
        if self._drag_mode == "place_bond":
            self._finalize_bond(event.modifiers())
            return
        if self._drag_mode == "place_ring":
            self._finalize_ring(event.modifiers())
            return
        if self._drag_mode == "place_chain":
            self._finalize_chain(event.modifiers())
            return

        if self._bracket_drag_start is not None and self.current_tool == "tool_brackets":
            rect = QRectF(self._bracket_drag_start, self._last_scene_pos).normalized()
            self._clear_bracket_preview()
            if rect.width() < 4.0 or rect.height() < 4.0:
                return
            items = [
                item
                for item in self.scene.items(rect)
                if isinstance(item, (AtomItem, BondItem))
            ]
            if items:
                bbox = items[0].sceneBoundingRect()
                for item in items[1:]:
                    bbox = bbox.united(item.sceneBoundingRect())
            else:
                bbox = rect
            cmd = AddBracketCommand(self, bbox, self.state.active_bracket_type)
            self.undo_stack.push(cmd)
            return

        if self._dragging_selection:
            if self._selection_move_handle is not None:
                self._selection_move_handle.setCursor(Qt.CursorShape.OpenHandCursor)
            if self._drag_has_moved:
                # Capture end states
                move_atoms = bool(self._drag_start_positions)
                move_text = bool(getattr(self, "_drag_start_text_positions", {}))
                
                if move_atoms and move_text:
                    self.undo_stack.beginMacro("Move selection")
                
                if move_atoms:
                    after = {
                        atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
                        for atom_id in self._drag_start_positions
                    }
                    cmd = MoveAtomsCommand(
                        self.model,
                        self,
                        self._drag_start_positions,
                        after,
                        skip_first_redo=True,
                    )
                    self.undo_stack.push(cmd)
                    
                if move_text:
                    before = self._drag_start_text_positions
                    after = {
                        item: (item.pos(), item.rotation())
                        for item in before.keys()
                    }
                    cmd = MoveTextItemsCommand(self, before, after)
                    # We already moved them, so don't re-execute redo immediately if command stack auto-redoes?
                    # QUndoStack.push() calls redo().
                    # MoveItemsCommand.redo() sets position.
                    # It's fine to set it again to same value.
                    self.undo_stack.push(cmd)
                    
                if move_atoms and move_text:
                    self.undo_stack.endMacro()

            self._dragging_selection = False
            self._drag_start_pos = None
            self._drag_start_positions = {}
            self._drag_start_text_positions = {}
            self._drag_has_moved = False
            return

        super().mouseReleaseEvent(event)
    
    def _delete_selection(
        self,
        atom_ids: set[int],
        bond_ids: set[int],
        arrow_items: Iterable = (),
        bracket_items: Iterable = (),
        text_items: Iterable = (),
    ) -> None:
        # Prioritize delete logic
        cmd = DeleteSelectionCommand(
            self.model,
            self,
            atom_ids,
            bond_ids,
            arrow_items=arrow_items,
            bracket_items=bracket_items,
            text_items=text_items,
        )
        self.undo_stack.push(cmd)
        self.scene.clearSelection()

    def _get_item_at(self, scene_pos: QPointF):
        for item in self.scene.items(scene_pos):
             if isinstance(item, (AtomItem, BondItem, ArrowItem, BracketItem, TextAnnotationItem)):
                return item
        return None

    def mouseDoubleClickEvent(self, event) -> None:
        scene_pos = self.mapToScene(event.pos())
        item = self._get_item_at(scene_pos)
        if isinstance(item, BondItem):
            bond = self.model.get_bond(item.bond_id)
            current = bond.length_px if bond.length_px is not None else self.state.bond_length
            value, ok = QInputDialog.getDouble(
                self,
                "Longitud de enlace",
                "Longitud (px, 0 para usar global):",
                float(current),
                0.0,
                999.0,
                1,
            )
            if ok:
                new_length = None if value <= 0 else float(value)
                cmd = ChangeBondLengthCommand(self.model, self, bond.id, new_length)
                self.undo_stack.push(cmd)
            return
        if isinstance(item, AtomItem):
            atom = self.model.get_atom(item.atom_id)
            value, anchor = self._prompt_atom_label(
                atom.element, atom_id=atom.id, initial=atom.element
            )
            if value is None:
                return
            cmd = ChangeAtomCommand(
                self.model,
                self,
                atom.id,
                value,
                anchor_override=anchor,
            )
            self.undo_stack.push(cmd)
            return
        super().mouseDoubleClickEvent(event)

    def wheelEvent(self, event: QWheelEvent) -> None:
        if event.angleDelta().y() > 0:
            self.zoom_in()
        else:
            self.zoom_out()

    def keyPressEvent(self, event) -> None:
        if event.key() in (Qt.Key.Key_Delete, Qt.Key.Key_Backspace):
            has_selected_arrows = any(
                isinstance(item, ArrowItem) for item in self.scene.selectedItems()
            )
            has_selected_brackets = any(
                isinstance(item, BracketItem) for item in self.scene.selectedItems()
            )
            if (
                self.state.selected_atoms
                or self.state.selected_bonds
                or has_selected_arrows
                or has_selected_brackets
            ):
                self.delete_selection()
                return
            if self._delete_hovered():
                return
        if self._handle_atom_text_entry(event):
            return
        if self._handle_select_all(event):
            return
        if self._handle_nudge(event):
            return
        if self._handle_hotkeys(event):
            return
        super().keyPressEvent(event)

    def _handle_atom_text_entry(self, event) -> bool:
        if event.modifiers() & (
            Qt.KeyboardModifier.ControlModifier
            | Qt.KeyboardModifier.AltModifier
            | Qt.KeyboardModifier.MetaModifier
        ):
            return False
        if not self.state.selected_atoms or self.state.selected_bonds:
            return False
        if len(self.state.selected_atoms) != 1:
            return False
        typed = event.text()
        if not typed or not typed.isalpha():
            return False

        atom_id = next(iter(self.state.selected_atoms))
        atom = self.model.get_atom(atom_id)
        seed = self._normalize_atom_label(typed) or typed
        value, anchor = self._prompt_atom_label(
            atom.element, atom_id=atom_id, initial=seed
        )
        if value is None:
            return True
        cmd = ChangeAtomCommand(
            self.model,
            self,
            atom_id,
            value,
            anchor_override=anchor,
        )
        self.undo_stack.push(cmd)
        return True

    def _handle_nudge(self, event) -> bool:
        key = event.key()
        if key not in (
            Qt.Key.Key_Left,
            Qt.Key.Key_Right,
            Qt.Key.Key_Up,
            Qt.Key.Key_Down,
        ):
            return False

        focus_item = self.scene.focusItem()
        if isinstance(focus_item, TextAnnotationItem) and (
            focus_item.textInteractionFlags()
            & Qt.TextInteractionFlag.TextEditorInteraction
        ):
            return False

        selected_atom_ids = self._selected_atom_ids_for_transform()
        selected_text_items = self._selected_text_items()
        if not selected_atom_ids and not selected_text_items:
            return False

        step = 1.0
        if event.modifiers() & Qt.KeyboardModifier.ShiftModifier:
            step = 10.0

        dx = 0.0
        dy = 0.0
        if key == Qt.Key.Key_Left:
            dx = -step
        elif key == Qt.Key.Key_Right:
            dx = step
        elif key == Qt.Key.Key_Up:
            dy = -step
        elif key == Qt.Key.Key_Down:
            dy = step

        if selected_atom_ids:
            before = {
                atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
                for atom_id in selected_atom_ids
                if atom_id in self.model.atoms
            }
            after = {
                atom_id: (x + dx, y + dy)
                for atom_id, (x, y) in before.items()
            }
            self.undo_stack.push(MoveAtomsCommand(self.model, self, before, after))

        if selected_text_items:
            before = {item: (item.pos(), item.rotation()) for item in selected_text_items}
            after = {
                item: (QPointF(pos.x() + dx, pos.y() + dy), rot)
                for item, (pos, rot) in before.items()
            }
            self.undo_stack.push(MoveTextItemsCommand(self, before, after))

        self._update_selection_overlay()
        return True

    def _select_all_items(self) -> None:
        self.scene.clearSelection()
        for item in self.scene.items():
            if isinstance(
                item,
                (AtomItem, BondItem, ArrowItem, BracketItem, TextAnnotationItem),
            ) and item.isVisible():
                item.setSelected(True)
        self._sync_selection_from_scene()

    def _handle_select_all(self, event) -> bool:
        if not (event.modifiers() & (Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.MetaModifier)):
            return False
        if event.key() not in (Qt.Key.Key_A, Qt.Key.Key_E):
            return False
        focus_item = self.scene.focusItem()
        if isinstance(focus_item, TextAnnotationItem) and (
            focus_item.textInteractionFlags()
            & Qt.TextInteractionFlag.TextEditorInteraction
        ):
            return False
        self._select_all_items()
        return True

    @staticmethod
    def _normalize_atom_label(text: str) -> str | None:
        cleaned = text.strip()
        if not cleaned:
            return None
        allowed = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+-()_^")
        if any(ch not in allowed for ch in cleaned):
            return None
        alias = FUNCTIONAL_GROUP_ALIASES.get(cleaned.lower())
        if alias:
            return alias
        if "_" in cleaned or "^" in cleaned:
            return cleaned
        if cleaned.isalpha() and len(cleaned) <= 2:
            if len(cleaned) == 1:
                normalized = cleaned.upper()
            else:
                normalized = cleaned[0].upper() + cleaned[1:].lower()
            if normalized in ELEMENT_SYMBOLS:
                return normalized
        return cleaned

    def _atom_label_items(self) -> list[str]:
        elements = ["C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "H"]
        return elements + FUNCTIONAL_GROUP_LABELS

    def _label_anchor_candidates(self, label: str) -> list[str]:
        cleaned = label.strip()
        if not cleaned:
            return []
        tokens = self._group_label_tokens(cleaned)
        if tokens:
            seen: set[str] = set()
            ordered: list[str] = []
            for symbol, _token in tokens:
                if symbol in seen:
                    continue
                seen.add(symbol)
                ordered.append(symbol)
            return ordered
        ordered: list[str] = []
        seen: set[str] = set()
        i = 0
        while i < len(cleaned):
            ch = cleaned[i]
            if not ch.isalpha():
                i += 1
                continue
            if ch.isupper():
                symbol = None
                if i + 1 < len(cleaned) and cleaned[i + 1].islower():
                    candidate = cleaned[i : i + 2]
                    if candidate in ELEMENT_SYMBOLS:
                        symbol = candidate
                        i += 2
                if symbol is None and ch in ELEMENT_SYMBOLS:
                    symbol = ch
                    i += 1
                if symbol is None:
                    i += 1
                    continue
                if symbol not in seen:
                    seen.add(symbol)
                    ordered.append(symbol)
                continue
            i += 1
        return ordered

    def _prompt_anchor_for_atom(self, atom_id: int) -> None:
        atom = self.model.get_atom(atom_id)
        if atom is None:
            return
        label = atom.element
        candidates = self._label_anchor_candidates(label)
        if not candidates:
            return
        items = ["Auto"] + candidates
        current = self._group_anchor_overrides.get(atom_id) or "Auto"
        index = items.index(current) if current in items else 0
        value, ok = QInputDialog.getItem(
            self,
            "Átomo de unión",
            "Átomo de unión:",
            items,
            index,
            True,
        )
        if not ok:
            return
        value = value.strip()
        anchor = None
        if value and value != "Auto":
            anchor = self._normalize_atom_label(value)
        self.set_anchor_override(atom_id, anchor)
        self._refresh_atom_label(atom_id)

    def _cycle_anchor_override(self, atom_id: int) -> bool:
        atom = self.model.get_atom(atom_id)
        if atom is None:
            return False
        if atom.element in ELEMENT_SYMBOLS:
            return False
        candidates = self._label_anchor_candidates(atom.element)
        if not candidates:
            return False
        current = self._group_anchor_overrides.get(atom_id)
        if current in candidates:
            idx = (candidates.index(current) + 1) % len(candidates)
        else:
            idx = 0
        self.set_anchor_override(atom_id, candidates[idx])
        self._refresh_atom_label(atom_id)
        return True

    def _prompt_atom_label(
        self, current: str, atom_id: Optional[int] = None, initial: str | None = None
    ) -> tuple[Optional[str], Optional[str]]:
        items = self._atom_label_items()
        seed = (initial or "").strip() or current
        if seed:
            if seed not in items:
                items = [seed] + items
            current_index = items.index(seed)
        else:
            current_index = 0
        anchor = self._group_anchor_overrides.get(atom_id) if atom_id is not None else None
        dialog = AtomLabelDialog(
            seed,
            anchor,
            items,
            ELEMENT_SYMBOLS,
            self,
        )
        if not dialog.exec():
            return None, None
        value, anchor_value = dialog.value()
        normalized = self._normalize_atom_label(value)
        if not normalized:
            return None, None
        anchor_normalized = None
        if anchor_value:
            anchor_normalized = self._normalize_atom_label(anchor_value)
        # If the user chose an anchor for a group label, reorder the label now
        # so the visible text reflects the chosen attachment atom.
        if (
            anchor_normalized
            and atom_id is not None
            and normalized not in ELEMENT_SYMBOLS
        ):
            normalized, _ = self._reflow_group_label(
                normalized, atom_id, anchor_normalized
            )
        return normalized, anchor_normalized

    def get_persistence_data(self) -> dict:
        """Collects canvas settings and non-structural items for serialization."""
        data = {
            "settings": {
                "paper_width": self.paper_width,
                "paper_height": self.paper_height,
                "show_grid": self.show_grid,
                "show_rulers": self.show_rulers,
                "bond_length": self.state.bond_length,
                "use_aromatic_circles": self.state.use_aromatic_circles,
                "show_carbons": self.state.show_implicit_carbons,
                "show_hydrogens": self.state.show_implicit_hydrogens,
                "use_element_colors": self.state.use_element_colors,
                "font_family": self.state.label_font_family,
                "font_size": self.state.label_font_size,
                "font_bold": self.state.label_font_bold,
                "font_italic": self.state.label_font_italic,
                "font_underline": self.state.label_font_underline,
            },
            "annotations": {
                "arrows": [],
                "brackets": [],
                "text_items": []
            },
            "label_anchors": {
                str(atom_id): anchor for atom_id, anchor in self._group_anchor_overrides.items()
            },
        }

        for item in self.scene.items():
            if isinstance(item, ArrowItem):
                if isinstance(item, PreviewArrowItem):
                    continue
                data["annotations"]["arrows"].append({
                    "kind": item.kind(),
                    "start": {"x": item.start_point().x(), "y": item.start_point().y()},
                    "end": {"x": item.end_point().x(), "y": item.end_point().y()}
                })
            elif isinstance(item, BracketItem):
                data["annotations"]["brackets"].append({
                    "kind": item._kind,
                    "rect": {
                        "x": item._base_rect.x(),
                        "y": item._base_rect.y(),
                        "w": item._base_rect.width(),
                        "h": item._base_rect.height()
                    },
                    "padding": item._padding
                })
            elif isinstance(item, TextAnnotationItem):
                data["annotations"]["text_items"].append({
                    "text": item.toPlainText(),
                    "html": item.toHtml(), # Store HTML for rich text (formatting, colors)
                    "x": item.pos().x(),
                    "y": item.pos().y(),
                    "rotation": item.rotation(),
                    "font": item.font().toString(),
                    "color": item.defaultTextColor().name()
                })
        return data

    def load_persistence_data(self, data: dict) -> None:
        """Restores canvas settings and non-structural items from a dictionary."""
        settings = data.get("settings", {})
        self.paper_width = settings.get("paper_width", DEFAULT_PAPER_WIDTH)
        self.paper_height = settings.get("paper_height", DEFAULT_PAPER_HEIGHT)
        self.set_show_grid(settings.get("show_grid", False))
        self.set_show_rulers(settings.get("show_rulers", False))
        self.state.bond_length = settings.get("bond_length", DEFAULT_BOND_LENGTH)
        self.state.use_aromatic_circles = settings.get("use_aromatic_circles", False)
        self.state.show_implicit_carbons = settings.get("show_carbons", False)
        self.state.show_implicit_hydrogens = settings.get("show_hydrogens", False)
        self.state.use_element_colors = settings.get("use_element_colors", True)
        self.state.label_font_family = settings.get("font_family", "Arial")
        self.state.label_font_size = settings.get("font_size", 11.0)
        self.state.label_font_bold = settings.get("font_bold", False)
        self.state.label_font_italic = settings.get("font_italic", False)
        self.state.label_font_underline = settings.get("font_underline", False)

        self._group_anchor_overrides = {
            int(atom_id): anchor
            for atom_id, anchor in data.get("label_anchors", {}).items()
            if anchor
        }

        # Re-apply paper size visual
        self._create_paper()

        # Restore Annotations
        annotations = data.get("annotations", {})
        
        for arrow_d in annotations.get("arrows", []):
            start = QPointF(arrow_d["start"]["x"], arrow_d["start"]["y"])
            end = QPointF(arrow_d["end"]["x"], arrow_d["end"]["y"])
            arrow = ArrowItem(start, end, kind=arrow_d["kind"])
            self.scene.addItem(arrow)

        for br_d in annotations.get("brackets", []):
            rect_d = br_d["rect"]
            rect = QRectF(rect_d["x"], rect_d["y"], rect_d["w"], rect_d["h"])
            bracket = BracketItem(rect, kind=br_d["kind"], padding=br_d.get("padding", 8.0))
            self.scene.addItem(bracket)

        for txt_d in annotations.get("text_items", []):
            text_item = TextAnnotationItem(txt_d["text"], txt_d["x"], txt_d["y"])
            if "html" in txt_d:
                text_item.setHtml(txt_d["html"])
            if "rotation" in txt_d:
                text_item.setRotation(txt_d["rotation"])
            if "font" in txt_d:
                font = QFont()
                font.fromString(txt_d["font"])
                text_item.setFont(font)
            if "color" in txt_d:
                text_item.setDefaultTextColor(QColor(txt_d["color"]))
            self.scene.addItem(text_item)

        # Full refresh to update atom visibility and circles
        self.refresh_atom_visibility()
        self.refresh_aromatic_circles()

    def zoom_in(self) -> None:
        if self._zoom_factor < self._max_zoom:
            self._zoom_factor *= 1.2
            self.scale(1.2, 1.2)

    def zoom_out(self) -> None:
        if self._zoom_factor > self._min_zoom:
            self._zoom_factor /= 1.2
            self.scale(1 / 1.2, 1 / 1.2)

    def clear_canvas(self) -> None:
        self._overlays_ready = False
        self._cancel_drag()
        self.undo_stack.blockSignals(True)
        
        # Nullify Python references to scene items before scene.clear()
        # so that subsequent recreation logic doesn't try to access deleted C++ objects.
        self.paper = None
        self._grid_minor_item = None
        self._grid_major_item = None
        self._selection_box = None
        self._selection_handle = None
        self._selection_move_handle = None
        self._select_preview_path = None
        self._select_preview_rect = None
        self._bracket_preview = None
        self._hover_atom_indicator = None
        self._hover_bond_indicator = None
        self._optimize_zone = None
        self._preview_bond_item = None
        self._preview_ring_item = None
        self._preview_chain_item = None
        self._preview_chain_label = None
        self._preview_arrow_item = None

        self.scene.clear()
        self.model.clear()
        self.undo_stack.clear()
        self.undo_stack.blockSignals(False)
        self.atom_items.clear()
        self.bond_items.clear()
        self.aromatic_circles.clear()
        self.arrow_items.clear()
        self.bracket_items.clear()
        self._implicit_h_overlays.clear()
        self._group_anchor_overrides.clear()
        self._ring_centers.clear()
        self._next_ring_id = 1
        self.state.selected_atoms.clear()
        self.state.selected_bonds.clear()
        self.bond_anchor_id = None
        self.hovered_atom_id = None
        self.hovered_bond_id = None
        
        self._create_paper()
        self._create_overlays()

    def _rebuild_items_from_model(self) -> None:
        """Create AtomItems and BondItems for everything currently in self.model."""
        self.atom_items.clear()
        self.bond_items.clear()
        
        # 1. Create AtomItems
        for atom in self.model.atoms.values():
            item = AtomItem(
                atom,
                radius=ATOM_HIT_RADIUS,
                show_carbon=self.state.show_implicit_carbons,
                show_hydrogen=self.state.show_implicit_hydrogens,
                label_font=QFont(self.state.label_font_family, int(self.state.label_font_size)),
                style=self.drawing_style,
                use_element_colors=self.state.use_element_colors
            )
            self.scene.addItem(item)
            self.atom_items[atom.id] = item
        
        # 2. Recalculate ring centers for double bond offsets
        self.refresh_ring_centers()
            
        # 3. Create BondItems with ring context
        for bond in self.model.bonds.values():
            a1 = self.model.atoms.get(bond.a1_id)
            a2 = self.model.atoms.get(bond.a2_id)
            if a1 and a2:
                item = BondItem(
                    bond, a1, a2,
                    render_aromatic_as_circle=self.state.use_aromatic_circles,
                    style=self.drawing_style
                )
                if bond.ring_id is not None:
                    item.set_ring_context(self._ring_centers.get(bond.ring_id))
                item.set_offset_sign(self._bond_offset_sign(bond))
                self._configure_bond_rendering(bond, item)
                # Explicitly update positions now that context and sign are set
                item.update_positions(a1, a2)
                self.scene.addItem(item)
                self.bond_items[bond.id] = item
        
        # 4. Global refresh: labels, visibility, shrinks, and implicit H overlays
        self.refresh_atom_visibility()
        self.refresh_aromatic_circles()

    def refresh_ring_centers(self) -> None:
        """Recalculate geometric centers for all rings in the model."""
        self._ring_centers.clear()
        ring_to_atoms = {}
        for bond in self.model.bonds.values():
            if bond.ring_id is not None:
                atoms = ring_to_atoms.setdefault(bond.ring_id, set())
                atoms.add(bond.a1_id)
                atoms.add(bond.a2_id)
        
        for ring_id, atom_ids in ring_to_atoms.items():
            sum_x = 0.0
            sum_y = 0.0
            count = 0
            for aid in atom_ids:
                atom = self.model.get_atom(aid)
                if atom:
                    sum_x += atom.x
                    sum_y += atom.y
                    count += 1
            if count > 0:
                self.register_ring_center(ring_id, (sum_x / count, sum_y / count))

    def add_text_item(self, item: TextAnnotationItem) -> None:
        if item.scene() is not self.scene:
            self.scene.addItem(item)
        item.setSelected(True)

    def remove_text_item(self, item: TextAnnotationItem) -> None:
        if item.scene() is self.scene:
            self.scene.removeItem(item)

    def readd_text_item(self, item: TextAnnotationItem) -> None:
        if item.scene() is not self.scene:
            self.scene.addItem(item)

    def delete_selection(self) -> None:
        selected_arrows = [
            item for item in self.scene.selectedItems() if isinstance(item, ArrowItem)
        ]
        selected_brackets = [
            item for item in self.scene.selectedItems() if isinstance(item, BracketItem)
        ]
        selected_text_items = [
            item for item in self.scene.selectedItems() if isinstance(item, TextAnnotationItem)
        ]
        if (
            not self.state.selected_atoms
            and not self.state.selected_bonds
            and not selected_arrows
            and not selected_brackets
            and not selected_text_items
        ):
            return
        self._delete_selection(
            set(self.state.selected_atoms),
            set(self.state.selected_bonds),
            selected_arrows,
            selected_brackets,
            selected_text_items,
        )

    def _delete_hovered(self) -> bool:
        if self.hovered_atom_id is not None:
            self._delete_selection({self.hovered_atom_id}, set())
            return True
        if self.hovered_bond_id is not None:
            self._delete_selection(set(), {self.hovered_bond_id})
            return True
        return False



    def copy_to_clipboard(self) -> None:
        selected_text_items = self._selected_text_items()
        selected_arrows = any(isinstance(item, ArrowItem) for item in self.scene.selectedItems())
        selected_brackets = any(isinstance(item, BracketItem) for item in self.scene.selectedItems())
        has_structure_selection = bool(self.state.selected_atoms or self.state.selected_bonds)
        only_text_selection = bool(selected_text_items) and not (
            has_structure_selection or selected_arrows or selected_brackets
        )
        if not self.model.atoms and not selected_text_items:
            return
        mime = QMimeData()
        if self.model.atoms and not only_text_selection:
            try:
                molfile = molgraph_to_molfile(self.model)
                smiles = molgraph_to_smiles(self.model)
                mime.setData("chemical/x-mdl-molfile", molfile.encode("utf-8"))
                if smiles:
                    mime.setText(smiles)
            except Exception:
                pass
        if only_text_selection:
            data = {
                "text_items": [
                    {
                        "text": item.toPlainText(),
                        "html": item.toHtml(),
                        "x": item.pos().x(),
                        "y": item.pos().y(),
                        "rotation": item.rotation(),
                        "font": item.font().toString(),
                        "color": item.defaultTextColor().name(),
                    }
                    for item in selected_text_items
                ]
            }
            mime.setData(
                "application/x-chemuson-text-items",
                json.dumps(data).encode("utf-8"),
            )

        has_selection = bool(self.state.selected_atoms or self.state.selected_bonds)
        has_selection = has_selection or bool(self.scene.selectedItems())
        image = self._render_scene_image(
            scale=CLIPBOARD_RENDER_SCALE,
            selected_only=has_selection,
        )
        if image is not None:
            buffer = QBuffer()
            buffer.open(QBuffer.OpenModeFlag.WriteOnly)
            image.save(buffer, "PNG")
            mime.setData("image/png", buffer.data())
            mime.setImageData(image)

        QApplication.clipboard().setMimeData(mime)

    def paste_from_clipboard(self) -> None:
        clipboard = QApplication.clipboard()
        mime = clipboard.mimeData()
        if mime is None:
            return

        if mime.hasFormat("application/x-chemuson-text-items"):
            try:
                payload = json.loads(bytes(mime.data("application/x-chemuson-text-items")).decode("utf-8"))
                self._paste_text_items(payload)
                return
            except Exception:
                pass

        if mime.hasFormat("chemical/x-mdl-molfile"):
            molfile = bytes(mime.data("chemical/x-mdl-molfile")).decode("utf-8", errors="ignore")
            try:
                graph = molfile_to_molgraph(molfile)
                self._insert_molgraph(graph)
                return
            except Exception:
                pass

        if mime.hasText():
            smiles = mime.text().strip()
            if smiles:
                try:
                    graph = smiles_to_molgraph(smiles)
                    self._insert_molgraph(graph)
                    return
                except Exception:
                    pass

        if mime.hasFormat("image/png") or mime.hasImage() or mime.hasFormat("image/svg+xml"):
            self._insert_image_from_clipboard(mime)

    def _paste_text_items(self, payload: dict) -> None:
        items = payload.get("text_items", [])
        if not items:
            return
        created: list[TextAnnotationItem] = []
        for txt_d in items:
            text_item = TextAnnotationItem(txt_d.get("text", ""), 0.0, 0.0)
            html = txt_d.get("html")
            if html:
                text_item.setHtml(html)
            if "rotation" in txt_d:
                text_item.setRotation(float(txt_d["rotation"]))
            if "font" in txt_d:
                font = QFont()
                font.fromString(txt_d["font"])
                text_item.setFont(font)
            if "color" in txt_d:
                text_item.setDefaultTextColor(QColor(txt_d["color"]))
            text_item.setPos(float(txt_d.get("x", 0.0)), float(txt_d.get("y", 0.0)))
            self.scene.addItem(text_item)
            created.append(text_item)

        bbox: Optional[QRectF] = None
        for item in created:
            rect = item.sceneBoundingRect()
            bbox = rect if bbox is None else bbox.united(rect)
        if bbox is None:
            return
        target = self._last_scene_pos
        if target is None:
            target = self.mapToScene(self.viewport().rect().center())
        dx = target.x() - bbox.left()
        dy = target.y() - bbox.top()
        for item in created:
            item.setPos(item.pos() + QPointF(dx, dy))

        self.scene.clearSelection()
        for item in created:
            item.setSelected(True)
        self._sync_selection_from_scene()

    def _render_scene_bounds(self, selected_only: bool = False) -> Optional[QRectF]:
        rect: Optional[QRectF] = None

        def extend(candidate: QRectF) -> None:
            nonlocal rect
            if not candidate.isValid() or candidate.isNull():
                return
            if rect is None:
                rect = candidate
            else:
                rect = rect.united(candidate)

        def extend_atom_bounds(atom_id: int) -> None:
            item = self.atom_items.get(atom_id)
            if item is None or item.scene() is not self.scene:
                return
            if item.pen().style() != Qt.PenStyle.NoPen or item.brush().style() != Qt.BrushStyle.NoBrush:
                extend(item.sceneBoundingRect())
            if item.label.isVisible():
                extend(item.label.sceneBoundingRect())
            if item.charge_label.isVisible():
                extend(item.charge_label.sceneBoundingRect())
            overlays = self._implicit_h_overlays.get(atom_id)
            if overlays:
                for line_item, text_item in overlays:
                    if line_item.scene() is self.scene and line_item.isVisible():
                        extend(line_item.sceneBoundingRect())
                    if text_item.scene() is self.scene and text_item.isVisible():
                        extend(text_item.sceneBoundingRect())

        if selected_only:
            for atom_id in self._selected_atom_ids_for_transform():
                extend_atom_bounds(atom_id)
            for bond_id in self.state.selected_bonds:
                item = self.bond_items.get(bond_id)
                if item is None or item.scene() is not self.scene:
                    continue
                extend(item.sceneBoundingRect())
            for item in self.scene.selectedItems():
                if isinstance(item, (ArrowItem, BracketItem, TextAnnotationItem)):
                    extend(item.sceneBoundingRect())
        else:
            for atom_id in self.atom_items.keys():
                extend_atom_bounds(atom_id)
            for item in self.bond_items.values():
                if item.scene() is not self.scene:
                    continue
                if item.isVisible():
                    extend(item.sceneBoundingRect())
            for item in self.aromatic_circles:
                if item.scene() is not self.scene:
                    continue
                if item.isVisible():
                    extend(item.sceneBoundingRect())
            for item in self.arrow_items:
                if item.scene() is not self.scene:
                    continue
                if item.isVisible():
                    extend(item.sceneBoundingRect())
            for item in self.bracket_items:
                if item.scene() is not self.scene:
                    continue
                if item.isVisible():
                    extend(item.sceneBoundingRect())
            for item in self.scene.items():
                if isinstance(item, TextAnnotationItem) and item.isVisible():
                    extend(item.sceneBoundingRect())

        if rect is None:
            return None
        pad = max(1.0, self.drawing_style.stroke_px)
        return rect.adjusted(-pad, -pad, pad, pad)

    def _hidden_render_items(self) -> list:
        items = [getattr(self, "paper", None)]
        overlay_attrs = [
            "_hover_atom_indicator",
            "_hover_bond_indicator",
            "_optimize_zone",
            "_preview_bond_item",
            "_preview_ring_item",
            "_preview_chain_item",
            "_preview_chain_label",
            "_preview_arrow_item",
            "_grid_minor_item",
            "_grid_major_item",
        ]
        for attr in overlay_attrs:
            items.append(getattr(self, attr, None))
        items.append(getattr(self, "_bracket_preview", None))
        items.append(getattr(self, "_select_preview_path", None))
        items.append(getattr(self, "_select_preview_rect", None))
        return [
            item
            for item in items
            if item is not None and item.scene() is self.scene
        ]

    def _with_hidden_render_items(self, render_fn):
        hidden = []
        for item in self._hidden_render_items():
            if item.isVisible():
                hidden.append(item)
                item.setVisible(False)

        selected_items = list(self.scene.selectedItems())
        saved_anchor = self.bond_anchor_id
        hovered_atom_id = self.hovered_atom_id
        hovered_bond_id = self.hovered_bond_id

        if hovered_atom_id in self.atom_items:
            self.atom_items[hovered_atom_id].set_hover(False)
        self.hovered_atom_id = None
        self.hovered_bond_id = None

        if selected_items or saved_anchor is not None:
            self.scene.blockSignals(True)
            self.scene.clearSelection()
            self.bond_anchor_id = None
            self.scene.blockSignals(False)
            self._sync_selection_from_scene()

        try:
            return render_fn()
        finally:
            if selected_items or saved_anchor is not None:
                self.scene.blockSignals(True)
                for item in selected_items:
                    item.setSelected(True)
                self.bond_anchor_id = saved_anchor
                self.scene.blockSignals(False)
                self._sync_selection_from_scene()

            self.hovered_atom_id = hovered_atom_id
            self.hovered_bond_id = hovered_bond_id
            if hovered_atom_id in self.atom_items:
                self.atom_items[hovered_atom_id].set_hover(True)

            for item in hidden:
                item.setVisible(True)

    def _with_hidden_unselected(self, render_fn):
        hidden = []
        selected = set(self.scene.selectedItems())
        for atom_id in self.state.selected_atoms:
            item = self.atom_items.get(atom_id)
            if item is not None:
                selected.add(item)
        for bond_id in self.state.selected_bonds:
            item = self.bond_items.get(bond_id)
            if item is not None:
                selected.add(item)
        for item in self.scene.items():
            parent = item.parentItem()
            if item in selected or (parent is not None and parent in selected):
                continue
            if item.isVisible():
                hidden.append(item)
                item.setVisible(False)
        try:
            return render_fn()
        finally:
            for item in hidden:
                item.setVisible(True)

    def _render_scene_image(
        self, scale: float = 1.0, selected_only: bool = False
    ) -> Optional[QImage]:
        rect = self._render_scene_bounds(selected_only=selected_only)
        if rect is None:
            return None
        scale = max(1.0, float(scale))
        width = max(1, math.ceil(rect.width() * scale))
        height = max(1, math.ceil(rect.height() * scale))

        def render():
            image = QImage(width, height, QImage.Format.Format_ARGB32)
            image.fill(Qt.GlobalColor.transparent)
            painter = QPainter(image)
            painter.setRenderHint(QPainter.RenderHint.Antialiasing)
            painter.setRenderHint(QPainter.RenderHint.TextAntialiasing)
            self.scene.render(painter, QRectF(0, 0, width, height), rect)
            painter.end()
            trimmed = self._trim_transparent_image(image)
            if trimmed is None:
                return None
            self._apply_image_dpi(trimmed, scale)
            return trimmed

        if selected_only:
            return self._with_hidden_unselected(
                lambda: self._with_hidden_render_items(render)
            )
        return self._with_hidden_render_items(render)

    def _apply_image_dpi(self, image: QImage, scale: float) -> None:
        base_dpi_x = self.logicalDpiX() or 96.0
        base_dpi_y = self.logicalDpiY() or 96.0
        dpi_x = base_dpi_x * scale
        dpi_y = base_dpi_y * scale
        dpm_x = max(1, round(dpi_x * 1000.0 / 25.4))
        dpm_y = max(1, round(dpi_y * 1000.0 / 25.4))
        image.setDotsPerMeterX(dpm_x)
        image.setDotsPerMeterY(dpm_y)

    def _trim_transparent_image(self, image: QImage) -> Optional[QImage]:
        width = image.width()
        height = image.height()
        if width <= 0 or height <= 0:
            return None
        left = width
        right = -1
        top = height
        bottom = -1
        for y in range(height):
            for x in range(width):
                if (image.pixel(x, y) >> 24) & 0xFF:
                    if x < left:
                        left = x
                    if x > right:
                        right = x
                    if y < top:
                        top = y
                    if y > bottom:
                        bottom = y
        if right < left or bottom < top:
            return None
        pad = 1
        left = max(0, left - pad)
        top = max(0, top - pad)
        right = min(width - 1, right + pad)
        bottom = min(height - 1, bottom + pad)
        return image.copy(QRect(left, top, right - left + 1, bottom - top + 1))

    def _render_scene_svg(self) -> Optional[bytes]:
        if QSvgGenerator is None:
            return None
        rect = self._render_scene_bounds()
        if rect is None:
            return None
        width = max(1, math.ceil(rect.width()))
        height = max(1, math.ceil(rect.height()))

        def render():
            buffer = QBuffer()
            buffer.open(QBuffer.OpenModeFlag.WriteOnly)
            generator = QSvgGenerator()
            generator.setOutputDevice(buffer)
            generator.setSize(QSize(width, height))
            generator.setViewBox(QRect(0, 0, width, height))
            painter = QPainter(generator)
            self.scene.render(painter, QRectF(0, 0, width, height), rect)
            painter.end()
            return bytes(buffer.data())

        return self._with_hidden_render_items(render)

    def _insert_molgraph(self, graph: MolGraph) -> None:
        if not graph.atoms:
            return
        xs = [atom.x for atom in graph.atoms.values()]
        ys = [atom.y for atom in graph.atoms.values()]
        min_x, max_x = min(xs), max(xs)
        min_y, max_y = min(ys), max(ys)
        center_x = (min_x + max_x) / 2
        center_y = (min_y + max_y) / 2
        target_x = self.paper_width / 2
        target_y = self.paper_height / 2
        dx = target_x - center_x
        dy = target_y - center_y

        self.undo_stack.beginMacro("Paste molecule")
        id_map: Dict[int, int] = {}
        for atom in graph.atoms.values():
            cmd = AddAtomCommand(
                self.model,
                self,
                atom.element,
                atom.x + dx,
                atom.y + dy,
                is_explicit=atom.is_explicit,
                auto_hydrogens=False,
            )
            self.undo_stack.push(cmd)
            if cmd.atom_id is not None:
                id_map[atom.id] = cmd.atom_id
        for bond in graph.bonds.values():
            a1 = id_map.get(bond.a1_id)
            a2 = id_map.get(bond.a2_id)
            if a1 is None or a2 is None:
                continue
            cmd = AddBondCommand(
                self.model,
                self,
                a1,
                a2,
                bond.order,
                bond.style,
                bond.stereo,
                is_aromatic=bond.is_aromatic,
            )
            self.undo_stack.push(cmd)
        self.undo_stack.endMacro()
        if any(bond.is_aromatic for bond in self.model.bonds.values()):
            self._kekulize_aromatic_bonds()

    def _insert_image_from_clipboard(self, mime: QMimeData) -> None:
        pixmap = None
        if mime.hasImage():
            image = mime.imageData()
            if isinstance(image, QImage):
                pixmap = QPixmap.fromImage(image)
            elif isinstance(image, QPixmap):
                pixmap = image
        if pixmap is None and mime.hasFormat("image/png"):
            pixmap = QPixmap()
            pixmap.loadFromData(mime.data("image/png"), "PNG")
        if pixmap is None and mime.hasFormat("image/svg+xml"):
            data = bytes(mime.data("image/svg+xml"))
            pixmap = QPixmap()
            pixmap.loadFromData(data)
        if pixmap is None or pixmap.isNull():
            return
        item = QGraphicsPixmapItem(pixmap)
        item.setZValue(20)
        item.setPos(
            self.paper_width / 2 - pixmap.width() / 2,
            self.paper_height / 2 - pixmap.height() / 2,
        )
        self.scene.addItem(item)

    def add_atom_item(self, atom) -> None:
        if atom.id in self.atom_items:
            return
        # Determine visibility based on state preferences
        show_c = self.state.show_implicit_carbons or atom.element != "C"
        show_h = self.state.show_implicit_hydrogens or atom.element != "H" or atom.is_explicit
        item = AtomItem(
            atom,
            show_carbon=show_c,
            show_hydrogen=show_h,
            label_font=self._label_font(),
            style=self.drawing_style,
            use_element_colors=self.state.use_element_colors,
        )
        self.scene.addItem(item)
        self.atom_items[atom.id] = item
        self._refresh_atom_label(atom.id)

    def update_atom_item(self, atom_id: int, x: float, y: float) -> None:
        item = self.atom_items.get(atom_id)
        if item is not None:
            item.setPos(x, y)

    def update_atom_item_element(
        self,
        atom_id: int,
        element: str,
        is_explicit: Optional[bool] = None,
    ) -> None:
        item = self.atom_items.get(atom_id)
        if item is not None:
            if is_explicit is None:
                is_explicit = self.model.get_atom(atom_id).is_explicit
            item.set_element(element, is_explicit=is_explicit)
            self._refresh_atom_label(atom_id)

    def update_atom_item_charge(self, atom_id: int, charge: int) -> None:
        item = self.atom_items.get(atom_id)
        if item is not None:
            item.set_charge(charge)

    def _label_font(self) -> QFont:
        size = self.state.label_font_size
        if size <= 0:
            size = 10.0
        font = QFont(self.state.label_font_family)
        font.setPointSizeF(size)
        font.setBold(self.state.label_font_bold)
        font.setItalic(self.state.label_font_italic)
        font.setUnderline(self.state.label_font_underline)
        return font

    def label_font(self) -> QFont:
        return self._label_font()

    def apply_label_font(self, font: QFont) -> None:
        size = font.pointSizeF()
        if size <= 0:
            size = font.pointSize()
        if size <= 0:
            size = 10.0
        self.state.label_font_family = font.family()
        self.state.label_font_size = float(size)
        self.state.label_font_bold = font.bold()
        self.state.label_font_italic = font.italic()
        self.state.label_font_underline = font.underline()
        self.refresh_label_fonts()

    def refresh_label_fonts(self) -> None:
        font = self._label_font()
        for item in self.atom_items.values():
            item.set_label_font(font)
        self.refresh_atom_labels()

    def set_use_element_colors(self, use_element_colors: bool) -> None:
        self.state.use_element_colors = use_element_colors
        for item in self.atom_items.values():
            item.set_use_element_colors(use_element_colors)
        self.refresh_atom_labels()

    def refresh_atom_labels(self, atom_ids: Optional[Iterable[int]] = None) -> None:
        if atom_ids is None:
            atom_ids = list(self.atom_items.keys())
        for atom_id in atom_ids:
            self._refresh_atom_label(atom_id)

    def _refresh_atom_label(self, atom_id: int) -> None:
        item = self.atom_items.get(atom_id)
        atom = self.model.atoms.get(atom_id)
        if item is None or atom is None:
            return
        label, anchor, offset = self._build_atom_label(atom)
        item.set_display_label(label, anchor)
        item.set_label_offset(offset)
        self._update_bond_label_shrinks({atom_id})

    def _build_atom_label(self, atom) -> tuple[str, Optional[str], QPointF]:
        label = atom.element
        anchor: Optional[str] = None
        anchor_override = self._group_anchor_overrides.get(atom.id)
        if atom.element in ELEMENT_SYMBOLS:
            if atom.element == "C" and not (
                self.state.show_implicit_carbons or atom.is_explicit
            ):
                return label, None, QPointF(0.0, 0.0)
            implicit_h = self._implicit_hydrogen_count(atom.id, atom.element)
            if implicit_h > 0 and atom.element != "H" and not self.state.show_implicit_hydrogens:
                h_text = "H" if implicit_h == 1 else f"H{implicit_h}"
                if self._prefer_prefix_h(atom.id):
                    label = f"{h_text}{atom.element}"
                    anchor = atom.element
                else:
                    label = f"{atom.element}{h_text}"
        else:
            label, anchor = self._reflow_group_label(label, atom.id, anchor_override)
            if anchor is None and anchor_override and anchor_override in label:
                anchor = anchor_override
        offset = self._label_offset(atom.id)
        return label, anchor, offset

    def _reflow_group_label(
        self, label: str, atom_id: int, anchor_override: Optional[str] = None
    ) -> tuple[str, Optional[str]]:
        cleaned = label.strip()
        if not cleaned:
            return label, None

        charge = ""
        if cleaned and cleaned[-1] in "+-":
            charge = cleaned[-1]
            cleaned = cleaned[:-1]

        tokens = self._group_label_tokens(cleaned)
        # If the user explicitly chose an anchor atom, force it to the front.
        # This keeps labels like "TBSO" consistent with an O-anchored bond.
        if anchor_override:
            anchor_first = True
        else:
            direction = self._label_open_direction(atom_id)
            anchor_first = direction.x() >= -0.2

        if anchor_override and anchor_override in cleaned:
            if not tokens or len(tokens) < 2:
                # Use bond direction to decide if anchor goes first or last
                cleaned = self._move_anchor_in_label(cleaned, anchor_override, anchor_first)
                return f"{cleaned}{charge}", anchor_override
            if not any(symbol == anchor_override for symbol, _token in tokens):
                # Anchor letter is inside an abbreviation token (e.g., TBSO -> B/S).
                cleaned = self._move_anchor_in_label(cleaned, anchor_override, anchor_first)
                return f"{cleaned}{charge}", anchor_override
        if not tokens or len(tokens) < 2:
            return label, None

        if anchor_override and any(symbol == anchor_override for symbol, _token in tokens):
            if anchor_first:
                while tokens[0][0] != anchor_override:
                    tokens = tokens[1:] + tokens[:1]
            else:
                while tokens[-1][0] != anchor_override:
                    tokens = tokens[1:] + tokens[:1]
            anchor_symbol = anchor_override
        else:
            anchor_symbol = tokens[0][0]
            if not anchor_first:
                # Move first token to end
                tokens = tokens[1:] + tokens[:1]
                anchor_symbol = tokens[-1][0]
        cleaned = "".join(token for _symbol, token in tokens)
        return f"{cleaned}{charge}", anchor_symbol

    def _move_anchor_in_label(self, label: str, anchor: str, anchor_first: bool) -> str:
        if not label or not anchor:
            return label
        if anchor_first:
            idx = label.find(anchor)
        else:
            idx = label.rfind(anchor)
        if idx < 0:
            return label
        before = label[:idx]
        after = label[idx + len(anchor):]
        core = before + after
        return f"{anchor}{core}" if anchor_first else f"{core}{anchor}"

    def get_anchor_override(self, atom_id: int) -> Optional[str]:
        return self._group_anchor_overrides.get(atom_id)

    def set_anchor_override(self, atom_id: int, anchor: Optional[str]) -> None:
        if anchor:
            self._group_anchor_overrides[atom_id] = anchor
        else:
            self._group_anchor_overrides.pop(atom_id, None)

    def _group_label_tokens(self, label: str) -> Optional[list[tuple[str, str]]]:
        tokens: list[tuple[str, str]] = []
        i = 0
        sorted_abbrs = sorted(list(ABBREVIATION_LABELS), key=len, reverse=True)
        
        while i < len(label):
            ch = label[i]
            symbol = None
            
            # 1. Try abbreviations first
            for abbr in sorted_abbrs:
                if label.startswith(abbr, i):
                    symbol = abbr
                    i += len(abbr)
                    break
            
            # 2. Try element symbols if no abbreviation matched
            if symbol is None:
                if not ch.isupper():
                    return None
                if i + 1 < len(label) and label[i + 1].islower():
                    candidate = label[i : i + 2]
                    if candidate in ELEMENT_SYMBOLS:
                        symbol = candidate
                        i += 2
                if symbol is None and ch in ELEMENT_SYMBOLS:
                    symbol = ch
                    i += 1
            
            if symbol is None:
                return None
                
            token = symbol
            # Collect numbers and script markers
            while i < len(label):
                if label[i].isdigit():
                    token += label[i]
                    i += 1
                    continue
                if label[i] in "_^":
                    marker = label[i]
                    token += marker
                    i += 1
                    while i < len(label) and label[i].isalnum():
                        token += label[i]
                        i += 1
                    continue
                break
            tokens.append((symbol, token))
            
        return tokens if tokens else None

    def _implicit_hydrogen_count(self, atom_id: int, element: str) -> int:
        if element not in IMPLICIT_H_ELEMENTS:
            return 0
        expected = VALENCE_MAP.get(element)
        if expected is None:
            return 0
        bond_order = 0
        for bond in self.model.bonds.values():
            if bond.a1_id != atom_id and bond.a2_id != atom_id:
                continue
            order = bond.display_order if bond.display_order is not None else bond.order
            bond_order += order
        return max(0, expected - bond_order)

    def _label_open_direction(self, atom_id: int) -> QPointF:
        def collect_vectors(skip_hydrogen: bool) -> list[tuple[float, float]]:
            atom = self.model.get_atom(atom_id)
            vectors: list[tuple[float, float]] = []
            for bond in self.model.bonds.values():
                if bond.a1_id == atom_id:
                    other = self.model.get_atom(bond.a2_id)
                elif bond.a2_id == atom_id:
                    other = self.model.get_atom(bond.a1_id)
                else:
                    continue
                if skip_hydrogen and other.element == "H":
                    continue
                dx = other.x - atom.x
                dy = other.y - atom.y
                length = math.hypot(dx, dy)
                if length <= 1e-6:
                    continue
                vectors.append((dx / length, dy / length))
            return vectors

        vectors = collect_vectors(skip_hydrogen=True)
        if not vectors:
            vectors = collect_vectors(skip_hydrogen=False)
        if not vectors:
            return QPointF(0.0, 0.0)
        sum_x = sum(v[0] for v in vectors)
        sum_y = sum(v[1] for v in vectors)
        if abs(sum_x) + abs(sum_y) < 1e-3:
            if len(vectors) == 1:
                sum_x, sum_y = -vectors[0][0], -vectors[0][1]
            else:
                return QPointF(0.0, 0.0)
        else:
            sum_x, sum_y = -sum_x, -sum_y
        length = math.hypot(sum_x, sum_y)
        if length <= 1e-6:
            return QPointF(0.0, 0.0)
        return QPointF(sum_x / length, sum_y / length)

    def _label_offset(self, atom_id: int) -> QPointF:
        atom = self.model.get_atom(atom_id)
        if any(ch.isalpha() for ch in atom.element):
            return QPointF(0.0, 0.0)
        direction = self._label_open_direction(atom_id)
        if direction == QPointF(0.0, 0.0):
            return QPointF(0.0, 0.0)
        size = self._label_font().pointSizeF()
        if size <= 0:
            size = 10.0
        offset = max(LABEL_OFFSET_MIN_PX, size * LABEL_OFFSET_SCALE)
        return QPointF(direction.x() * offset, direction.y() * offset)

    def _neighbor_angles_deg(self, atom_id: int) -> list[float]:
        atom = self.model.get_atom(atom_id)
        origin = QPointF(atom.x, atom.y)
        angles: list[float] = []
        for bond in self.model.bonds.values():
            if bond.a1_id == atom_id:
                other = self.model.get_atom(bond.a2_id)
            elif bond.a2_id == atom_id:
                other = self.model.get_atom(bond.a1_id)
            else:
                continue
            if other.element == "H":
                continue
            angles.append(angle_deg(origin, QPointF(other.x, other.y)))
        return angles

    def _atom_hybridization(self, atom_id: int) -> str:
        has_double = False
        has_triple = False
        has_aromatic = False
        for bond in self.model.bonds.values():
            if bond.a1_id != atom_id and bond.a2_id != atom_id:
                continue
            if bond.is_aromatic:
                has_aromatic = True
            order = bond.display_order if bond.display_order is not None else bond.order
            if order >= 3:
                has_triple = True
            elif order == 2:
                has_double = True
        if has_triple:
            return "sp"
        if has_double or has_aromatic:
            return "sp2"
        return "sp3"

    def _angle_of_vector(self, vx: float, vy: float) -> float:
        # Invert Y to convert from screen coordinates (Y down) to math coordinates (Y up)
        # This matches the convention used in geom.angle_deg and endpoint_from_angle_len
        return math.degrees(math.atan2(-vy, vx))

    def _normalize_angle(self, angle: float) -> float:
        value = angle % 360.0
        return value + 360.0 if value < 0 else value

    def _implicit_h_angles(self, atom_id: int, count: int) -> list[float]:
        """
        Calculate angles for implicit hydrogen atoms to be drawn.
        
        This function finds the best positions for implicit H atoms by:
        1. Finding all gaps between existing bonds
        2. Distributing H atoms into the largest gaps
        """
        if count <= 0:
            return []
        
        atom = self.model.get_atom(atom_id)
        origin = QPointF(atom.x, atom.y)
        
        # Collect all existing bond angles using angle_deg from geom.py
        # This function correctly handles screen coordinates (Y down -> math Y up)
        # Skip H atoms since we're computing where to place implicit H
        neighbor_angles: list[float] = []
        for bond in self.model.bonds.values():
            if bond.a1_id == atom_id:
                other = self.model.get_atom(bond.a2_id)
            elif bond.a2_id == atom_id:
                other = self.model.get_atom(bond.a1_id)
            else:
                continue
            # Skip existing hydrogen atoms
            if other.element == "H":
                continue
            # Use angle_deg from geom.py for consistent coordinate handling
            a = angle_deg(origin, QPointF(other.x, other.y))
            neighbor_angles.append(a)
        
        # No existing bonds - use default positions
        if not neighbor_angles:
            return [0.0, 120.0, 240.0, 300.0][:count]
        
        # Find gaps and place H atoms in the largest gaps
        return self._find_best_h_positions(neighbor_angles, count)
    
    def _find_best_h_positions(self, occupied_angles: list[float], count: int) -> list[float]:
        """
        Find the best positions for implicit H atoms using a bisector strategy.
        
        Args:
            occupied_angles: List of angles (in degrees) where bonds already exist
            count: Number of H atoms to place
            
        Returns:
            List of angles where H atoms should be placed
        """
        if count <= 0:
            return []
        
        if not occupied_angles:
            # Default positions (triangular/tetrahedral-ish)
            if count == 1: return [0.0]
            if count == 2: return [120.0, 240.0]
            return [0.0, 120.0, 240.0]
        
        # Local function to normalize angle to [0, 360) degrees
        def norm_deg(a: float) -> float:
            v = a % 360.0
            return v + 360.0 if v < 0 else v
        
        # Sort angles
        sorted_angles = sorted(norm_deg(a) for a in occupied_angles)
        n = len(sorted_angles)
        
        # Find the single largest gap
        max_gap = 0.0
        best_mid = 0.0
        
        for i in range(n):
            start = sorted_angles[i]
            end = sorted_angles[(i + 1) % n]
            
            if i == n - 1:
                gap_size = (360.0 - start) + end
            else:
                gap_size = end - start
            
            if gap_size > max_gap:
                max_gap = gap_size
                # Bisector of the gap
                best_mid = norm_deg(start + gap_size / 2.0)
        
        result_angles: list[float] = []
        
        # Distribute H atoms symmetrically around the bisector of the largest gap
        if count == 1:
            # Place directly in the middle
            result_angles.append(best_mid)
        elif count == 2:
            # Place in a "V" shape (±30 degrees from bisector)
            # This looks good for sp3 chains (tetrahedral projection)
            spread = 30.0
            if max_gap < 60.0: # If gap is very tight, reduce spread
                spread = max_gap / 3.0
            result_angles.append(norm_deg(best_mid - spread))
            result_angles.append(norm_deg(best_mid + spread))
        else:
            # Count >= 3 (e.g. terminal CH3)
            # Distribute in a fan shape (0, ±60 degrees) or similar
            # Standard "Trident" look: center, +60, -60
            result_angles.append(best_mid)
            result_angles.append(norm_deg(best_mid - 60.0))
            result_angles.append(norm_deg(best_mid + 60.0))
            
            # If we needed more than 3, we'd add more, but C is usually max 4 bonds
            if count > 3:
                # Fallback purely even distribution if strange valency
                start_fan = best_mid - (count - 1) * 30.0
                return [norm_deg(start_fan + i * 60.0) for i in range(count)]

        return result_angles

    def _clear_implicit_h_overlays(self, atom_ids: Optional[Iterable[int]] = None) -> None:
        if atom_ids is None:
            atom_ids = list(self._implicit_h_overlays.keys())
        for atom_id in list(atom_ids):
            overlays = self._implicit_h_overlays.pop(atom_id, [])
            for line_item, text_item in overlays:
                if line_item.scene() is self.scene:
                    self.scene.removeItem(line_item)
                if text_item.scene() is self.scene:
                    self.scene.removeItem(text_item)

    def _refresh_implicit_h_overlays(self, atom_ids: Optional[Iterable[int]] = None) -> None:
        if atom_ids is None:
            atom_ids = list(self.atom_items.keys())
        self._clear_implicit_h_overlays(atom_ids)
        if not self.state.show_implicit_hydrogens:
            return
        font = self._label_font()
        for atom_id in atom_ids:
            atom = self.model.atoms.get(atom_id)
            if atom is None or atom.element == "H":
                continue
            count = self._implicit_hydrogen_count(atom_id, atom.element)
            if count <= 0:
                continue
            bond_length = max(1.0, self._local_bond_length(atom_id))
            angles = self._implicit_h_angles(atom_id, count)
            overlays: list[tuple[QGraphicsLineItem, QGraphicsTextItem]] = []
            for angle_value in angles:
                # Calculate full position for text center
                pos = endpoint_from_angle_len(QPointF(0.0, 0.0), angle_value, bond_length)
                rad = math.radians(angle_value)
                ux = math.cos(rad)
                uy = -math.sin(rad)
                
                # Setup text item first to get its size
                text_item = QGraphicsTextItem("H")
                text_item.setFont(font)
                text_item.setDefaultTextColor(QColor("#000000"))
                text_item.setZValue(10)
                text_item.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
                text_item.setParentItem(self.atom_items[atom_id])
                
                rect = text_item.boundingRect()
                text_item.setPos(pos.x() - rect.width() / 2, pos.y() - rect.height() / 2)
                
                # Calculate shortened line end to avoid overlap with text
                # Shrink by roughly half the text width plus a small padding
                shrink = rect.width() / 2.0 + 3.0
                line_end_len = max(0.0, bond_length - shrink)

                # Also shrink the line start to avoid crossing atom labels (e.g., visible carbons).
                start_shrink = self._label_shrink_for_atom(atom_id, ux, uy)
                line_start_len = max(0.0, start_shrink)
                if line_end_len < line_start_len:
                    line_end_len = line_start_len
                line_start = endpoint_from_angle_len(QPointF(0.0, 0.0), angle_value, line_start_len)
                line_end = endpoint_from_angle_len(QPointF(0.0, 0.0), angle_value, line_end_len)
                
                line_item = QGraphicsLineItem(
                    line_start.x(),
                    line_start.y(),
                    line_end.x(),
                    line_end.y(),
                )
                pen = QPen(QColor(self.drawing_style.bond_color), self.drawing_style.stroke_px)
                pen.setCapStyle(self.drawing_style.cap_style)
                pen.setJoinStyle(self.drawing_style.join_style)
                line_item.setPen(pen)
                line_item.setZValue(-5)
                line_item.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
                line_item.setParentItem(self.atom_items[atom_id])

                overlays.append((line_item, text_item))
            if overlays:
                self._implicit_h_overlays[atom_id] = overlays

    def _local_bond_length(self, atom_id: int) -> float:
        atom = self.model.get_atom(atom_id)
        if atom is None:
            return float(self.state.bond_length)
        lengths: list[float] = []
        for bond in self.model.bonds.values():
            if bond.a1_id == atom_id:
                other_id = bond.a2_id
            elif bond.a2_id == atom_id:
                other_id = bond.a1_id
            else:
                continue
            other = self.model.get_atom(other_id)
            if other is None:
                continue
            length = math.hypot(other.x - atom.x, other.y - atom.y)
            if length > 1e-6:
                lengths.append(length)
        if lengths:
            return sum(lengths) / len(lengths)
        return float(self.state.bond_length)

    def _update_bond_label_shrinks(self, atom_ids: set[int]) -> None:
        for bond in self.model.bonds.values():
            if bond.a1_id not in atom_ids and bond.a2_id not in atom_ids:
                continue
            item = self.bond_items.get(bond.id)
            if item is None:
                continue
            self._configure_bond_rendering(bond, item)
            atom1 = self.model.get_atom(bond.a1_id)
            atom2 = self.model.get_atom(bond.a2_id)
            dx = atom2.x - atom1.x
            dy = atom2.y - atom1.y
            length = math.hypot(dx, dy)
            if length <= 1e-6:
                item.set_label_shrink(0.0, 0.0)
                item.update_positions(atom1, atom2)
                continue
            ux = dx / length
            uy = dy / length
            shrink_start = self._label_shrink_for_atom(bond.a1_id, ux, uy)
            shrink_end = self._label_shrink_for_atom(bond.a2_id, -ux, -uy)
            item.set_label_shrink(shrink_start, shrink_end)
            item.update_positions(atom1, atom2)
        self._refresh_implicit_h_overlays(atom_ids)

    def _label_shrink_for_atom(self, atom_id: int, ux: float, uy: float) -> float:
        item = self.atom_items.get(atom_id)
        if item is None or not item.label.isVisible():
            return 0.0
        rect = item.label.mapRectToParent(item.label.boundingRect())
        if rect.isNull() or rect.width() <= 0 or rect.height() <= 0:
            return 0.0
        # Increased padding to clear label characters (especially "C")
        # 6.0px provides a comfortable margin for standard font sizes.
        pad = 6.0
        atom = self.model.get_atom(atom_id)
        if atom is not None and atom.element in ELEMENT_SYMBOLS and atom.element not in {"C", "H"}:
            pad = 3.5
        rect = rect.adjusted(-pad, -pad, pad, pad)
        distance = self._ray_ellipse_distance(rect, ux, uy)
        return distance if distance is not None else 0.0

    def _configure_bond_rendering(self, bond: Bond, item: BondItem) -> None:
        prefer_full_length = self.state.show_implicit_carbons and self.state.show_implicit_hydrogens
        is_aromatic_ring = bond.is_aromatic and bond.ring_id is not None
        symmetric_double = bond.ring_id is None and not is_aromatic_ring
        item.set_multibond_rendering(prefer_full_length, symmetric_double)

    @staticmethod
    def _ray_ellipse_distance(rect: QRectF, ux: float, uy: float) -> Optional[float]:
        # Ray P = t * D, where D = (ux, uy). Origin (0,0).
        # Ellipse defined by rect: Center C, radii a, b.
        # Equation: ((x - cx)/a)^2 + ((y - cy)/b)^2 = 1
        # Substitute x = t*ux, y = t*uy
        
        a = rect.width() / 2.0
        b = rect.height() / 2.0
        if a <= 1e-9 or b <= 1e-9:
            return None
            
        cx = rect.center().x()
        cy = rect.center().y()
        
        # Coefficients for At^2 + Bt + K = 0
        # (1/a^2)*(t*ux - cx)^2 + (1/b^2)*(t*uy - cy)^2 - 1 = 0
        
        inv_a2 = 1.0 / (a * a)
        inv_b2 = 1.0 / (b * b)
        
        A = inv_a2 * ux * ux + inv_b2 * uy * uy
        B = -2.0 * (inv_a2 * ux * cx + inv_b2 * uy * cy)
        K = inv_a2 * cx * cx + inv_b2 * cy * cy - 1.0
        
        # Discriminant
        disc = B * B - 4 * A * K
        if disc < 0:
            return None
            
        sqrt_disc = math.sqrt(disc)
        
        # Solutions
        t1 = (-B - sqrt_disc) / (2 * A)
        t2 = (-B + sqrt_disc) / (2 * A)
        
        # We want the smallest positive t that is > 0
        # Since ray starts at 0,0 which is usually inside, one t is neg, one is pos.
        # We want the positive one.
        
        best = None
        if t1 > 1e-6:
            best = t1
        if t2 > 1e-6:
            if best is None or t2 < best:
                best = t2
                
        return best

    def _prefer_prefix_h(self, atom_id: int) -> bool:
        direction = self._label_open_direction(atom_id)
        return direction.x() < -0.2

    def remove_atom_item(self, atom_id: int) -> None:
        item = self.atom_items.pop(atom_id, None)
        if item is not None:
            self.scene.removeItem(item)
        self._clear_implicit_h_overlays([atom_id])
        self._group_anchor_overrides.pop(atom_id, None)
        if self.bond_anchor_id == atom_id:
            self.bond_anchor_id = None
        if self.hovered_atom_id == atom_id:
            self.hovered_atom_id = None

    def add_bond_item(self, bond) -> None:
        if bond.id in self.bond_items:
            return
        atom1 = self.model.get_atom(bond.a1_id)
        atom2 = self.model.get_atom(bond.a2_id)
        item = BondItem(
            bond,
            atom1,
            atom2,
            render_aromatic_as_circle=self.state.use_aromatic_circles,
            style=self.drawing_style,
        )
        if bond.ring_id is not None:
            item.set_ring_context(self._ring_centers.get(bond.ring_id))
        item.set_offset_sign(self._bond_offset_sign(bond))
        self._configure_bond_rendering(bond, item)
        item.update_positions(atom1, atom2)
        self.scene.addItem(item)
        self.bond_items[bond.id] = item
        self._refresh_atom_label(bond.a1_id)
        self._refresh_atom_label(bond.a2_id)

    def add_arrow_item(self, start: QPointF, end: QPointF, kind: str) -> ArrowItem:
        item = ArrowItem(start, end, kind=kind, style=self.drawing_style)
        self.scene.addItem(item)
        self.arrow_items.append(item)
        return item

    def readd_arrow_item(self, item: ArrowItem, start: QPointF, end: QPointF, kind: str) -> None:
        item.set_kind(kind)
        item.update_positions(start, end)
        if item.scene() is not self.scene:
            self.scene.addItem(item)
        if item not in self.arrow_items:
            self.arrow_items.append(item)

    def remove_arrow_item(self, item: ArrowItem) -> None:
        if item in self.arrow_items:
            self.arrow_items.remove(item)
        if item.scene() is self.scene:
            self.scene.removeItem(item)

    def add_bracket_item(self, rect: QRectF, kind: str) -> BracketItem:
        item = BracketItem(rect, kind=kind, style=self.drawing_style)
        self.scene.addItem(item)
        self.bracket_items.append(item)
        return item

    def readd_bracket_item(
        self,
        item: BracketItem,
        rect: QRectF,
        kind: str,
        padding: Optional[float] = None,
    ) -> None:
        item.set_rect(rect, padding=padding)
        item._kind = kind
        if item.scene() is not self.scene:
            self.scene.addItem(item)
        if item not in self.bracket_items:
            self.bracket_items.append(item)

    def remove_bracket_item(self, item: BracketItem) -> None:
        if item in self.bracket_items:
            self.bracket_items.remove(item)
        if item.scene() is self.scene:
            self.scene.removeItem(item)

    def _ensure_bracket_preview(self) -> QGraphicsRectItem:
        if self._bracket_preview is None:
            preview = QGraphicsRectItem()
            pen = QPen(QColor("#4A90D9"), 1.0, Qt.PenStyle.DashLine)
            preview.setPen(pen)
            preview.setBrush(QBrush(Qt.BrushStyle.NoBrush))
            preview.setZValue(40)
            self.scene.addItem(preview)
            self._bracket_preview = preview
        self._bracket_preview.setVisible(True)
        return self._bracket_preview

    def _clear_bracket_preview(self) -> None:
        if self._bracket_preview is not None:
            self.scene.removeItem(self._bracket_preview)
            self._bracket_preview = None
        self._bracket_drag_start = None

    def update_bond_item(self, bond_id: int) -> None:
        bond = self.model.get_bond(bond_id)
        atom1 = self.model.get_atom(bond.a1_id)
        atom2 = self.model.get_atom(bond.a2_id)
        item = self.bond_items.get(bond_id)
        if item is not None:
            if bond.ring_id is not None:
                item.set_ring_context(self._ring_centers.get(bond.ring_id))
            else:
                item.set_ring_context(None)
            item.set_offset_sign(self._bond_offset_sign(bond))
            self._configure_bond_rendering(bond, item)
            item.set_bond(bond, atom1, atom2)
            item.set_style(self.drawing_style, atom1, atom2)
        self._refresh_atom_label(bond.a1_id)
        self._refresh_atom_label(bond.a2_id)

    def update_bond_items_for_atoms(self, atom_ids: set[int]) -> None:
        for bond in self.model.bonds.values():
            if bond.a1_id in atom_ids or bond.a2_id in atom_ids:
                self.update_bond_item(bond.id)
        self._refresh_implicit_h_overlays(atom_ids)

    def remove_bond_item(self, bond_id: int) -> None:
        item = self.bond_items.pop(bond_id, None)
        if item is not None:
            a1_id = item.a1_id
            a2_id = item.a2_id
            self.scene.removeItem(item)
            self._refresh_atom_label(a1_id)
            self._refresh_atom_label(a2_id)

    def allocate_ring_id(self) -> int:
        ring_id = self._next_ring_id
        self._next_ring_id += 1
        return ring_id

    def register_ring_center(self, ring_id: int, center: tuple[float, float]) -> None:
        self._ring_centers[ring_id] = QPointF(center[0], center[1])

    def unregister_ring_center(self, ring_id: int) -> None:
        self._ring_centers.pop(ring_id, None)

    def _bond_offset_sign(self, bond) -> int:
        atom1 = self.model.get_atom(bond.a1_id)
        atom2 = self.model.get_atom(bond.a2_id)
        p1 = QPointF(atom1.x, atom1.y)
        p2 = QPointF(atom2.x, atom2.y)
        mid = QPointF((p1.x() + p2.x()) / 2, (p1.y() + p2.y()) / 2)
        dx = p2.x() - p1.x()
        dy = p2.y() - p1.y()
        length = math.hypot(dx, dy) or 1.0
        nx = -dy / length
        ny = dx / length

        def neighbor_score(atom_id: int) -> float:
            score = 0.0
            for other_bond in self.model.bonds.values():
                if other_bond.id == bond.id:
                    continue
                if other_bond.a1_id == atom_id:
                    other = self.model.get_atom(other_bond.a2_id)
                elif other_bond.a2_id == atom_id:
                    other = self.model.get_atom(other_bond.a1_id)
                else:
                    continue
                vx = other.x - mid.x()
                vy = other.y - mid.y()
                score += nx * vx + ny * vy
            return score

        score = neighbor_score(atom1.id) + neighbor_score(atom2.id)
        return -1 if score > 0 else 1

    # -------------------------------------------------------------------------
    # Visibility and Aromatic Circle Management
    # -------------------------------------------------------------------------
    def refresh_atom_visibility(self) -> None:
        """Update visibility of all atoms based on current state settings."""
        for atom_id, item in self.atom_items.items():
            atom = self.model.get_atom(atom_id)
            show_c = self.state.show_implicit_carbons or atom.element != "C"
            show_h = self.state.show_implicit_hydrogens or atom.element != "H" or atom.is_explicit
            item.set_visibility_flags(show_c, show_h)
            self._refresh_atom_label(atom_id)
        self._refresh_implicit_h_overlays()

    def _apply_text_settings(self, item: TextAnnotationItem, property_name: str = "all") -> None:
        """Apply text format settings to the item, ideally merging only what changed."""
        s = self._current_text_settings
        
        # Prepare a partial format for merging
        fmt = QTextCharFormat()
        
        if property_name in ("all", "family"):
            fmt.setFontFamilies([s["family"]])
        if property_name in ("all", "size"):
            fmt.setFontPointSize(s["size"])
        if property_name in ("all", "bold"):
            fmt.setFontWeight(QFont.Weight.Bold if s["bold"] else QFont.Weight.Normal)
        if property_name in ("all", "italic"):
            fmt.setFontItalic(s["italic"])
        if property_name in ("all", "underline"):
            fmt.setFontUnderline(s["underline"])
        if property_name in ("all", "sub", "sup"):
            if s["sub"]:
                fmt.setVerticalAlignment(QTextCharFormat.VerticalAlignment.AlignSubScript)
            elif s["sup"]:
                fmt.setVerticalAlignment(QTextCharFormat.VerticalAlignment.AlignSuperScript)
            else:
                fmt.setVerticalAlignment(QTextCharFormat.VerticalAlignment.AlignNormal)
        
        if property_name in ("all", "color"):
            fmt.setForeground(QBrush(s["color"]))
        
        # Get Cursor
        cursor = item.textCursor()
        is_editing = (item.textInteractionFlags() & Qt.TextInteractionFlag.TextEditorInteraction)
        
        if not is_editing:
            # Not editing -> Apply to whole document
            cursor.select(cursor.SelectionType.Document)
            # For whole document, also update default properties
            if property_name in ("all", "family", "size", "bold", "italic", "underline"):
                font = item.font()
                if property_name in ("all", "family"): font.setFamily(s["family"])
                if property_name in ("all", "size"): font.setPointSize(s["size"])
                if property_name in ("all", "bold"): font.setBold(s["bold"])
                if property_name in ("all", "italic"): font.setItalic(s["italic"])
                if property_name in ("all", "underline"): font.setUnderline(s["underline"])
                item.setFont(font)
            if property_name in ("all", "color"):
                item.setDefaultTextColor(s["color"])
        
        # Apply the format
        cursor.mergeCharFormat(fmt)
        item.setTextCursor(cursor)

    def update_text_format(self, family: str, size: int, b: bool, i: bool, u: bool, sub: bool, sup: bool, property_name: str = "all") -> None:
        """Update current text settings and apply to selected text items."""
        self._current_text_settings.update({
            "family": family, "size": size,
            "bold": b, "italic": i, "underline": u,
            "sub": sub, "sup": sup
        })
        
        # Items to update
        items_to_update = self.scene.selectedItems()
        
        # If a text item is in edit mode, it might assume it's handling input.
        # Ensure we capture it.
        focus_item = self.scene.focusItem()
        if isinstance(focus_item, TextAnnotationItem) and focus_item not in items_to_update:
             items_to_update.append(focus_item)

        for item in items_to_update:
            if isinstance(item, TextAnnotationItem):
                self._apply_text_settings(item, property_name)
            elif isinstance(item, AtomItem):
                pass

    def update_text_color(self, color: QColor) -> None:
        """Update current text color and apply to selection."""
        self._current_text_settings["color"] = color
        
        items_to_update = self.scene.selectedItems()
        focus_item = self.scene.focusItem()
        if isinstance(focus_item, TextAnnotationItem) and focus_item not in items_to_update:
            items_to_update.append(focus_item)

        for item in items_to_update:
            if isinstance(item, TextAnnotationItem):
                self._apply_text_settings(item, "color")

    def refresh_aromatic_circles(self) -> None:
        """Add/remove aromatic circle items based on current state."""
        # First, update all bond items to reflect the preference
        use_circles = self.state.use_aromatic_circles
        for bond_id, item in self.bond_items.items():
            if item.is_aromatic:
                item.set_render_aromatic_as_circle(use_circles)
                # Force redraw
                self.update_bond_item(bond_id)

        # Remove existing circles
        for circle in self.aromatic_circles:
            self.scene.removeItem(circle)
        self.aromatic_circles.clear()
        
        if not use_circles:
            return
        
        # Find aromatic rings and add circles
        rings = self._detect_aromatic_rings()
        for center_x, center_y, radius in rings:
            circle = AromaticCircleItem(center_x, center_y, radius)
            self.scene.addItem(circle)
            self.aromatic_circles.append(circle)

    def _detect_aromatic_rings(self) -> list:
        """
        Detect complete aromatic rings and return their centers.
        Returns list of (center_x, center_y, radius) tuples.
        """
        # Collect aromatic bonds
        aromatic_bonds = [b for b in self.model.bonds.values() if b.is_aromatic]
        if not aromatic_bonds:
            return []
        
        # Build adjacency for aromatic atoms only
        aromatic_atoms = set()
        adjacency: dict[int, list[int]] = {}
        for bond in aromatic_bonds:
            aromatic_atoms.add(bond.a1_id)
            aromatic_atoms.add(bond.a2_id)
            adjacency.setdefault(bond.a1_id, []).append(bond.a2_id)
            adjacency.setdefault(bond.a2_id, []).append(bond.a1_id)
        
        # Find rings using simple DFS for small rings (5-7 members)
        visited_rings = set()
        rings = []
        
        for start in aromatic_atoms:
            # Try to find rings starting from this atom
            ring = self._find_ring_from(start, adjacency, max_size=7)
            if ring:
                ring_key = tuple(sorted(ring))
                if ring_key not in visited_rings:
                    visited_rings.add(ring_key)
                    rings.append(ring)
        
        # Calculate centers and radii for each ring
        result = []
        for ring in rings:
            atoms = [self.model.get_atom(aid) for aid in ring]
            cx = sum(a.x for a in atoms) / len(atoms)
            cy = sum(a.y for a in atoms) / len(atoms)
            # Radius is about 60% of distance from center to atoms
            avg_dist = sum(((a.x - cx)**2 + (a.y - cy)**2)**0.5 for a in atoms) / len(atoms)
            radius = avg_dist * 0.55
            result.append((cx, cy, radius))
        
        return result

    def apply_drawing_style(self, style: DrawingStyle) -> None:
        """Apply a new drawing style and refresh items."""
        self.drawing_style = style
        self.state.bond_length = style.bond_length_px
        for item in self.atom_items.values():
            item.set_style(style)
        self.update_bond_items_for_atoms(set(self.model.atoms.keys()))
        for arrow in self.arrow_items:
            arrow.set_style(style)
        for bracket in self.bracket_items:
            bracket.set_style(style)
        self.refresh_atom_visibility()

    def _find_ring_from(self, start: int, adjacency: dict, max_size: int) -> Optional[list]:
        """Find a ring starting from the given atom using DFS."""
        # Simple BFS to find shortest cycle
        from collections import deque
        
        queue = deque([(start, [start])])
        visited = {start}
        
        while queue:
            current, path = queue.popleft()
            if len(path) > max_size:
                continue
            
            for neighbor in adjacency.get(current, []):
                if neighbor == start and len(path) >= 5:
                    return path
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append((neighbor, path + [neighbor]))
        
        return None

    def _get_atom_at(self, x: float, y: float) -> Optional[int]:
        for atom in self.model.atoms.values():
            dx = atom.x - x
            dy = atom.y - y
            if (dx * dx + dy * dy) < (ATOM_HIT_RADIUS * ATOM_HIT_RADIUS):
                return atom.id
        return None

    def _get_bond_at(self, scene_pos: QPointF) -> Optional[int]:
        for item in self.scene.items(scene_pos):
            if isinstance(item, BondItem):
                return item.bond_id
        return None



    def _sync_selection_from_scene(self) -> None:
        try:
            selected_items = list(self.scene.selectedItems())
        except RuntimeError:
            return
        selected_atoms = set()
        selected_bonds = set()
        for item in selected_items:
            if isinstance(item, AtomItem):
                if item.atom_id in self.model.atoms:
                    selected_atoms.add(item.atom_id)
            elif isinstance(item, BondItem):
                if item.bond_id in self.model.bonds:
                    selected_bonds.add(item.bond_id)

        selected_text_items = [item for item in selected_items if isinstance(item, TextAnnotationItem)]
        
        self.state.selected_atoms = selected_atoms
        self.state.selected_bonds = selected_bonds

        for atom_id, item in self.atom_items.items():
            item.set_selected(atom_id in selected_atoms or atom_id == self.bond_anchor_id)

        self._update_selection_overlay()

        # Emit selection signal
        details = {}
        if len(selected_atoms) == 1 and not selected_bonds and not selected_text_items:
            atom = self.model.get_atom(next(iter(selected_atoms)))
            details = {"type": "atom", "id": atom.id, "element": atom.element, "charge": atom.charge, "x": atom.x, "y": atom.y}
        elif len(selected_bonds) == 1 and not selected_atoms and not selected_text_items:
            bond_id = next(iter(selected_bonds))
            if bond_id in self.model.bonds:
                bond = self.model.get_bond(bond_id)
                details = {"type": "bond", "id": bond.id, "order": bond.order, "style": bond.style, "aromatic": bond.is_aromatic}
        elif len(selected_text_items) == 1 and not selected_atoms and not selected_bonds:
            item = selected_text_items[0]
            cursor = item.textCursor()
            fmt = cursor.charFormat()
            details = {
                "type": "text",
                "font": item.font(),
                "color": item.defaultTextColor(),
                "sub": fmt.verticalAlignment() == QTextCharFormat.VerticalAlignment.AlignSubScript,
                "sup": fmt.verticalAlignment() == QTextCharFormat.VerticalAlignment.AlignSuperScript
            }
            
        self.selection_changed.emit(len(selected_atoms), len(selected_bonds), len(selected_text_items), details)

    def _begin_selection_drag(
        self,
        start_pos: QPointF,
        free_select: bool,
        additive: bool,
    ) -> None:
        self._select_drag_mode = "free" if free_select else "rect"
        self._select_start_pos = start_pos
        self._select_additive = additive
        if free_select:
            path = QPainterPath(start_pos)
            self._select_path = path
            if self._select_preview_path is None:
                item = QGraphicsPathItem()
                pen = QPen(QColor("#4A90D9"), 1.2, Qt.PenStyle.DashLine)
                item.setPen(pen)
                item.setBrush(QBrush(Qt.BrushStyle.NoBrush))
                item.setZValue(40)
                self.scene.addItem(item)
                self._select_preview_path = item
            self._select_preview_path.setPath(path)
            self._select_preview_path.setVisible(True)
        else:
            if self._select_preview_rect is None:
                rect_item = QGraphicsRectItem()
                pen = QPen(QColor("#4A90D9"), 1.2, Qt.PenStyle.DashLine)
                rect_item.setPen(pen)
                rect_item.setBrush(QBrush(Qt.BrushStyle.NoBrush))
                rect_item.setZValue(40)
                self.scene.addItem(rect_item)
                self._select_preview_rect = rect_item
            self._select_preview_rect.setRect(QRectF(start_pos, start_pos))
            self._select_preview_rect.setVisible(True)

    def _update_selection_drag(self, scene_pos: QPointF) -> None:
        if self._select_drag_mode == "free" and self._select_path is not None:
            self._select_path.lineTo(scene_pos)
            if self._select_preview_path is not None:
                self._select_preview_path.setPath(self._select_path)
            return
        if self._select_drag_mode == "rect" and self._select_start_pos is not None:
            if self._select_preview_rect is not None:
                self._select_preview_rect.setRect(QRectF(self._select_start_pos, scene_pos).normalized())

    def _finalize_selection_drag(self) -> None:
        if self._select_drag_mode is None:
            return
        items: list = []
        if self._select_drag_mode == "free" and self._select_path is not None:
            path = QPainterPath(self._select_path)
            path.closeSubpath()
            candidates = [
                item
                for item in self.scene.items(path)
                if isinstance(
                    item,
                    (AtomItem, BondItem, ArrowItem, BracketItem, TextAnnotationItem),
                )
            ]
            items = [
                item
                for item in candidates
                if path.contains(item.sceneBoundingRect().center())
            ]
        elif self._select_drag_mode == "rect" and self._select_start_pos is not None:
            rect = QRectF(self._select_start_pos, self._last_scene_pos).normalized()
            items = [
                item
                for item in self.scene.items(rect)
                if isinstance(
                    item,
                    (AtomItem, BondItem, ArrowItem, BracketItem, TextAnnotationItem),
                )
            ]

        if not self._select_additive:
            self.scene.clearSelection()
        for item in items:
            item.setSelected(True)
        self._sync_selection_from_scene()
        self._clear_selection_drag()

    def _clear_selection_drag(self) -> None:
        self._select_drag_mode = None
        self._select_start_pos = None
        self._select_path = None
        self._select_additive = False
        if self._select_preview_path is not None:
            self.scene.removeItem(self._select_preview_path)
            self._select_preview_path = None
        if self._select_preview_rect is not None:
            self.scene.removeItem(self._select_preview_rect)
            self._select_preview_rect = None

    def _show_context_menu(
        self,
        scene_pos: QPointF,
        global_pos,
        clicked_item: Optional[QGraphicsItem] = None,
    ) -> None:
        menu = QMenu(self)
        has_selection = bool(
            self.state.selected_atoms
            or self.state.selected_bonds
            or self._selected_text_items()
            or any(isinstance(item, (ArrowItem, BracketItem)) for item in self.scene.selectedItems())
        )

        act_cut = menu.addAction("Cortar")
        act_copy = menu.addAction("Copiar")
        act_paste = menu.addAction("Pegar")
        menu.addSeparator()
        act_delete = menu.addAction("Eliminar")
        act_anchor = None
        if isinstance(clicked_item, AtomItem):
            atom = self.model.get_atom(clicked_item.atom_id)
            if atom is not None and atom.element not in ELEMENT_SYMBOLS:
                candidates = self._label_anchor_candidates(atom.element)
                if candidates:
                    menu.addSeparator()
                    act_anchor = menu.addAction("Elegir átomo de unión...")
        menu.addSeparator()
        act_select_all = menu.addAction("Seleccionar todo")

        act_cut.setEnabled(has_selection)
        act_copy.setEnabled(has_selection)
        act_delete.setEnabled(has_selection)

        action = menu.exec(global_pos)
        if action is None:
            return
        if act_anchor is not None and action == act_anchor and isinstance(clicked_item, AtomItem):
            self._prompt_anchor_for_atom(clicked_item.atom_id)
            return
        if action == act_copy:
            self.copy_to_clipboard()
            return
        if action == act_paste:
            self.paste_from_clipboard()
            return
        if action == act_cut:
            self.copy_to_clipboard()
            self.delete_selection()
            return
        if action == act_delete:
            self.delete_selection()
            return
        if action == act_select_all:
            self._select_all_items()


    def _ensure_selection_overlay(self) -> None:
        if self._selection_box is None:
            box = QGraphicsRectItem()
            pen = QPen(QColor("#4A90D9"), 0)
            box.setPen(pen)
            box.setBrush(QBrush(Qt.BrushStyle.NoBrush))
            box.setZValue(45)
            self.scene.addItem(box)
            self._selection_box = box
        if self._selection_handle is None:
            radius = SELECTION_HANDLE_RADIUS_PX
            handle = QGraphicsEllipseItem(0, 0, radius * 2.0, radius * 2.0)
            handle.setBrush(QBrush(QColor("#4A90D9")))
            handle.setPen(QPen(QColor("#4A90D9"), 0))
            handle.setZValue(46)
            handle.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIgnoresTransformations, True)
            self.scene.addItem(handle)
            self._selection_handle = handle

        if self._selection_move_handle is None:
            radius = SELECTION_HANDLE_RADIUS_PX
            handle = QGraphicsEllipseItem(0, 0, radius * 2.0, radius * 2.0)
            handle.setBrush(QBrush(QColor("#FFFFFF")))
            handle.setPen(QPen(QColor("#4A90D9"), 1))
            handle.setZValue(46)
            handle.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIgnoresTransformations, True)
            handle.setCursor(Qt.CursorShape.OpenHandCursor)
            self.scene.addItem(handle)
            self._selection_move_handle = handle

        if self._selection_scale_handle is None:
            size = SELECTION_HANDLE_RADIUS_PX * 2.0
            handle = QGraphicsRectItem(0, 0, size, size)
            handle.setBrush(QBrush(QColor("#4A90D9")))
            handle.setPen(QPen(QColor("#4A90D9"), 0))
            handle.setZValue(46)
            handle.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIgnoresTransformations, True)
            handle.setCursor(Qt.CursorShape.SizeFDiagCursor)
            self.scene.addItem(handle)
            self._selection_scale_handle = handle

    def _clear_selection_overlay(self) -> None:
        if self._selection_box is not None:
            self.scene.removeItem(self._selection_box)
            self._selection_box = None
        if self._selection_handle is not None:
            self.scene.removeItem(self._selection_handle)
            self._selection_handle = None
        if self._selection_move_handle is not None:
            self.scene.removeItem(self._selection_move_handle)
            self._selection_move_handle = None
        if self._selection_scale_handle is not None:
            self.scene.removeItem(self._selection_scale_handle)
            self._selection_scale_handle = None

    def _selected_items_bbox(self) -> Optional[QRectF]:
        rect: Optional[QRectF] = None

        def extend(candidate: QRectF) -> None:
            nonlocal rect
            if not candidate.isValid() or candidate.isNull():
                return
            rect = candidate if rect is None else rect.united(candidate)

        def extend_atom_bounds(atom_id: int) -> None:
            item = self.atom_items.get(atom_id)
            if item is None or item.scene() is not self.scene:
                return
            if item.pen().style() != Qt.PenStyle.NoPen or item.brush().style() != Qt.BrushStyle.NoBrush:
                extend(item.sceneBoundingRect())
            if item.label.isVisible():
                extend(item.label.sceneBoundingRect())
            if item.charge_label.isVisible():
                extend(item.charge_label.sceneBoundingRect())
            overlays = self._implicit_h_overlays.get(atom_id)
            if overlays:
                for line_item, text_item in overlays:
                    if line_item.scene() is self.scene and line_item.isVisible():
                        extend(line_item.sceneBoundingRect())
                    if text_item.scene() is self.scene and text_item.isVisible():
                        extend(text_item.sceneBoundingRect())

        for atom_id in self._selected_atom_ids_for_transform():
            extend_atom_bounds(atom_id)
        for bond_id in self.state.selected_bonds:
            item = self.bond_items.get(bond_id)
            if item is None:
                continue
            extend(item.sceneBoundingRect())
        for item in self.scene.selectedItems():
            # Include TextAnnotationItem
            if isinstance(item, (ArrowItem, BracketItem, TextAnnotationItem)):
                extend(item.sceneBoundingRect())
        return rect

    def _selected_atom_ids_for_transform(self) -> set[int]:
        atom_ids = set(self.state.selected_atoms)
        for bond_id in self.state.selected_bonds:
            if bond_id in self.model.bonds:
                bond = self.model.get_bond(bond_id)
                atom_ids.add(bond.a1_id)
                atom_ids.add(bond.a2_id)
        return atom_ids

    def _update_selection_overlay(self) -> None:
        bbox = self._selected_items_bbox()
        if bbox is None:
            if self._selection_box is not None:
                self._selection_box.setVisible(False)
            if self._selection_handle is not None:
                self._selection_handle.setVisible(False)
            if self._selection_move_handle is not None:
                self._selection_move_handle.setVisible(False)
            if self._selection_scale_handle is not None:
                self._selection_scale_handle.setVisible(False)
            return
        self._ensure_selection_overlay()
        padded = QRectF(bbox)
        pad = max(2.0, float(self.drawing_style.stroke_px))
        padded.adjust(-pad, -pad, pad, pad)
        if self._selection_box is not None:
            self._selection_box.setRect(padded)
            self._selection_box.setVisible(True)
        if self._selection_handle is not None:
            center_x = padded.center().x()
            handle_x = center_x - SELECTION_HANDLE_RADIUS_PX
            handle_y = padded.top() - SELECTION_HANDLE_OFFSET_PX - SELECTION_HANDLE_RADIUS_PX
            self._selection_handle.setPos(handle_x, handle_y)
            self._selection_handle.setVisible(True)

        if self._selection_move_handle is not None:
            center_x = padded.center().x()
            center_y = padded.center().y()
            handle_x = center_x - SELECTION_HANDLE_RADIUS_PX
            handle_y = center_y - SELECTION_HANDLE_RADIUS_PX
            self._selection_move_handle.setPos(handle_x, handle_y)
            self._selection_move_handle.setVisible(True)

        if self._selection_scale_handle is not None:
            handle_x = padded.right() - SELECTION_HANDLE_RADIUS_PX
            handle_y = padded.bottom() - SELECTION_HANDLE_RADIUS_PX
            self._selection_scale_handle.setPos(handle_x, handle_y)
            self._selection_scale_handle.setVisible(True)

    def _hit_selection_handle(self, scene_pos: QPointF) -> bool:
        if self._selection_handle is None or not self._selection_handle.isVisible():
            return False
        return self._selection_handle.sceneBoundingRect().contains(scene_pos)

    def _hit_selection_move_handle(self, scene_pos: QPointF) -> bool:
        if self._selection_move_handle is None or not self._selection_move_handle.isVisible():
            return False
        return self._selection_move_handle.sceneBoundingRect().contains(scene_pos)

    def _hit_selection_scale_handle(self, scene_pos: QPointF) -> bool:
        if self._selection_scale_handle is None or not self._selection_scale_handle.isVisible():
            return False
        return self._selection_scale_handle.sceneBoundingRect().contains(scene_pos)

    def _begin_rotation_drag(self, scene_pos: QPointF) -> None:
        if not self._selected_atom_ids_for_transform() and not self._selected_text_items():
            return
        bbox = self._selected_items_bbox()
        if bbox is None:
            return
        self._rotation_dragging = True
        self._rotation_center = bbox.center()
        dx = scene_pos.x() - self._rotation_center.x()
        dy = scene_pos.y() - self._rotation_center.y()
        self._rotation_start_angle = math.atan2(dy, dx)
        
        # Capture atoms
        atom_ids = self._selected_atom_ids_for_transform()
        self._rotation_start_positions = {
            atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
            for atom_id in atom_ids
            if atom_id in self.model.atoms
        }
        
        # Capture text items
        self._rotation_start_text_transforms = {} # item -> (pos, rotation)
        for item in self._selected_text_items():
            self._rotation_start_text_transforms[item] = (item.pos(), item.rotation())

    def _update_rotation_drag(self, scene_pos: QPointF) -> None:
        if self._rotation_center is None or self._rotation_start_angle is None:
            return
        dx = scene_pos.x() - self._rotation_center.x()
        dy = scene_pos.y() - self._rotation_center.y()
        current_angle = math.atan2(dy, dx)
        delta_angle = current_angle - self._rotation_start_angle
        degrees = math.degrees(delta_angle)

        self.rotate_selection_degrees(degrees, use_start_positions=True)

    def _finalize_rotation_drag(self) -> None:
        self._rotation_dragging = False
        self._rotation_center = None
        self._rotation_start_angle = None
        self._rotation_start_positions = {}
        self._rotation_start_text_transforms = {}
        self._update_selection_overlay()

    def _begin_scale_drag(self, scene_pos: QPointF) -> None:
        if not self._selected_atom_ids_for_transform() and not self._selected_text_items():
            return
        bbox = self._selected_items_bbox()
        if bbox is None or bbox.width() <= 1e-6 or bbox.height() <= 1e-6:
            return
        anchor = QPointF(bbox.left(), bbox.top())
        handle = QPointF(scene_pos.x(), scene_pos.y())
        start_len = math.hypot(handle.x() - anchor.x(), handle.y() - anchor.y())
        if start_len <= 1e-6:
            return

        self._scale_dragging = True
        self._scale_anchor = anchor
        self._scale_start_handle = handle
        self._scale_start_length = start_len

        atom_ids = self._selected_atom_ids_for_transform()
        self._scale_start_positions = {
            atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
            for atom_id in atom_ids
            if atom_id in self.model.atoms
        }

        self._scale_start_text_positions = {
            item: (item.pos(), item.rotation())
            for item in self._selected_text_items()
        }
        self._scale_has_moved = False

    def _update_scale_drag(self, scene_pos: QPointF) -> None:
        if not self._scale_dragging or self._scale_anchor is None:
            return
        if self._scale_start_length <= 1e-6:
            return
        dx = scene_pos.x() - self._scale_anchor.x()
        dy = scene_pos.y() - self._scale_anchor.y()
        current_len = math.hypot(dx, dy)
        scale = max(0.05, current_len / self._scale_start_length)
        if abs(scale - 1.0) > 1e-4:
            self._scale_has_moved = True

        if self._scale_start_positions:
            for atom_id, (x0, y0) in self._scale_start_positions.items():
                nx = self._scale_anchor.x() + (x0 - self._scale_anchor.x()) * scale
                ny = self._scale_anchor.y() + (y0 - self._scale_anchor.y()) * scale
                self.model.update_atom_position(atom_id, nx, ny)
                self.update_atom_item(atom_id, nx, ny)
            self.update_bond_items_for_atoms(set(self._scale_start_positions.keys()))

        if self._scale_start_text_positions:
            for item, (pos, rot) in self._scale_start_text_positions.items():
                nx = self._scale_anchor.x() + (pos.x() - self._scale_anchor.x()) * scale
                ny = self._scale_anchor.y() + (pos.y() - self._scale_anchor.y()) * scale
                item.setPos(QPointF(nx, ny))
                item.setRotation(rot)

        self._update_selection_overlay()

    def _finalize_scale_drag(self) -> None:
        if not self._scale_dragging:
            return

        if self._scale_has_moved:
            move_atoms = bool(self._scale_start_positions)
            move_text = bool(self._scale_start_text_positions)

            if move_atoms and move_text:
                self.undo_stack.beginMacro("Scale selection")

            if move_atoms:
                after = {
                    atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
                    for atom_id in self._scale_start_positions
                }
                cmd = MoveAtomsCommand(
                    self.model,
                    self,
                    self._scale_start_positions,
                    after,
                    skip_first_redo=True,
                )
                self.undo_stack.push(cmd)

            if move_text:
                before = self._scale_start_text_positions
                after = {item: (item.pos(), item.rotation()) for item in before.keys()}
                cmd = MoveTextItemsCommand(self, before, after)
                self.undo_stack.push(cmd)

            if move_atoms and move_text:
                self.undo_stack.endMacro()

        self._scale_dragging = False
        self._scale_anchor = None
        self._scale_start_handle = None
        self._scale_start_length = 0.0
        self._scale_start_positions = {}
        self._scale_start_text_positions = {}
        self._scale_has_moved = False
        self._update_selection_overlay()

    def _selected_text_items(self) -> list[TextAnnotationItem]:
        return [item for item in self.scene.selectedItems() if isinstance(item, TextAnnotationItem)]

    def rotate_selection_degrees(self, degrees: float, use_start_positions: bool = False) -> None:
        """Rotate selected items by degrees around their collective center."""
        selected_atom_ids = self._selected_atom_ids_for_transform()
        selected_text_items = self._selected_text_items()
        
        if not selected_atom_ids and not selected_text_items:
            return

        cx, cy = 0.0, 0.0
        
        # If dragging, usage rotation center
        if use_start_positions and self._rotation_center is not None:
             cx = self._rotation_center.x()
             cy = self._rotation_center.y()
             # We use the fixed center from drag start
        else:
             # Calculate collective center
            bbox = self._selected_items_bbox()
            if bbox is None:
                return
            center = bbox.center()
            cx, cy = center.x(), center.y()
            
        rad = math.radians(degrees)
        cos_a = math.cos(rad)
        sin_a = math.sin(rad)

        # Rotate Atoms
        if selected_atom_ids:
            new_positions = {}
            for aid in selected_atom_ids:
                if use_start_positions and aid in self._rotation_start_positions:
                    ox, oy = self._rotation_start_positions[aid]
                else:
                    atom = self.model.get_atom(aid)
                    ox, oy = atom.x, atom.y
                
                dx = ox - cx
                dy = oy - cy
                nx = dx * cos_a - dy * sin_a + cx
                ny = dx * sin_a + dy * cos_a + cy
                new_positions[aid] = (nx, ny)

            # If this is live drag, update directly without undo stack (handled at end)
            # Actually, canvas usually updates model directly during drag and pushes single command at end?
            # Existing _update_rotation_drag called rotate_selection_degrees.
            # But the original implementation of rotate_selection_degrees PUSHED COMMANDS.
            # We should probably separate "preview" from "commit".
            # For this refactor, let's assume this method updates positions efficiently.
            
            # Use MoveAtomsCommand for undo stack if valid drop?
            # Or just update positions if use_start_positions is True (dragging).
            if use_start_positions:
                for aid, (nx, ny) in new_positions.items():
                    self.model.update_atom_position(aid, nx, ny)
                    self.update_atom_item(aid, nx, ny)
                self.update_bond_items_for_atoms(selected_atom_ids)
            else:
                before = {
                    atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
                    for atom_id in selected_atom_ids
                }
                after = new_positions
                self.undo_stack.push(
                    MoveAtomsCommand(self.model, self, before, after)
                )

        # Rotate Text Items
        if selected_text_items:
             for item in selected_text_items:
                if use_start_positions and hasattr(self, "_rotation_start_text_transforms") and item in self._rotation_start_text_transforms:
                    start_pos, start_rot = self._rotation_start_text_transforms[item]
                    
                    # Logic: We are rotating around (cx, cy) by 'degrees' relative to START pose.
                    # QGraphicsItem rotation is around its transform origin.
                    
                    # 1. Calculate new center position (orbit)
                    # Item's scene pos is top-left usually? TextItem might be different.
                    # Let's use mapToScene(center) if we stored it?
                    # Getting intricate. Simplified:
                    # Treat item.pos() as the anchor.
                    
                    # Better: define item center relative to scene.
                    # Start Logic:
                    # We stored start_pos (item.pos) and start_rot.
                    
                    # We need the item's center at start state.
                    # But calculating that is hard if we don't store it.
                    # Let's assume start_pos is good enough for now, or improve storage.
                    
                    # Let's do delta rotation on CURRENT state? 
                    # use_start_positions implies absolute rotation from start.
                    
                    # Re-implement correctly for visual drag:
                    
                    # Calculate the item's center in scene coordinates based on its START position and rotation
                    # This is complex because QGraphicsTextItem's bounding rect changes with rotation.
                    # For simplicity during drag, let's apply the rotation to the item's current state
                    # relative to the rotation center, and then reset its rotation.
                    
                    # This is a tricky part. The current implementation of _update_rotation_drag
                    # calls this function with the *total* delta_angle from the start.
                    # So, we need to apply the rotation from the *original* state.
                    
                    # Reset to start state
                    item.setPos(start_pos)
                    item.setRotation(start_rot)
                    
                    # Now apply the full rotation from the start state
                    # Set origin to center for rotation
                    br = item.boundingRect()
                    center = br.center()
                    item.setTransformOriginPoint(center)
                    
                    # 1. Rotate in place
                    item.setRotation(start_rot + degrees) # Apply total rotation
                    
                    # 2. Orbit if we are rotating around a collective center that isn't just this item's center
                    item_scene_center = item.mapToScene(center)
                    # Only orbit if the item's center is not the rotation center (i.e., it's not a single item rotation around its own center)
                    if math.hypot(item_scene_center.x() - cx, item_scene_center.y() - cy) > 1.0:
                        # Calculate new center position after orbiting
                        dx_orbit = item_scene_center.x() - cx
                        dy_orbit = item_scene_center.y() - cy
                        nx_center_orbit = dx_orbit * cos_a - dy_orbit * sin_a + cx
                        ny_center_orbit = dx_orbit * sin_a + dy_orbit * cos_a + cy
                        
                        # Move item so its center is at nx_center_orbit, ny_center_orbit
                        offset = QPointF(nx_center_orbit, ny_center_orbit) - item_scene_center
                        item.moveBy(offset.x(), offset.y())
                else:
                    # Increment rotation (toolbar)
                    current_rot = item.rotation()
                    # Rotate around center of bbox? Or item center?
                    # Users usually expect in-place rotation if handled via toolbar unless group.
                    # For consistency with atoms, rotate around group center.
                    
                    item_center = item.sceneBoundingRect().center()
                    dx = item_center.x() - cx
                    dy = item_center.y() - cy
                    
                    nx_center = dx * cos_a - dy * sin_a + cx
                    ny_center = dx * sin_a + dy * cos_a + cy
                    
                    # Move item center to new position
                    # Diff
                    offset = QPointF(nx_center, ny_center) - item_center
                    item.moveBy(offset.x(), offset.y())
                    
                    # Add rotation
                    item.setTransformOriginPoint(item.boundingRect().center())
                    item.setRotation(current_rot + degrees)

        # Update overlay at the end
        if use_start_positions:
             self._update_selection_overlay()

    def flip_selection_horizontal(self) -> None:
        if not self._selected_atom_ids_for_transform() and not self._selected_text_items():
            return
        bbox = self._selected_items_bbox()
        if bbox is None:
            return
        cx = bbox.center().x()
        
        # Atoms
        atom_ids = self._selected_atom_ids_for_transform()
        # ... existing atom flip logic ...
        if atom_ids:
            before = {
                atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
                for atom_id in atom_ids
                if atom_id in self.model.atoms
            }
            after = {
                atom_id: (cx - (x0 - cx), y0)
                for atom_id, (x0, y0) in before.items()
            }
            cmd = MoveAtomsCommand(self.model, self, before, after)
            self.undo_stack.push(cmd)
            
        # Text Items
        # Flip position relative to cx
        for item in self._selected_text_items():
            # Flip center X
            br = item.sceneBoundingRect()
            center = br.center()
            new_cx = cx - (center.x() - cx)
            # Move
            item.moveBy(new_cx - center.x(), 0)
            # Potentially flip content? Usually text isn't mirrored.
            
        self._update_selection_overlay()

    def flip_selection_vertical(self) -> None:
        if not self._selected_atom_ids_for_transform() and not self._selected_text_items():
            return
        bbox = self._selected_items_bbox()
        if bbox is None:
            return
        cy = bbox.center().y()
        atom_ids = self._selected_atom_ids_for_transform()
        before = {
            atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
            for atom_id in atom_ids
            if atom_id in self.model.atoms
        }
        after = {
            atom_id: (x0, cy - (y0 - cy))
            for atom_id, (x0, y0) in before.items()
        }
        cmd = MoveAtomsCommand(self.model, self, before, after)
        self.undo_stack.push(cmd)
        self._update_selection_overlay()

    def _set_bond_anchor(self, atom_id: int, reset_angle: bool = False) -> None:
        if self.bond_anchor_id is not None and self.bond_anchor_id in self.atom_items:
            if self.bond_anchor_id not in self.state.selected_atoms:
                self.atom_items[self.bond_anchor_id].set_selected(False)
        self.bond_anchor_id = atom_id
        if atom_id in self.atom_items:
            self.atom_items[atom_id].set_selected(True)
        if reset_angle:
            self._bond_last_angle = None
            self._bond_zigzag_sign = 1

    def _clear_bond_anchor(self) -> None:
        if self.bond_anchor_id is not None and self.bond_anchor_id in self.atom_items:
            if self.bond_anchor_id not in self.state.selected_atoms:
                self.atom_items[self.bond_anchor_id].set_selected(False)
        self.bond_anchor_id = None
        self._bond_last_angle = None
        self._bond_zigzag_sign = 1

    def _create_or_update_bond(
        self,
        a1_id: int,
        a2_id: int,
        order: int,
        style: BondStyle,
        stereo: BondStereo,
        is_aromatic: bool,
    ) -> None:
        existing = self.model.find_bond_between(a1_id, a2_id)
        if existing is not None:
            if (
                self.state.active_bond_mode == "increment"
                and not is_aromatic
                and style == BondStyle.PLAIN
                and stereo == BondStereo.NONE
            ):
                new_order = 1 if existing.order >= 3 else existing.order + 1
                new_style = existing.style
                new_stereo = existing.stereo
                new_is_aromatic = False
            else:
                new_order = order
                new_style = style
                new_stereo = stereo
                new_is_aromatic = is_aromatic
            cmd = ChangeBondCommand(
                self.model,
                self,
                existing.id,
                new_order=new_order,
                new_style=new_style,
                new_stereo=new_stereo,
                new_is_aromatic=new_is_aromatic,
            )
            self.undo_stack.push(cmd)
            if new_is_aromatic:
                self._kekulize_aromatic_bonds(seed_atoms={a1_id, a2_id})
            return

        cmd = AddBondCommand(
            self.model,
            self,
            a1_id,
            a2_id,
            order,
            style,
            stereo,
            is_aromatic=is_aromatic,
        )
        self.undo_stack.push(cmd)
        if is_aromatic:
            self._kekulize_aromatic_bonds(seed_atoms={a1_id, a2_id})

    def _add_bond_between(self, a1_id: int, a2_id: int) -> None:
        order = self.state.active_bond_order
        style = self.state.active_bond_style
        stereo = self.state.active_bond_stereo
        is_aromatic = self.state.active_bond_aromatic
        if style != BondStyle.PLAIN or is_aromatic:
            order = 1
        self._create_or_update_bond(a1_id, a2_id, order, style, stereo, is_aromatic)

    def _cycle_bond_order(self, bond_id: int) -> None:
        bond = self.model.get_bond(bond_id)
        if bond.style != BondStyle.PLAIN or bond.stereo != BondStereo.NONE:
            cmd = ChangeBondCommand(
                self.model,
                self,
                bond_id,
                new_order=1,
                new_style=BondStyle.PLAIN,
                new_stereo=BondStereo.NONE,
                new_is_aromatic=False,
            )
            self.undo_stack.push(cmd)
            return
        new_order = 1 if bond.order >= 3 else bond.order + 1
        cmd = ChangeBondCommand(
            self.model,
            self,
            bond_id,
            new_order=new_order,
            new_is_aromatic=False,
        )
        self.undo_stack.push(cmd)

    def _apply_bond_style(self, bond_id: int) -> None:
        order = self.state.active_bond_order
        style = self.state.active_bond_style
        stereo = self.state.active_bond_stereo
        is_aromatic = self.state.active_bond_aromatic
        if style != BondStyle.PLAIN:
            order = 1
        if is_aromatic:
            order = 1
            style = BondStyle.PLAIN
            stereo = BondStereo.NONE
        cmd = ChangeBondCommand(
            self.model,
            self,
            bond_id,
            new_order=order,
            new_style=style,
            new_stereo=stereo,
            new_is_aromatic=is_aromatic,
        )
        self.undo_stack.push(cmd)
        if is_aromatic:
            bond = self.model.get_bond(bond_id)
            self._kekulize_aromatic_bonds(seed_atoms={bond.a1_id, bond.a2_id})

    def _create_first_bond(self, scene_pos: QPointF, modifiers: Qt.KeyboardModifiers) -> None:
        self.undo_stack.beginMacro("Add bond")
        atom_cmd = AddAtomCommand(
            self.model,
            self,
            self.state.default_element,
            scene_pos.x(),
            scene_pos.y(),
            expected_bonds=1,
        )
        self.undo_stack.push(atom_cmd)
        anchor_id = atom_cmd.atom_id
        if anchor_id is not None:
            self._set_bond_anchor(anchor_id, reset_angle=True)
            angle = self._default_bond_angle(anchor_id)
            self._add_bond_with_new_atom(anchor_id, angle)
        self.undo_stack.endMacro()

    def _extend_bond_from_anchor(self, scene_pos: QPointF, modifiers: Qt.KeyboardModifiers) -> None:
        anchor_id = self.bond_anchor_id
        if anchor_id is None:
            return
        anchor = self.model.get_atom(anchor_id)
        dx = scene_pos.x() - anchor.x
        dy = scene_pos.y() - anchor.y
        dist = math.hypot(dx, dy)
        step = self._bond_environment_step(anchor_id)
        angle_from_click = math.atan2(dy, dx) if dist > 1e-6 else 0.0
        existing_angles = self._get_anchor_bond_angles(anchor_id)
        if not (modifiers & Qt.KeyboardModifier.AltModifier) and len(existing_angles) <= 1:
            angle = self._default_bond_angle(anchor_id)
        else:
            near_anchor = dist < self.state.bond_length * 0.6
            aligned_with_existing = False
            if existing_angles:
                nearest = min(self._angle_distance(angle_from_click, a) for a in existing_angles)
                aligned_with_existing = nearest < math.radians(30)

            if near_anchor or aligned_with_existing:
                angle = self._default_bond_angle(anchor_id)
            else:
                angle = angle_from_click
                if not (modifiers & Qt.KeyboardModifier.AltModifier):
                    angle = self._snap_angle_to_environment(angle, anchor_id, step)
                self._reset_zigzag()

        self.undo_stack.beginMacro("Add bond")
        self._add_bond_with_new_atom(anchor_id, angle)
        self.undo_stack.endMacro()

    def _add_bond_with_new_atom(self, anchor_id: int, angle: float) -> None:
        anchor = self.model.get_atom(anchor_id)
        new_x = anchor.x + self.state.bond_length * math.cos(angle)
        new_y = anchor.y + self.state.bond_length * math.sin(angle)
        atom_cmd = AddAtomCommand(
            self.model,
            self,
            self.state.default_element,
            new_x,
            new_y,
            expected_bonds=1,
        )
        self.undo_stack.push(atom_cmd)
        new_atom_id = atom_cmd.atom_id
        if new_atom_id is None:
            return
        self._add_bond_between(anchor_id, new_atom_id)
        self._record_bond_angle(angle)
        self._set_bond_anchor(new_atom_id, reset_angle=False)

    def _record_bond_angle_between(self, a1_id: int, a2_id: int) -> None:
        a1 = self.model.get_atom(a1_id)
        a2 = self.model.get_atom(a2_id)
        angle = math.atan2(a2.y - a1.y, a2.x - a1.x)
        self._record_bond_angle(angle)

    def _record_bond_angle(self, angle: float) -> None:
        self._bond_last_angle = self._normalize_angle(angle)

    def _default_bond_angle(self, anchor_id: int) -> float:
        step = self._bond_environment_step(anchor_id)
        if self._bond_last_angle is not None:
            angle = self._bond_last_angle + step * self._bond_zigzag_sign
            if step < math.pi:
                self._bond_zigzag_sign *= -1
            return self._normalize_angle(angle)

        existing = self._get_anchor_bond_angles(anchor_id)
        if not existing:
            return 0.0
        candidates = []
        for base in existing:
            candidates.append(base + step)
            candidates.append(base - step)
        return self._best_separated_angle(candidates, existing)

    def _reset_zigzag(self) -> None:
        self._bond_zigzag_sign = 1

    def _snap_angle_to_environment(self, angle: float, anchor_id: int, step: float) -> float:
        if step <= 0:
            return angle
        existing = self._get_anchor_bond_angles(anchor_id)
        if not existing:
            return self._snap_angle(angle, step)
        candidates = []
        for base in existing:
            snapped = base + round((angle - base) / step) * step
            candidates.append(snapped)
        best = min(candidates, key=lambda a: self._angle_distance(a, angle))
        return self._normalize_angle(best)

    def _snap_angle(self, angle: float, step: float) -> float:
        if step <= 0:
            return angle
        return self._normalize_angle(round(angle / step) * step)

    def _bond_environment_step(self, anchor_id: int) -> float:
        max_order = max(1, self.state.active_bond_order)
        for bond in self.model.bonds.values():
            if bond.a1_id == anchor_id or bond.a2_id == anchor_id:
                order = 2 if bond.is_aromatic else bond.order
                max_order = max(max_order, order)
        if max_order >= 3:
            return math.pi
        if max_order == 2:
            return 2 * math.pi / 3
        return math.radians(self._sp3_display_angle_deg())

    def _get_anchor_bond_angles(self, anchor_id: int) -> List[float]:
        angles = []
        anchor = self.model.get_atom(anchor_id)
        for bond in self.model.bonds.values():
            if bond.a1_id == anchor_id:
                other = self.model.get_atom(bond.a2_id)
            elif bond.a2_id == anchor_id:
                other = self.model.get_atom(bond.a1_id)
            else:
                continue
            angle = math.atan2(other.y - anchor.y, other.x - anchor.x)
            angles.append(self._normalize_angle(angle))
        return angles

    def _best_separated_angle(self, candidates: List[float], existing: List[float]) -> float:
        best = None
        best_sep = -1.0
        for candidate in candidates:
            sep = min(self._angle_distance(candidate, a) for a in existing)
            if sep > best_sep:
                best_sep = sep
                best = candidate
        return self._normalize_angle(best) if best is not None else 0.0

    def _normalize_angle(self, angle: float) -> float:
        return (angle + math.pi * 2) % (math.pi * 2)

    def _angle_distance(self, a: float, b: float) -> float:
        diff = (a - b + math.pi) % (2 * math.pi) - math.pi
        return abs(diff)

    def _angle_snap_step_deg(self) -> float:
        if not self.state.fixed_angles:
            return 0.0
        step = float(self.state.angle_step_deg)
        return step if step > 0 else 0.0

    def _snap_angles_to_grid(self, angles: Iterable[float]) -> list[float]:
        step = self._angle_snap_step_deg()
        if step <= 0:
            return list(angles)
        snapped = [snap_angle_deg(angle, step) for angle in angles]
        seen: set[float] = set()
        deduped: list[float] = []
        for angle in snapped:
            key = round(angle, 6)
            if key in seen:
                continue
            seen.add(key)
            deduped.append(angle)
        return deduped

    def _sp3_display_angle_deg(self) -> float:
        step = self._angle_snap_step_deg()
        if step <= 0:
            return SP3_BOND_ANGLE_DEG
        return snap_angle_deg(SP3_BOND_ANGLE_DEG, step)

    # -------------------------------------------------------------------------
    # Hover + Drag State
    # -------------------------------------------------------------------------
    def _pick_hover_target(self, scene_pos: QPointF) -> tuple[Optional[int], Optional[int]]:
        atoms = [(atom.id, atom.x, atom.y) for atom in self.model.atoms.values()]
        atom_id = closest_atom(scene_pos, atoms, HOVER_ATOM_RADIUS)
        if atom_id is not None:
            return atom_id, None

        bonds = []
        for bond in self.model.bonds.values():
            a1 = self.model.get_atom(bond.a1_id)
            a2 = self.model.get_atom(bond.a2_id)
            bonds.append((bond.id, QPointF(a1.x, a1.y), QPointF(a2.x, a2.y)))
        bond_id = closest_bond(scene_pos, bonds, HOVER_BOND_DISTANCE)
        return None, bond_id

    def _update_hover(self, scene_pos: QPointF) -> None:
        atom_id, bond_id = self._pick_hover_target(scene_pos)

        if self.hovered_atom_id is not None and self.hovered_atom_id in self.atom_items:
            if self.hovered_atom_id not in self.state.selected_atoms:
                self.atom_items[self.hovered_atom_id].set_hover(False)

        self.hovered_atom_id = atom_id
        self.hovered_bond_id = bond_id if atom_id is None else None

        if self.hovered_atom_id is not None and self.hovered_atom_id in self.atom_items:
            if self.hovered_atom_id not in self.state.selected_atoms:
                self.atom_items[self.hovered_atom_id].set_hover(True)
            atom = self.model.get_atom(self.hovered_atom_id)
            self._hover_atom_indicator.update_position(atom.x, atom.y)
            self._hover_bond_indicator.hide_indicator()
            return

        self._hover_atom_indicator.hide_indicator()

        if self.hovered_bond_id is not None:
            bond = self.model.get_bond(self.hovered_bond_id)
            a1 = self.model.get_atom(bond.a1_id)
            a2 = self.model.get_atom(bond.a2_id)
            self._hover_bond_indicator.update_for_bond(
                QPointF(a1.x, a1.y),
                QPointF(a2.x, a2.y),
            )
        else:
            self._hover_bond_indicator.hide_indicator()

    def _cancel_drag(self) -> None:
        self._drag_mode = "none"
        self._drag_anchor = None
        self._ring_last_vertices = None
        self._chain_last_points = None
        self._drag_free_orientation = False
        self._bond_drag_start = None
        self._rotation_dragging = False
        self._scale_dragging = False
        self._rotation_center = None
        self._rotation_start_angle = None
        self._rotation_start_positions = {}
        self._clear_selection_drag()
        if self._overlays_ready:
            self._optimize_zone.hide_zone()
            self._preview_bond_item.hide_preview()
            self._preview_ring_item.hide_preview()
            self._preview_chain_item.hide_preview()
            self._preview_chain_label.hide_label()
            self._preview_arrow_item.hide_preview()

    def _on_undo_stack_changed(self, _index: int) -> None:
        saved_atoms = set(self.state.selected_atoms)
        saved_bonds = set(self.state.selected_bonds)
        self._cancel_drag()
        self._sync_scene_with_model()
        for atom_id in saved_atoms:
            item = self.atom_items.get(atom_id)
            if item is not None:
                item.setSelected(True)
        for bond_id in saved_bonds:
            item = self.bond_items.get(bond_id)
            if item is not None:
                item.setSelected(True)
        self._sync_selection_from_scene()
        self._kekulize_aromatic_bonds()
        self.show_valence_errors(self.model.validate())
        self._update_hover(self._last_scene_pos)
        self._update_selection_overlay()

    def _sync_scene_with_model(self) -> None:
        for item in list(self.scene.items()):
            if isinstance(item, (AtomItem, BondItem)):
                self.scene.removeItem(item)

        for atom_id in list(self.atom_items.keys()):
            self.remove_atom_item(atom_id)
        for bond_id in list(self.bond_items.keys()):
            self.remove_bond_item(bond_id)

        self._ring_centers.clear()
        ring_atoms: dict[int, set[int]] = {}
        for bond in self.model.bonds.values():
            if bond.ring_id is None:
                continue
            ring_atoms.setdefault(bond.ring_id, set()).update({bond.a1_id, bond.a2_id})
        for ring_id, atom_ids in ring_atoms.items():
            xs = [self.model.get_atom(aid).x for aid in atom_ids]
            ys = [self.model.get_atom(aid).y for aid in atom_ids]
            center = (sum(xs) / len(xs), sum(ys) / len(ys))
            self.register_ring_center(ring_id, center)

        for atom in self.model.atoms.values():
            self.add_atom_item(atom)
        for bond in self.model.bonds.values():
            self.add_bond_item(bond)

        self.refresh_atom_visibility()
        self.refresh_aromatic_circles()

    def show_valence_errors(self, errors: Iterable[int]) -> None:
        error_ids = set(errors)
        for atom_id, item in self.atom_items.items():
            item.set_valence_error(atom_id in error_ids)

    def _get_anchor_bond_angles_deg(self, anchor_id: int) -> list[float]:
        angles = []
        anchor = self.model.get_atom(anchor_id)
        p0 = QPointF(anchor.x, anchor.y)
        for bond in self.model.bonds.values():
            if bond.a1_id == anchor_id:
                other = self.model.get_atom(bond.a2_id)
            elif bond.a2_id == anchor_id:
                other = self.model.get_atom(bond.a1_id)
            else:
                continue
            p1 = QPointF(other.x, other.y)
            angles.append(angle_deg(p0, p1))
        return angles

    def _bond_geometry(self, anchor_id: int, bond_order: int, is_aromatic: bool) -> str:
        has_triple = bond_order >= 3
        has_double = bond_order == 2 or is_aromatic
        for bond in self.model.bonds.values():
            if bond.a1_id != anchor_id and bond.a2_id != anchor_id:
                continue
            if bond.order >= 3:
                has_triple = True
            if bond.order == 2 or bond.is_aromatic:
                has_double = True
        if has_triple:
            return "sp"
        if has_double:
            return "sp2"
        return "sp3"

    def _incoming_angle_deg(self, anchor_id: int) -> Optional[float]:
        angles = self._get_anchor_bond_angles_deg(anchor_id)
        if len(angles) == 1:
            return angles[0]
        return None

    def _direction_collision_metrics(
        self, anchor_id: Optional[int], p0: QPointF, p1: QPointF
    ) -> tuple[int, float, float]:
        intersections = 0
        min_atom_dist = float("inf")
        min_bond_dist = float("inf")
        atom_threshold = self.state.bond_length * MIN_ATOM_DIST_SCALE
        bond_threshold = self.state.bond_length * MIN_BOND_DIST_SCALE
        excluded_atoms: set[int] = set()
        if anchor_id is not None:
            excluded_atoms.add(anchor_id)
        target_id = closest_atom(
            p1,
            [(atom.id, atom.x, atom.y) for atom in self.model.atoms.values()],
            ATOM_HIT_RADIUS,
        )
        if target_id is not None:
            excluded_atoms.add(target_id)

        for atom in self.model.atoms.values():
            if atom.id in excluded_atoms:
                continue
            dist = math.hypot(p1.x() - atom.x, p1.y() - atom.y)
            min_atom_dist = min(min_atom_dist, dist)
            if dist < atom_threshold:
                intersections += 1

        for bond in self.model.bonds.values():
            if bond.a1_id in excluded_atoms or bond.a2_id in excluded_atoms:
                continue
            a1 = self.model.get_atom(bond.a1_id)
            a2 = self.model.get_atom(bond.a2_id)
            pA = QPointF(a1.x, a1.y)
            pB = QPointF(a2.x, a2.y)
            if segments_intersect(p0, p1, pA, pB):
                intersections += 1
            bond_dist = segment_min_distance(p0, p1, pA, pB)
            min_bond_dist = min(min_bond_dist, bond_dist)
            if bond_dist < bond_threshold:
                intersections += 1

        return intersections, min_atom_dist, min_bond_dist

    def _direction_score(
        self,
        anchor_id: Optional[int],
        p0: QPointF,
        p1: QPointF,
        angle_deg_value: float,
        mouse_angle_deg: float,
        preferred_deg: list[float],
    ) -> float:
        intersections, min_atom_dist, min_bond_dist = self._direction_collision_metrics(
            anchor_id, p0, p1
        )
        atom_threshold = self.state.bond_length * MIN_ATOM_DIST_SCALE
        bond_threshold = self.state.bond_length * MIN_BOND_DIST_SCALE

        score = angle_distance_deg(angle_deg_value, mouse_angle_deg)
        if preferred_deg:
            if min(angle_distance_deg(angle_deg_value, pref) for pref in preferred_deg) <= 15.0:
                score -= 15.0
        score += intersections * 100.0
        if min_atom_dist < atom_threshold:
            score += (atom_threshold - min_atom_dist) * 5.0
        if min_bond_dist < bond_threshold:
            score += (bond_threshold - min_bond_dist) * 5.0
        return score

    def _direction_is_valid(
        self,
        anchor_id: Optional[int],
        p0: QPointF,
        p1: QPointF,
    ) -> bool:
        intersections, min_atom_dist, min_bond_dist = self._direction_collision_metrics(
            anchor_id, p0, p1
        )
        atom_threshold = self.state.bond_length * MIN_ATOM_DIST_SCALE
        bond_threshold = self.state.bond_length * MIN_BOND_DIST_SCALE
        if intersections > 0:
            return False
        if min_atom_dist < atom_threshold:
            return False
        if min_bond_dist < bond_threshold:
            return False
        return True

    def _pick_bond_direction_deg(
        self,
        anchor: QPointF,
        anchor_id: Optional[int],
        mouse_angle_deg: float,
        bond_order: int,
        is_aromatic: bool,
        length: float,
        apply_collisions: bool = True,
        allow_length_boost: bool = True,
    ) -> tuple[float, float]:
        existing_angles = self._get_anchor_bond_angles_deg(anchor_id) if anchor_id else []
        geometry = (
            self._bond_geometry(anchor_id, bond_order, is_aromatic)
            if anchor_id is not None
            else geometry_for_bond(bond_order, is_aromatic, [])
        )
        candidates = candidate_directions_deg(geometry, existing_angles, mouse_angle_deg)
        if self.state.fixed_angles:
            candidates = self._snap_angles_to_grid(candidates)
        candidates = filter_occupied_angles_deg(
            candidates, existing_angles, ANGLE_OCCUPIED_TOLERANCE_DEG
        )
        if not candidates:
            candidates = candidate_directions_deg(geometry, [], mouse_angle_deg)
            if self.state.fixed_angles:
                candidates = self._snap_angles_to_grid(candidates)

        preferred: list[float] = []
        if geometry == "sp3" and anchor_id is not None:
            incoming = self._incoming_angle_deg(anchor_id)
            if incoming is not None:
                sp3_angle = self._sp3_display_angle_deg()
                preferred = [
                    (incoming + sp3_angle) % 360.0,
                    (incoming - sp3_angle) % 360.0,
                ]
                if self.state.fixed_angles:
                    preferred = self._snap_angles_to_grid(preferred)

        if not apply_collisions:
            picked = pick_closest_direction_deg(candidates, mouse_angle_deg, preferred)
            return (picked if picked is not None else mouse_angle_deg), length

        length_candidates = (length, length * COLLISION_LENGTH_BOOST) if allow_length_boost else (length,)
        valid_choices: list[tuple[float, float, float]] = []
        for length_candidate in length_candidates:
            for angle_candidate in candidates:
                p1 = endpoint_from_angle_len(anchor, angle_candidate, length_candidate)
                if not self._direction_is_valid(anchor_id, anchor, p1):
                    continue
                score = self._direction_score(
                    anchor_id,
                    anchor,
                    p1,
                    angle_candidate,
                    mouse_angle_deg,
                    preferred,
                )
                valid_choices.append((score, angle_candidate, length_candidate))

        if valid_choices:
            _, best_angle, best_length = min(valid_choices, key=lambda item: item[0])
            return best_angle, best_length

        best_score = None
        best_angle = None
        best_length = length
        for length_candidate in length_candidates:
            for angle_candidate in candidates:
                p1 = endpoint_from_angle_len(anchor, angle_candidate, length_candidate)
                score = self._direction_score(
                    anchor_id,
                    anchor,
                    p1,
                    angle_candidate,
                    mouse_angle_deg,
                    preferred,
                )
                if best_score is None or score < best_score:
                    best_score = score
                    best_angle = angle_candidate
                    best_length = length_candidate

        return (best_angle if best_angle is not None else mouse_angle_deg), best_length

    def _select_preferred_angle(
        self, anchor_id: int, cursor_angle: float, bond_order: int, is_aromatic: bool
    ) -> float:
        anchor = self.model.get_atom(anchor_id)
        p0 = QPointF(anchor.x, anchor.y)
        angle, _ = self._pick_bond_direction_deg(
            p0,
            anchor_id,
            cursor_angle,
            bond_order,
            is_aromatic,
            self.state.bond_length,
            apply_collisions=False,
            allow_length_boost=False,
        )
        return angle

    def _find_overlapping_bond(self, p0: QPointF, p1: QPointF) -> Optional[int]:
        for bond in self.model.bonds.values():
            a1 = self.model.get_atom(bond.a1_id)
            a2 = self.model.get_atom(bond.a2_id)
            b0 = QPointF(a1.x, a1.y)
            b1 = QPointF(a2.x, a2.y)
            if segments_nearly_equal(p0, p1, b0, b1, BOND_OVERLAP_TOLERANCE_PX):
                return bond.id
        return None

    def _snap_ring_vertex(self, pos: QPointF, excluded: set[int]) -> Optional[int]:
        candidates = [
            (atom.id, atom.x, atom.y)
            for atom in self.model.atoms.values()
            if atom.id not in excluded
        ]
        return closest_atom(pos, candidates, ATOM_HIT_RADIUS)

    def _build_ring_vertex_defs(
        self, vertices: List[QPointF], anchor_type: str, anchor_id: Optional[int]
    ) -> List[Tuple[Optional[int], float, float]]:
        vertex_defs: List[Tuple[Optional[int], float, float]] = []
        used_ids: set[int] = set()

        if anchor_type == "bond" and anchor_id is not None:
            bond = self.model.get_bond(anchor_id)
            used_ids.update({bond.a1_id, bond.a2_id})
            vertex_defs.append((bond.a1_id, vertices[0].x(), vertices[0].y()))
            vertex_defs.append((bond.a2_id, vertices[1].x(), vertices[1].y()))
            remaining_vertices = vertices[2:]
        elif anchor_type == "atom" and anchor_id is not None:
            used_ids.add(anchor_id)
            vertex_defs.append((anchor_id, vertices[0].x(), vertices[0].y()))
            remaining_vertices = vertices[1:]
        else:
            remaining_vertices = vertices

        for v in remaining_vertices:
            snapped_id = self._snap_ring_vertex(v, used_ids)
            if snapped_id is not None:
                used_ids.add(snapped_id)
                vertex_defs.append((snapped_id, v.x(), v.y()))
            else:
                vertex_defs.append((None, v.x(), v.y()))

        return vertex_defs

    # -------------------------------------------------------------------------
    # Bond Placement
    # -------------------------------------------------------------------------
    def _bond_drag_distance(self, scene_pos: QPointF) -> float:
        if self._bond_drag_start is None:
            return 0.0
        return (scene_pos - self._bond_drag_start).manhattanLength()

    def _update_bond_drag_state(self, scene_pos: QPointF) -> None:
        self._drag_free_orientation = (
            self._bond_drag_start is not None
            and self._bond_drag_distance(scene_pos) >= BOND_DRAG_THRESHOLD_PX
        )

    def _should_use_default_bond_angle(
        self, modifiers: Qt.KeyboardModifiers, scene_pos: Optional[QPointF] = None
    ) -> bool:
        if self._drag_mode != "place_bond":
            return False
        if modifiers & Qt.KeyboardModifier.AltModifier:
            return False
        if scene_pos is None:
            scene_pos = self._last_scene_pos
        if self._bond_drag_distance(scene_pos) >= BOND_DRAG_THRESHOLD_PX:
            return False
        return True

    def _current_bond_geometry(self) -> tuple[int, bool]:
        order = self.state.active_bond_order
        is_aromatic = self.state.active_bond_aromatic
        if self.state.active_bond_style != BondStyle.PLAIN or is_aromatic:
            order = 1
        return order, is_aromatic

    def _peek_default_bond_angle(self, anchor_id: Optional[int]) -> float:
        if anchor_id is None:
            return 0.0
        existing = self._get_anchor_bond_angles(anchor_id)
        step = self._bond_environment_step(anchor_id)
        if self._bond_last_angle is not None and existing:
            tolerance = math.radians(BOND_LAST_ANGLE_TOLERANCE_DEG)
            if min(self._angle_distance(self._bond_last_angle, a) for a in existing) <= tolerance:
                angle = self._bond_last_angle + step * self._bond_zigzag_sign
                return self._normalize_angle(angle)
        if not existing:
            return 0.0
        candidates = []
        for base in existing:
            candidates.append(base + step)
            candidates.append(base - step)
        return self._best_separated_angle(candidates, existing)

    def _compute_default_bond_endpoint(
        self, anchor: QPointF, anchor_id: Optional[int]
    ) -> QPointF:
        order, aromatic = self._current_bond_geometry()
        default_angle_deg = math.degrees(self._peek_default_bond_angle(anchor_id))
        length = self.state.bond_length
        theta, final_length = self._pick_bond_direction_deg(
            anchor,
            anchor_id,
            default_angle_deg,
            order,
            aromatic,
            length,
            apply_collisions=True,
            allow_length_boost=self.state.fixed_lengths,
        )
        p1 = endpoint_from_angle_len(anchor, theta, final_length)
        snap_id = closest_atom(
            p1,
            [(atom.id, atom.x, atom.y) for atom in self.model.atoms.values()],
            ATOM_HIT_RADIUS,
        )
        if snap_id is not None and (anchor_id is None or snap_id != anchor_id):
            atom = self.model.get_atom(snap_id)
            return QPointF(atom.x, atom.y)
        return p1

    def _update_bond_angle_state(
        self,
        p0: QPointF,
        p1: QPointF,
        used_default: bool,
        anchor_id: Optional[int],
    ) -> None:
        angle = math.atan2(p1.y() - p0.y(), p1.x() - p0.x())
        self._bond_last_angle = self._normalize_angle(angle)
        if not used_default:
            self._bond_zigzag_sign = 1
            return
        step = (
            self._bond_environment_step(anchor_id)
            if anchor_id is not None
            else math.radians(self._sp3_display_angle_deg())
        )
        if step < math.pi:
            self._bond_zigzag_sign *= -1

    def _begin_place_bond(self, anchor_id: Optional[int], scene_pos: QPointF) -> None:
        self._drag_mode = "place_bond"
        self._drag_free_orientation = False
        self._bond_drag_start = scene_pos
        if anchor_id is None:
            self._drag_anchor = {"type": "free", "id": None, "pos": scene_pos}
            self._optimize_zone.hide_zone()
        else:
            anchor = self.model.get_atom(anchor_id)
            self._drag_anchor = {"type": "atom", "id": anchor_id, "pos": QPointF(anchor.x, anchor.y)}
            radius = self.state.bond_length * OPTIMIZE_ZONE_SCALE
            self._optimize_zone.set_radius(radius)
            self._optimize_zone.update_center(anchor.x, anchor.y)
        self._update_bond_preview(scene_pos, Qt.KeyboardModifier.NoModifier)

    def _update_bond_preview(self, scene_pos: QPointF, modifiers: Qt.KeyboardModifiers) -> None:
        if self._drag_anchor is None:
            return
        if self._drag_anchor["id"] is None:
            p0 = QPointF(self._drag_anchor["pos"].x(), self._drag_anchor["pos"].y())
        else:
            anchor = self.model.get_atom(self._drag_anchor["id"])
            p0 = QPointF(anchor.x, anchor.y)
        if self._should_use_default_bond_angle(modifiers, scene_pos):
            p1 = self._compute_default_bond_endpoint(p0, self._drag_anchor["id"])
        else:
            p1 = self._compute_bond_endpoint(p0, scene_pos, modifiers)
        self._preview_bond_item.update_line(p0, p1)

    def _finalize_bond(self, modifiers: Qt.KeyboardModifiers) -> None:
        if self._drag_anchor is None:
            self._cancel_drag()
            return
        anchor_id = self._drag_anchor["id"]
        if anchor_id is None:
            p0 = QPointF(self._drag_anchor["pos"].x(), self._drag_anchor["pos"].y())
        else:
            anchor = self.model.get_atom(anchor_id)
            p0 = QPointF(anchor.x, anchor.y)
        use_default = self._should_use_default_bond_angle(modifiers, self._last_scene_pos)
        if use_default:
            p1 = self._compute_default_bond_endpoint(p0, anchor_id)
        else:
            p1 = self._compute_bond_endpoint(p0, self._last_scene_pos, modifiers)

        atom_positions = [(atom.id, atom.x, atom.y) for atom in self.model.atoms.values()]
        target_id = closest_atom(self._last_scene_pos, atom_positions, ATOM_HIT_RADIUS)
        if target_id is None:
            target_id = closest_atom(p1, atom_positions, ATOM_HIT_RADIUS)
        if target_id == anchor_id:
            target_id = None

        if target_id is None:
            overlapping_bond = self._find_overlapping_bond(p0, p1)
            if overlapping_bond is not None:
                if self.state.active_bond_mode == "increment" and not self.state.active_bond_aromatic:
                    self._cycle_bond_order(overlapping_bond)
                else:
                    self._apply_bond_style(overlapping_bond)
                self._cancel_drag()
                return
            bond_hit = self._get_bond_at(self._last_scene_pos)
            if bond_hit is not None:
                if self.state.active_bond_mode == "increment" and not self.state.active_bond_aromatic:
                    self._cycle_bond_order(bond_hit)
                else:
                    self._apply_bond_style(bond_hit)
                self._cancel_drag()
                return

        order = self.state.active_bond_order
        style = self.state.active_bond_style
        stereo = self.state.active_bond_stereo
        is_aromatic = self.state.active_bond_aromatic
        if style != BondStyle.PLAIN or is_aromatic:
            order = 1

        final_p1 = p1
        if target_id is not None:
            target_atom = self.model.get_atom(target_id)
            final_p1 = QPointF(target_atom.x, target_atom.y)

        if anchor_id is None:
            self.undo_stack.beginMacro("Add bond")
            anchor_cmd = AddAtomCommand(
                self.model,
                self,
                self.state.default_element,
                p0.x(),
                p0.y(),
                expected_bonds=1,
            )
            self.undo_stack.push(anchor_cmd)
            anchor_id = anchor_cmd.atom_id
            if target_id is not None:
                self._create_or_update_bond(
                    anchor_id,
                    target_id,
                    order,
                    style,
                    stereo,
                    is_aromatic,
                )
            else:
                cmd = AddBondCommand(
                    self.model,
                    self,
                    anchor_id,
                    None,
                    order,
                    style,
                    stereo,
                    is_aromatic=is_aromatic,
                    new_atom_element=self.state.default_element,
                    new_atom_pos=(p1.x(), p1.y()),
                )
                self.undo_stack.push(cmd)
                if is_aromatic:
                    self._kekulize_aromatic_bonds()
            self.undo_stack.endMacro()
            self._update_bond_angle_state(p0, final_p1, use_default, anchor_id)
        else:
            if target_id is not None:
                self._create_or_update_bond(
                    anchor_id,
                    target_id,
                    order,
                    style,
                    stereo,
                    is_aromatic,
                )
            else:
                cmd = AddBondCommand(
                    self.model,
                    self,
                    anchor_id,
                    None,
                    order,
                    style,
                    stereo,
                    is_aromatic=is_aromatic,
                    new_atom_element=self.state.default_element,
                    new_atom_pos=(p1.x(), p1.y()),
                )
                self.undo_stack.push(cmd)
                if is_aromatic:
                    self._kekulize_aromatic_bonds()
            self._update_bond_angle_state(p0, final_p1, use_default, anchor_id)
        self._cancel_drag()

    def _compute_bond_endpoint(
        self,
        anchor: QPointF,
        scene_pos: QPointF,
        modifiers: Qt.KeyboardModifiers,
        bond_order: Optional[int] = None,
        is_aromatic: Optional[bool] = None,
    ) -> QPointF:
        dx = scene_pos.x() - anchor.x()
        dy = scene_pos.y() - anchor.y()
        dist = math.hypot(dx, dy)

        use_alt = bool(modifiers & Qt.KeyboardModifier.AltModifier)
        use_shift = bool(modifiers & Qt.KeyboardModifier.ShiftModifier)
        use_optimize = False
        if self._optimize_zone.isVisible():
            use_optimize = dist <= self._optimize_zone.radius()

        cursor_theta = angle_deg(anchor, scene_pos)
        free_drag = self._drag_mode == "place_bond" and self._drag_free_orientation
        if not self.state.fixed_angles or use_alt or free_drag:
            theta = cursor_theta
            length = self.state.bond_length if self.state.fixed_lengths and not use_shift else max(5.0, dist)
            p1 = endpoint_from_angle_len(anchor, theta, length)
            snap_id = closest_atom(
                p1,
                [(atom.id, atom.x, atom.y) for atom in self.model.atoms.values()],
                ATOM_HIT_RADIUS,
            )
            if snap_id is not None and (self._drag_anchor is None or snap_id != self._drag_anchor["id"]):
                atom = self.model.get_atom(snap_id)
                return QPointF(atom.x, atom.y)
            return p1

        order = bond_order if bond_order is not None else self.state.active_bond_order
        aromatic_flag = is_aromatic if is_aromatic is not None else self.state.active_bond_aromatic
        style = (
            self.state.active_bond_style if bond_order is None and is_aromatic is None else BondStyle.PLAIN
        )
        if style != BondStyle.PLAIN or aromatic_flag:
            order = 1

        if use_optimize and self._drag_anchor["id"] is not None:
            cursor_theta = self._select_preferred_angle(
                self._drag_anchor["id"], cursor_theta, order, aromatic_flag
            )

        length = self.state.bond_length if self.state.fixed_lengths and not use_shift else max(5.0, dist)
        theta, final_length = self._pick_bond_direction_deg(
            anchor,
            self._drag_anchor["id"] if self._drag_anchor else None,
            cursor_theta,
            order,
            aromatic_flag,
            length,
            apply_collisions=True,
            allow_length_boost=self.state.fixed_lengths and not use_shift,
        )
        p1 = endpoint_from_angle_len(anchor, theta, final_length)
        snap_id = closest_atom(
            p1,
            [(atom.id, atom.x, atom.y) for atom in self.model.atoms.values()],
            ATOM_HIT_RADIUS,
        )
        if snap_id is not None and (self._drag_anchor is None or snap_id != self._drag_anchor["id"]):
            atom = self.model.get_atom(snap_id)
            return QPointF(atom.x, atom.y)
        return p1

    # -------------------------------------------------------------------------
    # Ring Placement
    # -------------------------------------------------------------------------
    def _begin_place_ring(self, anchor_type: str, anchor_id: Optional[int], scene_pos: QPointF) -> None:
        self._drag_mode = "place_ring"
        self._drag_anchor = {"type": anchor_type, "id": anchor_id, "pos": scene_pos}
        if anchor_type == "atom":
            atom = self.model.get_atom(anchor_id)
            radius = self.state.bond_length * OPTIMIZE_ZONE_SCALE
            self._optimize_zone.set_radius(radius)
            self._optimize_zone.update_center(atom.x, atom.y)
        elif anchor_type == "bond":
            bond = self.model.get_bond(anchor_id)
            a1 = self.model.get_atom(bond.a1_id)
            a2 = self.model.get_atom(bond.a2_id)
            mid = QPointF((a1.x + a2.x) / 2, (a1.y + a2.y) / 2)
            radius = self.state.bond_length * OPTIMIZE_ZONE_SCALE
            self._optimize_zone.set_radius(radius)
            self._optimize_zone.update_center(mid.x(), mid.y())
        else:
            self._optimize_zone.hide_zone()
        self._update_ring_preview(scene_pos, Qt.KeyboardModifier.NoModifier)

    def _update_ring_preview(self, scene_pos: QPointF, modifiers: Qt.KeyboardModifiers) -> None:
        vertices = self._compute_ring_vertices(scene_pos, modifiers)
        self._ring_last_vertices = vertices
        self._preview_ring_item.update_polygon(vertices)

    def _finalize_ring(self, modifiers: Qt.KeyboardModifiers) -> None:
        vertices = self._ring_last_vertices or []
        if not vertices or self._drag_anchor is None:
            self._cancel_drag()
            return

        ring_size = max(3, int(self.state.active_ring_size))
        aromatic = self.state.active_ring_aromatic

        vertex_defs = self._build_ring_vertex_defs(
            vertices, self._drag_anchor["type"], self._drag_anchor["id"]
        )

        edge_defs: List[Tuple[int, int, int, BondStyle, BondStereo, bool]] = []
        if aromatic:
            edge_defs = self._build_aromatic_edges(vertex_defs, ring_size)
        else:
            for i in range(ring_size):
                j = (i + 1) % ring_size
                edge_defs.append((i, j, 1, BondStyle.PLAIN, BondStereo.NONE, aromatic))

        cmd = AddRingCommand(
            self.model,
            self,
            vertex_defs,
            edge_defs,
            element=self.state.default_element,
        )
        self.undo_stack.push(cmd)
        if aromatic:
            self._kekulize_aromatic_bonds()
        self._cancel_drag()

    def _compute_ring_vertices(
        self, scene_pos: QPointF, modifiers: Qt.KeyboardModifiers
    ) -> List[QPointF]:
        if self._drag_anchor is None:
            return []

        ring_size = max(3, int(self.state.active_ring_size))
        anchor_type = self._drag_anchor["type"]

        if anchor_type == "free":
            center = self._drag_anchor["pos"]
            angle0 = angle_deg(center, scene_pos)
            radius = self.state.bond_length / (2 * math.sin(math.pi / ring_size))
            vertices = []
            step = 360.0 / ring_size
            for i in range(ring_size):
                theta = angle0 + step * i
                vertices.append(endpoint_from_angle_len(center, theta, radius))
            return vertices

        if anchor_type == "bond":
            bond = self.model.get_bond(self._drag_anchor["id"])
            a1 = self.model.get_atom(bond.a1_id)
            a2 = self.model.get_atom(bond.a2_id)
            p1 = QPointF(a1.x, a1.y)
            p2 = QPointF(a2.x, a2.y)
            direction = bond_side(p1, p2, scene_pos)
            return self._regular_polygon_from_edge(p1, p2, ring_size, direction)

        anchor = self.model.get_atom(self._drag_anchor["id"])
        p1 = QPointF(anchor.x, anchor.y)
        p2 = self._compute_bond_endpoint(
            p1,
            scene_pos,
            modifiers,
            bond_order=1,
            is_aromatic=self.state.active_ring_aromatic,
        )
        direction = bond_side(p1, p2, scene_pos)
        return self._regular_polygon_from_edge(p1, p2, ring_size, direction)

    def _regular_polygon_from_edge(
        self,
        p1: QPointF,
        p2: QPointF,
        ring_size: int,
        direction: int,
    ) -> List[QPointF]:
        dx = p2.x() - p1.x()
        dy = p2.y() - p1.y()
        length = math.hypot(dx, dy)
        if length < 1e-6:
            return []

        step = 360.0 / ring_size
        angle0 = angle_deg(p1, p2)
        vertices = [p1, p2]
        current = p2
        current_angle = angle0
        for _ in range(ring_size - 2):
            current_angle = (current_angle + direction * step) % 360.0
            next_point = endpoint_from_angle_len(current, current_angle, length)
            vertices.append(next_point)
            current = next_point
        return vertices

    # -------------------------------------------------------------------------
    # Chain Placement
    # -------------------------------------------------------------------------
    def _begin_place_chain(self, anchor_id: Optional[int], scene_pos: QPointF) -> None:
        self._drag_mode = "place_chain"
        if anchor_id is None:
            self._drag_anchor = {"type": "free", "id": None, "pos": QPointF(scene_pos)}
            self._optimize_zone.hide_zone()
        else:
            anchor = self.model.get_atom(anchor_id)
            self._drag_anchor = {"type": "atom", "id": anchor_id, "pos": QPointF(anchor.x, anchor.y)}
            radius = self.state.bond_length * OPTIMIZE_ZONE_SCALE
            self._optimize_zone.set_radius(radius)
            self._optimize_zone.update_center(anchor.x, anchor.y)
        self._update_chain_preview(scene_pos, Qt.KeyboardModifier.NoModifier)

    def _update_chain_preview(self, scene_pos: QPointF, modifiers: Qt.KeyboardModifiers) -> None:
        if self._drag_anchor is None:
            return
        anchor_id = self._drag_anchor["id"]
        if anchor_id is None:
            p0 = QPointF(self._drag_anchor["pos"])
        else:
            anchor = self.model.get_atom(anchor_id)
            p0 = QPointF(anchor.x, anchor.y)
        dx = scene_pos.x() - p0.x()
        dy = scene_pos.y() - p0.y()
        dist = math.hypot(dx, dy)

        if dist < 1.0:
            points = [p0]
            self._chain_last_points = points
            self._preview_chain_item.update_polyline(points)
            self._preview_chain_label.hide_label()
            return

        raw_angle = angle_deg(p0, scene_pos)
        use_alt = bool(modifiers & Qt.KeyboardModifier.AltModifier)
        use_shift = bool(modifiers & Qt.KeyboardModifier.ShiftModifier)
        use_optimize = False
        if self._optimize_zone.isVisible():
            use_optimize = dist <= self._optimize_zone.radius()

        if use_shift:
            if abs(dx) >= abs(dy):
                raw_angle = 0.0 if dx >= 0 else 180.0
            else:
                raw_angle = 90.0 if -(dy) >= 0 else 270.0

        if not self.state.fixed_angles or use_alt:
            base_angle = raw_angle
        else:
            cursor_angle = raw_angle
            if use_optimize and anchor_id is not None:
                cursor_angle = self._select_preferred_angle(anchor_id, cursor_angle, 1, False)
            base_angle, _ = self._pick_bond_direction_deg(
                p0,
                anchor_id,
                cursor_angle,
                1,
                False,
                self.state.bond_length,
                apply_collisions=True,
                allow_length_boost=self.state.fixed_lengths and not use_shift,
            )

        n = max(1, min(CHAIN_MAX_BONDS, int(round(dist / self.state.bond_length))))
        points = [p0]
        geometry = self._bond_geometry(anchor_id, 1, False)
        zigzag = geometry == "sp3" and not use_alt and self.state.fixed_angles
        sp3_angle = self._sp3_display_angle_deg()
        zigzag_delta = (180.0 - sp3_angle) / 2.0
        if zigzag:
            incoming = (
                self._incoming_angle_deg(anchor_id)
                if anchor_id is not None
                else None
            )
            if incoming is not None:
                target_angles = [
                    (incoming + sp3_angle) % 360.0,
                    (incoming - sp3_angle) % 360.0,
                ]
                axis_options = [(target - zigzag_delta) % 360.0 for target in target_angles]
                axis_options = self._snap_angles_to_grid(axis_options)
                base_angle = min(axis_options, key=lambda ang: angle_distance_deg(ang, base_angle))
        current = p0
        for i in range(1, n + 1):
            if zigzag:
                step_angle = base_angle + (zigzag_delta if i % 2 == 1 else -zigzag_delta)
            else:
                step_angle = base_angle
            next_point = endpoint_from_angle_len(current, step_angle, self.state.bond_length)
            points.append(next_point)
            current = next_point

        self._chain_last_points = points
        self._preview_chain_item.update_polyline(points)
        self._preview_chain_label.update_label(str(n), points[-1] + QPointF(0, -10))

    def _finalize_chain(self, modifiers: Qt.KeyboardModifiers) -> None:
        if not self._chain_last_points or self._drag_anchor is None:
            self._cancel_drag()
            return
        anchor_id = self._drag_anchor["id"]
        new_positions = [(p.x(), p.y()) for p in self._chain_last_points[1:]]
        if not new_positions:
            self._cancel_drag()
            return
        anchor_pos = None
        if anchor_id is None:
            anchor_point = QPointF(self._drag_anchor["pos"])
            anchor_pos = (anchor_point.x(), anchor_point.y())
        cmd = AddChainCommand(
            self.model,
            self,
            anchor_id,
            new_positions,
            element=self.state.default_element,
            anchor_position=anchor_pos,
        )
        self.undo_stack.push(cmd)
        self._cancel_drag()

    # -------------------------------------------------------------------------
    # Hotkeys
    # -------------------------------------------------------------------------
    def _handle_hotkeys(self, event) -> bool:
        key = event.key()
        if Qt.Key.Key_0 <= key <= Qt.Key.Key_9:
            digit = key - Qt.Key.Key_0
            if digit == 0:
                return False
            if self.hovered_atom_id is not None:
                if event.modifiers() & Qt.KeyboardModifier.ShiftModifier:
                    if digit >= 3:
                        self._create_ring_hotkey(self.hovered_atom_id, digit)
                        return True
                else:
                    self._create_chain_hotkey(self.hovered_atom_id, digit)
                    return True
            if self.hovered_bond_id is not None:
                if event.modifiers() & Qt.KeyboardModifier.ShiftModifier:
                    if digit >= 3:
                        self._create_ring_from_bond_hotkey(self.hovered_bond_id, digit)
                        return True
                else:
                    self._set_bond_order_hotkey(self.hovered_bond_id, digit)
                    return True
        return False

    def _create_chain_hotkey(self, anchor_id: int, n: int) -> None:
        anchor = self.model.get_atom(anchor_id)
        p0 = QPointF(anchor.x, anchor.y)
        angles = self._get_anchor_bond_angles_deg(anchor_id)
        theta = choose_optimal_direction(angles)
        positions = []
        for i in range(1, n + 1):
            p = endpoint_from_angle_len(p0, theta, self.state.bond_length * i)
            positions.append((p.x(), p.y()))
        cmd = AddChainCommand(self.model, self, anchor_id, positions, self.state.default_element)
        self.undo_stack.push(cmd)

    def _create_ring_hotkey(self, anchor_id: int, ring_size: int) -> None:
        anchor = self.model.get_atom(anchor_id)
        p0 = QPointF(anchor.x, anchor.y)
        angles = self._get_anchor_bond_angles_deg(anchor_id)
        theta = choose_optimal_direction(angles)
        p1 = endpoint_from_angle_len(p0, theta, self.state.bond_length)
        direction = bond_side(p0, p1, self._last_scene_pos)
        vertices = self._regular_polygon_from_edge(p0, p1, ring_size, direction)
        self._commit_ring_vertices(vertices, anchor_type="atom", anchor_id=anchor_id)

    def _create_ring_from_bond_hotkey(self, bond_id: int, ring_size: int) -> None:
        bond = self.model.get_bond(bond_id)
        a1 = self.model.get_atom(bond.a1_id)
        a2 = self.model.get_atom(bond.a2_id)
        p1 = QPointF(a1.x, a1.y)
        p2 = QPointF(a2.x, a2.y)
        direction = bond_side(p1, p2, self._last_scene_pos)
        vertices = self._regular_polygon_from_edge(p1, p2, ring_size, direction)
        self._commit_ring_vertices(vertices, anchor_type="bond", anchor_id=bond_id)

    def _commit_ring_vertices(
        self, vertices: List[QPointF], anchor_type: str, anchor_id: int
    ) -> None:
        if not vertices:
            return
        ring_size = len(vertices)
        aromatic = self.state.active_ring_aromatic
        vertex_defs = self._build_ring_vertex_defs(vertices, anchor_type, anchor_id)

        edge_defs: List[Tuple[int, int, int, BondStyle, BondStereo, bool]] = []
        if aromatic:
            edge_defs = self._build_aromatic_edges(vertex_defs, ring_size)
        else:
            for i in range(ring_size):
                j = (i + 1) % ring_size
                edge_defs.append((i, j, 1, BondStyle.PLAIN, BondStereo.NONE, aromatic))

        cmd = AddRingCommand(
            self.model,
            self,
            vertex_defs,
            edge_defs,
            element=self.state.default_element,
        )
        self.undo_stack.push(cmd)
        if aromatic:
            self._kekulize_aromatic_bonds()

    def _set_bond_order_hotkey(self, bond_id: int, order: int) -> None:
        order = max(1, min(3, order))
        cmd = ChangeBondCommand(
            self.model,
            self,
            bond_id,
            new_order=order,
            new_style=BondStyle.PLAIN,
            new_stereo=BondStereo.NONE,
            new_is_aromatic=False,
        )
        self.undo_stack.push(cmd)

    def _kekulize_aromatic_bonds(self, seed_atoms: Optional[Iterable[int]] = None) -> None:
        aromatic_bonds = [bond for bond in self.model.bonds.values() if bond.is_aromatic]
        if not aromatic_bonds:
            return

        component_atoms: Optional[set[int]] = None
        if seed_atoms is not None:
            adjacency: dict[int, list[int]] = {}
            for bond in aromatic_bonds:
                adjacency.setdefault(bond.a1_id, []).append(bond.a2_id)
                adjacency.setdefault(bond.a2_id, []).append(bond.a1_id)
            seeds = [atom_id for atom_id in seed_atoms if atom_id in adjacency]
            if not seeds:
                return
            component_atoms = set()
            stack = list(seeds)
            while stack:
                node = stack.pop()
                if node in component_atoms:
                    continue
                component_atoms.add(node)
                for neighbor in adjacency.get(node, []):
                    if neighbor not in component_atoms:
                        stack.append(neighbor)

        display_orders = kekulize_display_orders(self.model, seed_atoms=component_atoms)

        if display_orders is None:
            try:
                import networkx as nx
            except Exception:
                nx = None

            if nx is not None:
                graph = nx.Graph()
                bond_by_pair: dict[frozenset[int], int] = {}
                for bond in aromatic_bonds:
                    if component_atoms is not None and (
                        bond.a1_id not in component_atoms or bond.a2_id not in component_atoms
                    ):
                        continue
                    pair = frozenset({bond.a1_id, bond.a2_id})
                    bond_by_pair[pair] = bond.id
                    weight = 10_000_000 - (min(pair) * 1000 + max(pair))
                    graph.add_edge(bond.a1_id, bond.a2_id, weight=weight)
                matching = nx.max_weight_matching(graph, maxcardinality=True, weight="weight")
                display_orders = {bond_id: 1 for bond_id in bond_by_pair.values()}
                for u, v in matching:
                    bond_id = bond_by_pair.get(frozenset({u, v}))
                    if bond_id is not None:
                        display_orders[bond_id] = 2

        if display_orders is None:
            display_orders = {}
            used_atoms: set[int] = set()
            sorted_bonds = sorted(
                aromatic_bonds,
                key=lambda bond: (min(bond.a1_id, bond.a2_id), max(bond.a1_id, bond.a2_id)),
            )
            for bond in sorted_bonds:
                if component_atoms is not None and (
                    bond.a1_id not in component_atoms or bond.a2_id not in component_atoms
                ):
                    continue
                display_orders[bond.id] = 1
            for bond in sorted_bonds:
                if component_atoms is not None and (
                    bond.a1_id not in component_atoms or bond.a2_id not in component_atoms
                ):
                    continue
                if bond.a1_id in used_atoms or bond.a2_id in used_atoms:
                    continue
                display_orders[bond.id] = 2
                used_atoms.add(bond.a1_id)
                used_atoms.add(bond.a2_id)

        external_double_atoms: set[int] = set()
        for bond in self.model.bonds.values():
            if bond.is_aromatic:
                continue
            if bond.order >= 2:
                external_double_atoms.add(bond.a1_id)
                external_double_atoms.add(bond.a2_id)

        for bond in aromatic_bonds:
            if component_atoms is not None and (
                bond.a1_id not in component_atoms or bond.a2_id not in component_atoms
            ):
                continue
            new_order = display_orders.get(bond.id, 1)
            if bond.a1_id in external_double_atoms or bond.a2_id in external_double_atoms:
                new_order = 1
            if bond.display_order != new_order or bond.order != new_order:
                bond.display_order = new_order
                bond.order = new_order
                self.update_bond_item(bond.id)

    def _build_aromatic_edges(
        self,
        vertex_defs: List[Tuple[Optional[int], float, float]],
        ring_size: int,
    ) -> List[Tuple[int, int, int, BondStyle, BondStereo, bool]]:
        if ring_size % 2 != 0:
            return [
                (i, (i + 1) % ring_size, 1, BondStyle.PLAIN, BondStereo.NONE, True)
                for i in range(ring_size)
            ]

        constraints: dict[int, int] = {}
        for i in range(ring_size):
            j = (i + 1) % ring_size
            a_id = vertex_defs[i][0]
            b_id = vertex_defs[j][0]
            if a_id is None or b_id is None:
                continue
            existing = self.model.find_bond_between(a_id, b_id)
            if existing is not None:
                constraints[i] = 2 if existing.order >= 2 else 1

        def atom_has_double(atom_id: Optional[int], other_id: Optional[int]) -> bool:
            if atom_id is None:
                return False
            for bond in self.model.bonds.values():
                if bond.order < 2:
                    continue
                if bond.a1_id == atom_id and bond.a2_id != other_id:
                    return True
                if bond.a2_id == atom_id and bond.a1_id != other_id:
                    return True
            return False

        def score_phase(phase: int) -> int:
            score = 0
            for idx, desired in constraints.items():
                current = 2 if ((idx + phase) % 2 == 0) else 1
                if current != desired:
                    score += 1
            return score

        phase = 0 if score_phase(0) <= score_phase(1) else 1
        edges: List[Tuple[int, int, int, BondStyle, BondStereo, bool]] = []
        for i in range(ring_size):
            j = (i + 1) % ring_size
            order = 2 if ((i + phase) % 2 == 0) else 1
            a_id = vertex_defs[i][0]
            b_id = vertex_defs[j][0]
            if order == 2 and (atom_has_double(a_id, b_id) or atom_has_double(b_id, a_id)):
                order = 1
            edges.append((i, j, order, BondStyle.PLAIN, BondStereo.NONE, True))
        return edges

    def _begin_drag(self, scene_pos: QPointF) -> None:
        if not self._selected_atom_ids_for_transform() and not self._selected_text_items():
            return
        self._dragging_selection = True
        self._drag_start_pos = scene_pos
        
        # Capture atoms
        atom_ids = self._selected_atom_ids_for_transform()
        self._drag_start_positions = {
            atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
            for atom_id in atom_ids
            if atom_id in self.model.atoms
        }
        
        # Capture text items
        self._drag_start_text_positions = {}
        for item in self._selected_text_items():
            self._drag_start_text_positions[item] = (item.pos(), item.rotation())
            
        self._drag_has_moved = False

    def _is_on_paper(self, x: float, y: float) -> bool:
        return (
            PAPER_MARGIN <= x <= self.paper_width - PAPER_MARGIN
            and PAPER_MARGIN <= y <= self.paper_height - PAPER_MARGIN
        )

    def _update_scene_rect(self) -> None:
        margin = 100
        self.scene.setSceneRect(
            -margin,
            -margin,
            self.paper_width + 2 * margin,
            self.paper_height + 2 * margin,
        )
