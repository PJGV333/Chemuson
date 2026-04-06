"""Clase pública del canvas compuesta por mixins de responsabilidad."""
from __future__ import annotations

from ._shared import (
    AromaticCircleItem,
    ArrowItem,
    AtomItem,
    BondItem,
    BracketItem,
    CHEMDOODLE_LIKE,
    ChemState,
    CompositeDiagramItem,
    DEFAULT_PAPER_HEIGHT,
    DEFAULT_PAPER_WIDTH,
    Dict,
    DrawingStyle,
    EnergyDiagramItem,
    GRID_MAJOR_STEP_PX,
    GRID_MINOR_STEP_PX,
    HoverAtomIndicatorItem,
    HoverBondIndicatorItem,
    ImageAnnotationItem,
    List,
    MolGraph,
    NUMBERING_AUTO_TEXT_ROLE,
    NUMBERING_KEY_ROLE,
    NUMBERING_KIND_ROLE,
    NUMBERING_TEXT_ROLE,
    NumberedStructure,
    OptimizeZoneItem,
    Optional,
    OrbitalAnnotationItem,
    PAPER_ITEM_ROLE,
    PreviewArrowItem,
    PreviewBondItem,
    PreviewChainItem,
    PreviewChainLabelItem,
    PreviewRingItem,
    QBrush,
    QColor,
    QFont,
    QGraphicsDropShadowEffect,
    QGraphicsEllipseItem,
    QGraphicsItem,
    QGraphicsLineItem,
    QGraphicsPathItem,
    QGraphicsRectItem,
    QGraphicsScene,
    QGraphicsTextItem,
    QGraphicsView,
    QPainter,
    QPainterPath,
    QPen,
    QPoint,
    QPointF,
    QRectF,
    QTimer,
    QUndoStack,
    Qt,
    RULER_MAJOR_STEP_PX,
    RULER_MINOR_STEP_PX,
    RULER_THICKNESS_PX,
    TextAnnotationItem,
    Tuple,
    WavyAnchorItem,
    compute_atom_numbers,
    compute_structure_numbers,
    math,
    pyqtSignal,
)
from .canvas_input import CanvasInputMixin
from .canvas_selection import CanvasSelectionMixin
from .canvas_text import CanvasTextMixin
from .canvas_render import CanvasRenderMixin
from .canvas_structure import CanvasStructureMixin

class ChemusonCanvas(
    CanvasInputMixin,
    CanvasSelectionMixin,
    CanvasTextMixin,
    CanvasRenderMixin,
    CanvasStructureMixin,
    QGraphicsView,
):
    """Lienzo principal para dibujar moleculas en paginas."""

    selection_changed = pyqtSignal(int, int, int, dict)

    def __init__(self, parent=None) -> None:
        """Inicializa la instancia.

        Args:
            parent: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        super().__init__(parent)

        self.scene = QGraphicsScene()
        self.setScene(self.scene)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.paper_width = DEFAULT_PAPER_WIDTH
        self.paper_height = DEFAULT_PAPER_HEIGHT
        self.model = MolGraph()
        self.state = ChemState()
        # Preferencias de nomenclatura avanzada (fase 4/6) por documento.
        self.name_advanced_enabled = True
        self.name_rdkit_isolated = True
        self.undo_stack = QUndoStack(self)
        self.drawing_style: DrawingStyle = CHEMDOODLE_LIKE
        self._ring_centers: dict[int, QPointF] = {}
        self._next_ring_id = 1
        self._electron_dots: set[TextAnnotationItem] = set()
        self._wavy_anchors: set[WavyAnchorItem] = set()

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
        self._flex_bond_pending: Optional[dict] = None
        self._arrow_start_pos: Optional[QPointF] = None
        self._arrow_end_pos: Optional[QPointF] = None
        self._arrow_curve_adjust_mode = False
        self._arrow_curve_factor = ArrowItem.CURVE_FACTOR_DEFAULT
        self.arrow_items: list[ArrowItem] = []
        self._bracket_drag_start: Optional[QPointF] = None
        self._bracket_preview: Optional[QGraphicsRectItem] = None
        self.bracket_items: list[BracketItem] = []
        self.energy_diagram_items: list[EnergyDiagramItem] = []
        self.plate_items: list = []
        self.semantic_diagram_items: list[CompositeDiagramItem] = []
        self.orbital_items: list[OrbitalAnnotationItem] = []
        self.image_items: list[ImageAnnotationItem] = []
        
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
        self._drag_start_text_positions: Dict[TextAnnotationItem, Tuple[QPointF, float]] = {}
        self._drag_start_arrow_positions: Dict[ArrowItem, Tuple[QPointF, QPointF]] = {}
        self._drag_start_bracket_rects: Dict[BracketItem, QRectF] = {}
        self._drag_start_energy_diagram_snapshots: Dict[EnergyDiagramItem, Tuple[QPointF, float, float, float]] = {}
        self._drag_start_semantic_diagram_snapshots: Dict[CompositeDiagramItem, Tuple[QPointF, float, float, float]] = {}
        self._drag_start_orbital_snapshots: Dict[OrbitalAnnotationItem, Tuple[QPointF, QPointF]] = {}
        self._drag_start_image_snapshots: Dict[ImageAnnotationItem, Tuple[QPointF, float, float, float]] = {}
        self._drag_start_selection_bbox: Optional[QRectF] = None
        self._drag_affected_bond_ids: set[int] = set()
        self._drag_affects_ring_centers = False
        self._drag_has_moved = False
        self._interaction_mouse_grabbed = False
        self._select_drag_mode: Optional[str] = None
        self._select_start_pos: Optional[QPointF] = None
        self._select_path: Optional[QPainterPath] = None
        self._select_preview_path: Optional[QGraphicsPathItem] = None
        self._select_preview_rect: Optional[QGraphicsRectItem] = None
        self._select_additive = False
        self._panning = False
        self._space_panning = False
        self._pan_last_pos: Optional[QPoint] = None

        self._drag_mode = "none"
        self._drag_anchor: Optional[dict] = None
        self._ring_last_vertices: Optional[List[QPointF]] = None
        self._chain_last_points: Optional[List[QPointF]] = None
        self._pending_template_graph: Optional[MolGraph] = None
        self._pending_template_label: Optional[str] = None
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
        self._rotation_start_arrow_positions: Dict[ArrowItem, Tuple[QPointF, QPointF]] = {}
        self._rotation_start_energy_diagram_snapshots: Dict[EnergyDiagramItem, Tuple[QPointF, float, float, float]] = {}
        self._rotation_start_semantic_diagram_snapshots: Dict[CompositeDiagramItem, Tuple[QPointF, float, float, float]] = {}
        self._rotation_start_orbital_snapshots: Dict[OrbitalAnnotationItem, Tuple[QPointF, QPointF]] = {}
        self._rotation_start_image_snapshots: Dict[ImageAnnotationItem, Tuple[QPointF, float, float, float]] = {}
        self._is_rotating_3d = False
        self._rotation_3d_start_view_pos: Optional[QPoint] = None
        self._rotation_3d_before_positions: Dict[int, Tuple[float, float]] = {}
        self._rotation_3d_click_atom_id: Optional[int] = None
        self._rotation_3d_has_moved = False
        self._rotation_3d_ref_atom_ids: tuple[int, ...] = tuple()
        self._rotation_3d_ref_positions: Dict[int, Tuple[float, float]] = {}
        self._rotation_3d_pitch_deg = 0.0
        self._rotation_3d_yaw_deg = 0.0
        self._rotation_3d_drag_start_pitch_deg = 0.0
        self._rotation_3d_drag_start_yaw_deg = 0.0
        self._fragment_pivot_atom_id: Optional[int] = None
        self._scale_dragging = False
        self._scale_anchor: Optional[QPointF] = None
        self._scale_start_handle: Optional[QPointF] = None
        self._scale_start_length = 0.0
        self._scale_start_positions: Dict[int, Tuple[float, float]] = {}
        self._scale_start_text_positions: Dict[TextAnnotationItem, Tuple[QPointF, float]] = {}
        self._scale_start_text_styles: Dict[TextAnnotationItem, Tuple[str, float]] = {}
        self._scale_start_arrow_positions: Dict[ArrowItem, Tuple[QPointF, QPointF]] = {}
        self._scale_start_arrow_strokes: Dict[ArrowItem, Optional[float]] = {}
        self._scale_start_bracket_rects: Dict[BracketItem, QRectF] = {}
        self._scale_start_bracket_styles: Dict[BracketItem, Tuple[float, Optional[float]]] = {}
        self._scale_start_energy_diagram_snapshots: Dict[EnergyDiagramItem, Tuple[QPointF, float, float, float]] = {}
        self._scale_start_semantic_diagram_snapshots: Dict[CompositeDiagramItem, Tuple[QPointF, float, float, float]] = {}
        self._scale_start_orbital_snapshots: Dict[OrbitalAnnotationItem, Tuple[QPointF, QPointF]] = {}
        self._scale_start_image_snapshots: Dict[ImageAnnotationItem, Tuple[QPointF, float, float, float]] = {}
        self._scale_start_atom_label_scales: Dict[int, Optional[float]] = {}
        self._scale_start_atom_sphere_radii: Dict[int, Tuple[Optional[float], float]] = {}
        self._scale_start_bond_strokes: Dict[int, Tuple[Optional[float], float]] = {}
        self._scale_start_bond_lengths: Dict[int, Optional[float]] = {}
        self._scale_start_text_effective_widths: Dict[TextAnnotationItem, float] = {}
        self._scale_text_reflow_only = False
        self._scale_has_moved = False
        self._implicit_h_overlays: Dict[int, list[tuple[QGraphicsLineItem, QGraphicsTextItem]]] = {}
        self._group_anchor_overrides: Dict[int, str] = {}
        self._numbering_overlay_items: list[TextAnnotationItem] = []
        self._numbering_atom_items: dict[int, TextAnnotationItem] = {}
        self._numbering_structure_items: dict[int, TextAnnotationItem] = {}
        self._numbering_atom_base_pos: dict[int, tuple[float, float]] = {}
        self._numbering_structure_base_pos: dict[int, tuple[float, float]] = {}
        self._atom_numbering: dict[int, int] = {}
        self._structure_numbering: list[NumberedStructure] = []
        self._suspend_numbering_refresh = False
        self._validation_batch_depth = 0
        self._validation_dirty = False
        self._force_full_scene_sync_on_undo_change = False

        self._setup_view()
        self._pending_initial_center = True
        self._create_paper()
        self._create_overlays()

        self.scene.selectionChanged.connect(self._sync_selection_from_scene)
        self.undo_stack.indexChanged.connect(self._on_undo_stack_changed)
        self.setMouseTracking(True)
        self.viewport().setMouseTracking(True)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)

    @property
    def graph(self) -> MolGraph:
        """Alias for model for compatibility."""
        return self.model

    def _setup_view(self) -> None:
        """Método auxiliar para  setup view.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self.setBackgroundBrush(QBrush(QColor("#E0E0E0")))
        self.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.setRenderHint(QPainter.RenderHint.TextAntialiasing)

        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)

        self._update_scene_rect()

        self.setDragMode(QGraphicsView.DragMode.NoDrag)

        self._zoom_factor = 1.0
        self._min_zoom = 0.25
        self._max_zoom = 4.0

    def set_show_rulers(self, enabled: bool) -> None:
        """Actualiza el estado de show rulers.

        Args:
            enabled: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self.show_rulers = enabled
        self.viewport().update()

    def set_show_grid(self, enabled: bool) -> None:
        """Actualiza el estado de show grid.

        Args:
            enabled: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self.show_grid = enabled
        if self._grid_minor_item is not None:
            self._grid_minor_item.setVisible(enabled)
        if self._grid_major_item is not None:
            self._grid_major_item.setVisible(enabled)
        self.viewport().update()

    @staticmethod
    def _normalize_numbering_mode(mode: str) -> str:
        """Normaliza el modo de numeración a un valor soportado."""
        value = str(mode or "").strip().lower()
        if value in {"atoms", "structures", "both"}:
            return value
        return "atoms"

    def set_numbering_enabled(self, enabled: bool) -> None:
        """Activa/desactiva la numeración visual del lienzo."""
        self.state.numbering_enabled = bool(enabled)
        self.recompute_numbering()

    def set_numbering_mode(self, mode: str) -> None:
        """Configura modo de numeración: atoms, structures o both."""
        self.state.numbering_mode = self._normalize_numbering_mode(mode)
        self.recompute_numbering()

    def set_numbering_include_in_export(self, enabled: bool) -> None:
        """Define si la numeración se incluye en exportaciones."""
        self.state.numbering_include_in_export = bool(enabled)

    def set_numbering_style(
        self,
        *,
        font_size: Optional[float] = None,
        offset_x: Optional[float] = None,
        offset_y: Optional[float] = None,
        circle: Optional[bool] = None,
        color: Optional[str] = None,
        background: Optional[str] = None,
    ) -> None:
        """Actualiza estilo visual de numeración y refresca overlay."""
        if font_size is not None:
            self.state.numbering_font_size = max(1.0, float(font_size))
        if offset_x is not None:
            self.state.numbering_offset_x = float(offset_x)
        if offset_y is not None:
            self.state.numbering_offset_y = float(offset_y)
        if circle is not None:
            self.state.numbering_circle = bool(circle)
        if color:
            self.state.numbering_color = str(color)
        if background:
            self.state.numbering_background = str(background)
        # Si se actualiza estilo explícitamente, aplicarlo en todos los números existentes.
        if font_size is not None or color:
            font = self._numbering_font()
            color_value = QColor(self.state.numbering_color)
            for item in self._numbering_overlay_items:
                if item is None or item.scene() is not self.scene:
                    continue
                if font_size is not None:
                    item.setFont(font)
                if color:
                    item.setDefaultTextColor(color_value)
        self.recompute_numbering(force_reset=False)

    def recompute_numbering(self, force_reset: bool = False) -> None:
        """Recalcula y redibuja overlays de numeración."""
        if self._suspend_numbering_refresh:
            return
        self._cleanup_numbering_overlays()
        if not self.state.numbering_enabled or not self.model.atoms:
            self._clear_numbering_overlays()
            self._atom_numbering = {}
            self._structure_numbering = []
            return

        font = self._numbering_font()
        color = QColor(self.state.numbering_color)
        mode = self._normalize_numbering_mode(self.state.numbering_mode)
        self.state.numbering_mode = mode

        required_atom_ids: set[int] = set()
        required_structure_ids: set[int] = set()

        if mode in {"atoms", "both"}:
            self._atom_numbering = compute_atom_numbers(self.model)
            off_x = float(self.state.numbering_offset_x)
            off_y = float(self.state.numbering_offset_y)
            for atom_id, number in sorted(self._atom_numbering.items(), key=lambda item: item[0]):
                atom = self.model.atoms.get(atom_id)
                if atom is None:
                    continue
                required_atom_ids.add(int(atom_id))
                item = self._numbering_atom_items.get(int(atom_id))
                if item is None or item.scene() is not self.scene:
                    item = self._create_numbering_item(kind="atom", key=int(atom_id), auto_text=str(number))
                    self._numbering_atom_items[int(atom_id)] = item
                    self._numbering_overlay_items.append(item)
                self._update_numbering_item(
                    item=item,
                    auto_text=str(number),
                    anchor_x=float(atom.x) + off_x,
                    anchor_y=float(atom.y) + off_y,
                    base_positions=self._numbering_atom_base_pos,
                    key=int(atom_id),
                    font=font,
                    color=color,
                    force_reset=force_reset,
                )
        else:
            self._atom_numbering = {}

        if mode in {"structures", "both"}:
            self._structure_numbering = compute_structure_numbers(self.model)
            for numbered_structure in self._structure_numbering:
                structure_id = int(numbered_structure.structure_id)
                required_structure_ids.add(structure_id)
                item = self._numbering_structure_items.get(structure_id)
                auto_text = f"({numbered_structure.number})"
                if item is None or item.scene() is not self.scene:
                    item = self._create_numbering_item(
                        kind="structure",
                        key=structure_id,
                        auto_text=auto_text,
                    )
                    self._numbering_structure_items[structure_id] = item
                    self._numbering_overlay_items.append(item)
                self._update_numbering_item(
                    item=item,
                    auto_text=auto_text,
                    anchor_x=float(numbered_structure.x),
                    anchor_y=float(numbered_structure.y),
                    base_positions=self._numbering_structure_base_pos,
                    key=structure_id,
                    font=font,
                    color=color,
                    force_reset=force_reset,
                )
        else:
            self._structure_numbering = []

        for atom_id in list(self._numbering_atom_items.keys()):
            if atom_id not in required_atom_ids:
                self._drop_numbering_item(self._numbering_atom_items.pop(atom_id, None))
                self._numbering_atom_base_pos.pop(atom_id, None)
        for structure_id in list(self._numbering_structure_items.keys()):
            if structure_id not in required_structure_ids:
                self._drop_numbering_item(self._numbering_structure_items.pop(structure_id, None))
                self._numbering_structure_base_pos.pop(structure_id, None)

    def _numbering_font(self) -> QFont:
        """Fuente por defecto para números de numeración."""
        font = QFont(self.state.label_font_family)
        font_size = float(self.state.numbering_font_size)
        if font_size <= 0.0:
            font_size = 10.0
        font.setPointSizeF(font_size)
        font.setBold(False)
        font.setItalic(False)
        font.setUnderline(False)
        return font

    def _create_numbering_item(self, *, kind: str, key: int, auto_text: str) -> TextAnnotationItem:
        """Crea un item de texto para numeración editable/movible."""
        item = TextAnnotationItem(auto_text, 0.0, 0.0)
        item.setData(NUMBERING_TEXT_ROLE, True)
        item.setData(NUMBERING_KIND_ROLE, str(kind))
        item.setData(NUMBERING_KEY_ROLE, int(key))
        item.setData(NUMBERING_AUTO_TEXT_ROLE, str(auto_text))
        item.setFont(self._numbering_font())
        item.setDefaultTextColor(QColor(self.state.numbering_color))
        item.setOpacity(self.canvas_default_opacity())
        item.setZValue(35)
        self.scene.addItem(item)
        return item

    def _update_numbering_item(
        self,
        *,
        item: TextAnnotationItem,
        auto_text: str,
        anchor_x: float,
        anchor_y: float,
        base_positions: dict[int, tuple[float, float]],
        key: int,
        font: QFont,
        color: QColor,
        force_reset: bool,
    ) -> None:
        """Actualiza texto y posición de un item de numeración."""
        previous_auto = str(item.data(NUMBERING_AUTO_TEXT_ROLE) or "")
        current_text = item.toPlainText().strip()

        # Mantiene texto manual del usuario, salvo recalculo forzado.
        if force_reset or not current_text or current_text == previous_auto:
            if item.toPlainText() != auto_text:
                item.setPlainText(auto_text)
        item.setData(NUMBERING_AUTO_TEXT_ROLE, str(auto_text))

        if force_reset:
            item.setFont(font)
            item.setDefaultTextColor(color)

        rect = item.boundingRect()
        base_x = float(anchor_x) - rect.width() * 0.5
        base_y = float(anchor_y) - rect.height() * 0.5

        if force_reset:
            offset_x = 0.0
            offset_y = 0.0
        else:
            previous_base = base_positions.get(int(key))
            if previous_base is None:
                offset_x = 0.0
                offset_y = 0.0
            else:
                offset_x = float(item.pos().x()) - float(previous_base[0])
                offset_y = float(item.pos().y()) - float(previous_base[1])

        item.setPos(base_x + offset_x, base_y + offset_y)
        base_positions[int(key)] = (base_x, base_y)

    def _drop_numbering_item(self, item: Optional[TextAnnotationItem]) -> None:
        """Elimina un item de numeración de la escena y colecciones."""
        if item is None:
            return
        if item in self._numbering_overlay_items:
            self._numbering_overlay_items.remove(item)
        if item.scene() is self.scene:
            self.scene.removeItem(item)

    def _cleanup_numbering_overlays(self) -> None:
        """Limpia referencias a items de numeración borrados externamente."""
        self._numbering_overlay_items = [
            item
            for item in self._numbering_overlay_items
            if item is not None and item.scene() is self.scene
        ]
        for atom_id, item in list(self._numbering_atom_items.items()):
            if item is None or item.scene() is not self.scene:
                self._numbering_atom_items.pop(atom_id, None)
                self._numbering_atom_base_pos.pop(atom_id, None)
        for structure_id, item in list(self._numbering_structure_items.items()):
            if item is None or item.scene() is not self.scene:
                self._numbering_structure_items.pop(structure_id, None)
                self._numbering_structure_base_pos.pop(structure_id, None)

    def _clear_numbering_overlays(self) -> None:
        """Elimina todos los overlays de numeración."""
        for item in list(self._numbering_overlay_items):
            self._drop_numbering_item(item)
        self._numbering_overlay_items.clear()
        self._numbering_atom_items.clear()
        self._numbering_structure_items.clear()
        self._numbering_atom_base_pos.clear()
        self._numbering_structure_base_pos.clear()

    def paintEvent(self, event) -> None:
        """Método auxiliar para paintEvent.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        super().paintEvent(event)
        if not self.show_rulers:
            return
        painter = QPainter(self.viewport())
        painter.setRenderHint(QPainter.RenderHint.TextAntialiasing)
        if self.show_rulers:
            self._paint_rulers(painter)
        painter.end()

    def _paint_rulers(self, painter: QPainter) -> None:
        """Método auxiliar para  paint rulers.

        Args:
            painter: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para leaveEvent.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        super().leaveEvent(event)

    def _create_paper(self) -> None:
        """Método auxiliar para  create paper.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        for item in self._paper_scene_items():
            self.scene.removeItem(item)
        self.paper = QGraphicsRectItem(0, 0, self.paper_width, self.paper_height)
        self.paper.setBrush(QBrush(Qt.GlobalColor.white))
        self.paper.setPen(QPen(QColor("#CCCCCC"), 1))
        self.paper.setZValue(-10)
        self.paper.setData(PAPER_ITEM_ROLE, True)

        shadow = QGraphicsDropShadowEffect()
        shadow.setBlurRadius(20)
        shadow.setColor(QColor(0, 0, 0, 80))
        shadow.setOffset(5, 5)
        self.paper.setGraphicsEffect(shadow)

        self.scene.addItem(self.paper)
        self._create_grid()
        self.center_on_paper()

    def _is_paper_scene_item(self, item: QGraphicsItem) -> bool:
        """Detecta la hoja de trabajo, incluyendo duplicados huérfanos legacy."""
        if item is None or item.scene() is not self.scene:
            return False
        try:
            if bool(item.data(PAPER_ITEM_ROLE)):
                return True
        except Exception:
            return False
        if not isinstance(item, QGraphicsRectItem):
            return False
        if item.parentItem() is not None or item.zValue() > -9.5:
            return False
        rect = item.rect()
        if (
            abs(float(rect.x())) > 0.1
            or abs(float(rect.y())) > 0.1
            or rect.width() <= 0.0
            or rect.height() <= 0.0
        ):
            return False
        try:
            brush = item.brush()
            pen = item.pen()
            brush_color = brush.color().name().lower()
            pen_color = pen.color().name().lower()
        except Exception:
            return False
        return brush_color == "#ffffff" and pen_color == "#cccccc"

    def _paper_scene_items(self) -> list[QGraphicsRectItem]:
        """Devuelve todas las hojas visibles/huérfanas presentes en la escena."""
        return [
            item
            for item in self.scene.items()
            if isinstance(item, QGraphicsRectItem) and self._is_paper_scene_item(item)
        ]

    def _create_grid(self) -> None:
        """Método auxiliar para  create grid.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  create overlays.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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

    def center_on_paper(self) -> None:
        """Método auxiliar para center on paper.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self.centerOn(self.paper_width / 2, self.paper_height / 2)

    def set_paper_size(self, width: int, height: int) -> None:
        """Actualiza el estado de paper size.

        Args:
            width: Descripción del parámetro.
            height: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        self.center_on_paper()
        self._update_selection_overlay()
        self.viewport().update()

    def begin_validation_batch(self) -> None:
        """Agrupa validaciones sucesivas para ejecutarlas una sola vez al final."""
        self._validation_batch_depth += 1

    def end_validation_batch(self) -> None:
        """Cierra un lote de validación y ejecuta la validación pendiente si toca."""
        if self._validation_batch_depth <= 0:
            self._validation_batch_depth = 0
            return
        self._validation_batch_depth -= 1
        if self._validation_batch_depth == 0 and self._validation_dirty:
            self._validation_dirty = False
            self.validate_structure()

    def request_structure_validation(self) -> list[int]:
        """Valida ahora o difiere la validación si se está en un batch."""
        if self._validation_batch_depth > 0:
            self._validation_dirty = True
            return []
        self._validation_dirty = False
        return self.validate_structure()

    def request_full_scene_sync(self) -> None:
        """Marca que el siguiente cambio de undo stack debe resincronizar toda la escena."""
        self._force_full_scene_sync_on_undo_change = True

    def zoom_in(self) -> None:
        """Método auxiliar para zoom in.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._zoom_factor < self._max_zoom:
            self._zoom_factor *= 1.2
            self.scale(1.2, 1.2)
            self._update_scene_rect()
            self._update_selection_overlay()

    def zoom_out(self) -> None:
        """Método auxiliar para zoom out.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._zoom_factor > self._min_zoom:
            self._zoom_factor /= 1.2
            self.scale(1 / 1.2, 1 / 1.2)
            self._update_scene_rect()
            self._update_selection_overlay()

    def _show_status_message(self, message: str, timeout_ms: int = 7000) -> None:
        """Muestra un mensaje transitorio si la ventana padre expone barra de estado."""
        if not message:
            return
        try:
            window = self.window()
        except RuntimeError:
            return
        if window is None or not hasattr(window, "statusBar"):
            return
        try:
            status_bar = window.statusBar()
        except Exception:
            return
        if status_bar is None:
            return
        try:
            status_bar.showMessage(str(message), int(timeout_ms))
        except Exception:
            return

    def _update_scene_rect(self) -> None:
        # NOTE: The sceneRect drives scroll bar ranges. If we compute its margin
        # purely in view pixels, zooming out makes the visible scene area
        # (in scene units) much larger than the margin, which can clamp the view
        # and make parts of the paper unreachable after resizes.
        """Método auxiliar para  update scene rect.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        viewport_rect = self.viewport().rect()
        if viewport_rect.isEmpty():
            view_w = float(self.viewport().width())
            view_h = float(self.viewport().height())
        else:
            view_poly = self.mapToScene(viewport_rect)
            view_bounds = view_poly.boundingRect()
            view_w = float(view_bounds.width())
            view_h = float(view_bounds.height())

        margin = max(100.0, view_w, view_h)

        # Preserve current view center so changing the scene rect doesn't jump
        # the document to a corner.
        current_center = self.mapToScene(viewport_rect.center()) if not viewport_rect.isEmpty() else None
        self.scene.setSceneRect(
            -margin,
            -margin,
            float(self.paper_width) + 2.0 * margin,
            float(self.paper_height) + 2.0 * margin,
        )
        if current_center is not None:
            self.centerOn(current_center)

    def resizeEvent(self, event) -> None:
        """Método auxiliar para resizeEvent.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        super().resizeEvent(event)
        self._update_scene_rect()
        self._update_selection_overlay()
        # Compare in scene units so zoom level is respected.
        viewport_rect = self.viewport().rect()
        if viewport_rect.isEmpty():
            return
        view_bounds = self.mapToScene(viewport_rect).boundingRect()
        if view_bounds.width() >= self.paper_width and view_bounds.height() >= self.paper_height:
            self.center_on_paper()

    def showEvent(self, event) -> None:
        """Método auxiliar para showEvent.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        super().showEvent(event)
        if not getattr(self, "_pending_initial_center", False):
            return
        self._pending_initial_center = False
        QTimer.singleShot(0, self._center_canvas_initial)

    def _center_canvas_initial(self) -> None:
        """Método auxiliar para  center canvas initial.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self.viewport().rect().isEmpty():
            self._pending_initial_center = True
            return
        self._update_scene_rect()
        self.center_on_paper()
        self._update_selection_overlay()
