"""
Lienzo principal de Chemuson (QGraphicsView/QGraphicsScene).

Gestiona el dibujo de moléculas, herramientas de edición, anotaciones,
persistencia visual y renderización para exportaciones.
"""
from __future__ import annotations

import base64
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
    QDialog,
    QInputDialog,
    QColorDialog,
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
from PyQt6.QtCore import Qt, QPoint, QPointF, QRectF, QRect, QSize, QBuffer, QMimeData, QTimer, pyqtSignal

try:
    from PyQt6.QtSvg import QSvgGenerator
except Exception:  # Optional Qt module at runtime
    QSvgGenerator = None

from core.model import Bond, BondStyle, BondStereo, ChemState, MolGraph
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
    WavyAnchorItem,
    TextAnnotationItem,
    ABBREVIATION_LABELS,
)
from gui.style import CHEMDOODLE_LIKE, DrawingStyle
from gui.numbering import NumberedStructure, compute_atom_numbers, compute_structure_numbers
from gui.templates import build_haworth_template, build_chair_template
from gui.geom import (
    angle_deg,
    snap_angle_deg,
    endpoint_from_angle_len,
    project_3d_rotation,
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
    AddTextItemCommand,
    AddWavyAnchorCommand,
    ChangeAtomCommand,
    ChangeBondCommand,
    ChangeBondLengthCommand,
    ChangeBondStrokeCommand,
    ChangeBondColorCommand,
    ChangeChargeCommand,
    ChangeNoImplicitCommand,
    ChangeCoordinationSphereStyleCommand,
    SetCoordinationCenterCommand,
    DeleteSelectionCommand,
    MoveAtomsCommand,
    MoveTextItemsCommand,
    MoveArrowItemsCommand,
    MoveBracketItemsCommand,
)
from gui.dialogs import AtomLabelDialog, TrackballRotationDialog
from chemname.molview import MolView
from chemname.rings import find_rings_simple, ring_bonds
from chemio.rdkit_io import (
    molgraph_to_molfile,
    molgraph_to_smiles,
    molfile_to_molgraph,
    smiles_to_molgraph,
    kekulize_display_orders,
)
from chemname import iupac_name, NameOptions


# Dimensiones de papel (aprox. A4 en px)
DEFAULT_PAPER_WIDTH = 800
DEFAULT_PAPER_HEIGHT = 1000
PAPER_MARGIN = 40  # Márgenes donde no se deben colocar átomos
ATOM_HIT_RADIUS = 16.0
DEFAULT_BOND_LENGTH = 40.0
HOVER_ATOM_RADIUS = 16.0
HOVER_BOND_DISTANCE = 10.0
OPTIMIZE_ZONE_SCALE = 1.2
# Límite superior de segmentos para herramienta de cadena (zig-zag).
# Se mantiene un tope alto para evitar recortes prematuros en cadenas largas.
CHAIN_MAX_BONDS = 120
ANGLE_OCCUPIED_TOLERANCE_DEG = 20.0
# Tolerancia específica para sp3: evita descartar candidatos tetraédricos.
SP3_OCCUPIED_TOLERANCE_DEG = 8.0
MIN_ATOM_DIST_SCALE = 0.65
MIN_BOND_DIST_SCALE = 0.2
COLLISION_LENGTH_BOOST = 1.2
BOND_OVERLAP_TOLERANCE_PX = 5.0
BOND_DRAG_THRESHOLD_PX = 6.0
BOND_LAST_ANGLE_TOLERANCE_DEG = 20.0
CLIPBOARD_RENDER_SCALE = 5.0
TRACKBALL_ROTATION_DEG_PER_PIXEL = 1.0
TRACKBALL_MAX_TILT_DEG = 60.0
TRACKBALL_REFERENCE_MATCH_TOLERANCE_PX = 1.5
# Etiquetas rápidas para grupos funcionales comunes.
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
ANALYSIS_MARGIN_PX = 14.0
ANALYSIS_MIN_PEAK_PERCENT = 1.0
ANALYSIS_MAX_PEAKS = 10
ANALYSIS_DIST_KEEP = 400

ATOMIC_WEIGHTS = {
    "H": 1.00794,
    "C": 12.0107,
    "N": 14.0067,
    "O": 15.9994,
    "S": 32.065,
    "P": 30.973762,
    "F": 18.998403,
    "Cl": 35.453,
    "Br": 79.904,
    "I": 126.90447,
    "Si": 28.0855,
    "B": 10.811,
    "Se": 78.96,
    "Li": 6.941,
    "Na": 22.98977,
    "K": 39.0983,
    "Mg": 24.305,
}

MONOISOTOPIC_MASSES = {
    "H": 1.007825032,
    "C": 12.0,
    "N": 14.003074005,
    "O": 15.99491462,
    "S": 31.972071,
    "P": 30.973761998,
    "F": 18.998403163,
    "Cl": 34.96885268,
    "Br": 78.9183376,
    "I": 126.904468,
    "Si": 27.9769265,
    "B": 11.009305,
    "Se": 79.916521,
    "Li": 7.016003,
    "Na": 22.98976928,
    "K": 38.96370668,
    "Mg": 23.9850417,
}

ISOTOPE_ABUNDANCES = {
    "H": [(1.007825032, 0.999885), (2.014101778, 0.000115)],
    "C": [(12.0, 0.9893), (13.003354835, 0.0107)],
    "N": [(14.003074005, 0.99632), (15.000108898, 0.00368)],
    "O": [
        (15.99491462, 0.99757),
        (16.999131757, 0.00038),
        (17.999159613, 0.00205),
    ],
    "S": [
        (31.972071, 0.9499),
        (32.971458, 0.0075),
        (33.967867, 0.0425),
        (35.967081, 0.0001),
    ],
    "P": [(30.973761998, 1.0)],
    "F": [(18.998403163, 1.0)],
    "Cl": [(34.96885268, 0.7576), (36.96590260, 0.2424)],
    "Br": [(78.9183376, 0.5069), (80.916291, 0.4931)],
    "I": [(126.904468, 1.0)],
    "Si": [(27.9769265, 0.92223), (28.9764947, 0.04685), (29.9737701, 0.03092)],
    "B": [(10.012937, 0.199), (11.009305, 0.801)],
    "Se": [
        (73.922476, 0.0089),
        (75.919214, 0.0937),
        (76.919915, 0.0763),
        (77.917309, 0.2377),
        (79.916522, 0.4961),
        (81.916700, 0.0873),
    ],
    "Li": [(6.015123, 0.075), (7.016004, 0.925)],
    "Na": [(22.98976928, 1.0)],
    "K": [(38.96370668, 0.9326), (39.96399848, 0.000117), (40.96182576, 0.0673)],
    "Mg": [(23.9850417, 0.7899), (24.9858369, 0.10), (25.9825929, 0.1101)],
}
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
SELECTION_ROTATE_OFFSET_PX = 13.0
SELECTION_MOVE_OFFSET_PX = 0.0
GRID_MINOR_STEP_PX = 20
GRID_MAJOR_STEP_PX = 100
SELECTION_BOX_PADDING_PX = 6.0
SELECTION_HANDLE_RADIUS_PX = 6.0
SELECTION_HANDLE_OFFSET_PX = 18.0

HETERO_ELECTRON_ATOMS = {"N", "O", "S", "P", "F", "Cl", "Br", "I", "Se", "B"}
ELECTRON_SLOT_ANGLES = [0, 45, 90, 135, 180, 225, 270, 315]
ELECTRON_SLOT_TOLERANCE_DEG = 18.0
ELECTRON_PAIR_SPREAD_PX = 6.5
ELECTRON_DOT_ROLE = 1001
ELECTRON_ANCHOR_ROLE = 1002
ELECTRON_SLOT_ROLE = 1003
ELECTRON_SIDE_ROLE = 1004
ELECTRON_SCALE_ROLE = 1005
WAVY_ANCHOR_ROLE = 2001
WAVY_ANCHOR_ANGLE_ROLE = 2002
WAVY_ANCHOR_LENGTH_ROLE = 2003
WAVY_ANCHOR_BOND_ROLE = 2004
AROMATIC_CIRCLE_ATOMS_ROLE = 3001
NUMBERING_TEXT_ROLE = 4001
NUMBERING_KIND_ROLE = 4002
NUMBERING_KEY_ROLE = 4003
NUMBERING_AUTO_TEXT_ROLE = 4004

SYMBOL_TEXT_TOOLS = {
    "tool_symbol_plus": {"text": "+", "scale": 1.0, "anchor": True},
    "tool_symbol_minus": {"text": "-", "scale": 1.0, "anchor": True},
    "tool_symbol_radical": {"text": "•", "scale": 1.2, "anchor": True, "rotate": False, "electrons": 1},
    "tool_symbol_lone_pair": {"text": "··", "scale": 1.1, "anchor": True, "rotate": True, "electrons": 2},
    "tool_symbol_radical_cation": {"text": "•+", "scale": 1.1, "anchor": True, "rotate": False},
    "tool_symbol_radical_anion": {"text": "•-", "scale": 1.1, "anchor": True, "rotate": False},
    "tool_symbol_partial_plus": {"text": "δ+", "scale": 1.0, "anchor": True, "rotate": False},
    "tool_symbol_partial_minus": {"text": "δ-", "scale": 1.0, "anchor": True, "rotate": False},
}


class ChemusonCanvas(QGraphicsView):
    """
    Lienzo principal para dibujar moléculas en páginas.
    Usa QGraphicsView/QGraphicsScene con una hoja centrada.
    """
    
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
        self._drag_has_moved = False
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
        self._scale_dragging = False
        self._scale_anchor: Optional[QPointF] = None
        self._scale_start_handle: Optional[QPointF] = None
        self._scale_start_length = 0.0
        self._scale_start_positions: Dict[int, Tuple[float, float]] = {}
        self._scale_start_text_positions: Dict[TextAnnotationItem, Tuple[QPointF, float]] = {}
        self._scale_start_arrow_positions: Dict[ArrowItem, Tuple[QPointF, QPointF]] = {}
        self._scale_start_bracket_rects: Dict[BracketItem, QRectF] = {}
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

    def set_current_tool(self, tool_id: str) -> None:
        """Actualiza el estado de current tool.

        Args:
            tool_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        tool_id = tool_id or "tool_none"
        if tool_id.startswith("atom_"):
            self.state.default_element = tool_id.split("_", 1)[1]
            tool_id = "tool_atom"
        elif tool_id.startswith("coord_"):
            tool_id = "tool_coordination_center"
        elif tool_id.startswith("bond_"):
            tool_id = "tool_bond"
        elif tool_id == "tool_coordination_bond":
            self.state.active_bond_order = 1
            self.state.active_bond_style = BondStyle.COORDINATION
            self.state.active_bond_stereo = BondStereo.NONE
            self.state.active_bond_mode = "set"
            self.state.active_bond_aromatic = False
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
        if self._pending_template_graph is not None:
            self.cancel_template_insert_mode()
        self._clear_bond_anchor()
        self._cancel_drag()
        self._arrow_start_pos = None
        self._arrow_end_pos = None
        self._arrow_curve_adjust_mode = False
        self._arrow_curve_factor = ArrowItem.CURVE_FACTOR_DEFAULT
        if hasattr(self, "_preview_arrow_item") and self._preview_arrow_item is not None:
            self._preview_arrow_item.hide_preview()
        self._clear_bracket_preview()

        if tool_id in {"tool_select", "tool_select_lasso", "tool_none"}:
            self.setCursor(Qt.CursorShape.ArrowCursor)
            self.setDragMode(QGraphicsView.DragMode.NoDrag)
        else:
            self.setCursor(Qt.CursorShape.CrossCursor)
            self.setDragMode(QGraphicsView.DragMode.NoDrag)

    def set_active_bond(self, bond_spec: dict) -> None:
        """Actualiza el estado de active bond.

        Args:
            bond_spec: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Actualiza el estado de active ring.

        Args:
            ring_spec: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self.state.active_ring_size = ring_spec.get("size", 6)
        self.state.active_ring_aromatic = ring_spec.get("aromatic", False)
        self.state.active_ring_template = ring_spec.get("template")
        self.state.active_ring_anomeric = ring_spec.get("anomeric")
        self.set_current_tool("tool_ring")

    def set_active_element(self, element: str) -> None:
        """Actualiza el estado de active element.

        Args:
            element: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self.state.default_element = element
        self.set_current_tool("tool_atom")

    @staticmethod
    def _curve_factor_from_pointer(start: QPointF, end: QPointF, pointer: QPointF) -> float:
        """Calcula curvatura normalizada a partir de la distancia perpendicular."""
        dx = end.x() - start.x()
        dy = end.y() - start.y()
        length = math.hypot(dx, dy)
        if length <= 1e-6:
            return ArrowItem.CURVE_FACTOR_DEFAULT
        nx = -dy / length
        ny = dx / length
        mid_x = (start.x() + end.x()) * 0.5
        mid_y = (start.y() + end.y()) * 0.5
        offset = (pointer.x() - mid_x) * nx + (pointer.y() - mid_y) * ny
        # Higher sensitivity so users can reach strong curvature comfortably.
        factor = offset / (length * 0.65)
        return max(ArrowItem.CURVE_FACTOR_MIN, min(ArrowItem.CURVE_FACTOR_MAX, factor))

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

    def mousePressEvent(self, event) -> None:
        """Método auxiliar para mousePressEvent.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if event.button() == Qt.MouseButton.MiddleButton or (
            self._space_panning and event.button() == Qt.MouseButton.LeftButton
        ):
            self._panning = True
            self._pan_last_pos = event.pos()
            self.setCursor(Qt.CursorShape.ClosedHandCursor)
            event.accept()
            return
        scene_pos = self.mapToScene(event.pos())
        self._last_scene_pos = scene_pos

        if self._drag_mode == "flex_adjust":
            if event.button() == Qt.MouseButton.LeftButton:
                self._commit_flex_bond()
                return
            if event.button() == Qt.MouseButton.RightButton:
                self._cancel_drag()
                return

        if not self._is_on_paper(scene_pos.x(), scene_pos.y()) and self.current_tool not in {"tool_select", "tool_select_lasso"}:
            super().mousePressEvent(event)
            return

        if event.button() == Qt.MouseButton.RightButton:
            clicked_item = self._get_item_at(scene_pos)
            if isinstance(
                clicked_item,
                (AtomItem, BondItem, ArrowItem, BracketItem, TextAnnotationItem),
            ):
                if not (event.modifiers() & (Qt.KeyboardModifier.ShiftModifier | Qt.KeyboardModifier.ControlModifier)):
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

        if self.current_tool == "tool_none":
            return

        if self._pending_template_graph is not None:
            clicked_atom_id, _clicked_bond_id = self._pick_hover_target(scene_pos)
            if not self._is_on_paper(scene_pos.x(), scene_pos.y()):
                return
            self._insert_molgraph_at(
                self._pending_template_graph,
                scene_pos,
                attach_to_atom_id=clicked_atom_id,
            )
            self.cancel_template_insert_mode()
            return

        if self.current_tool in {"tool_select", "tool_select_lasso"}:
            if self._space_panning:
                return
            clicked_item = self._get_item_at(scene_pos)
            blocked_modifiers = (
                Qt.KeyboardModifier.ShiftModifier
                | Qt.KeyboardModifier.ControlModifier
                | Qt.KeyboardModifier.MetaModifier
            )
            if (
                event.modifiers() & Qt.KeyboardModifier.AltModifier
                and not (event.modifiers() & blocked_modifiers)
                and self._begin_3d_rotation_drag(
                    scene_pos,
                    clicked_item.atom_id if isinstance(clicked_item, AtomItem) else None,
                )
            ):
                return
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
                additive = bool(
                    event.modifiers()
                    & (Qt.KeyboardModifier.ShiftModifier | Qt.KeyboardModifier.ControlModifier)
                )
                free_select = self.current_tool == "tool_select_lasso" or bool(
                    event.modifiers() & Qt.KeyboardModifier.ControlModifier
                )
                self._begin_selection_drag(scene_pos, free_select, additive)
                return
            if event.modifiers() & (Qt.KeyboardModifier.ShiftModifier | Qt.KeyboardModifier.ControlModifier):
                clicked_item.setSelected(not clicked_item.isSelected())
                self._sync_selection_from_scene()
                return
            else:
                if not clicked_item.isSelected():
                    self.scene.clearSelection()
                    clicked_item.setSelected(True)
            self._sync_selection_from_scene()
            if clicked_item.isSelected():
                if isinstance(clicked_item, TextAnnotationItem) and (
                    clicked_item.textInteractionFlags()
                    & Qt.TextInteractionFlag.TextEditorInteraction
                ):
                    return
                if isinstance(clicked_item, (AtomItem, TextAnnotationItem, ArrowItem, BracketItem)):
                    self._begin_drag(scene_pos)
            return

        clicked_atom_id, clicked_bond_id = self._pick_hover_target(scene_pos)

        if self.current_tool == "tool_rotate_3d_precise":
            self._prompt_precise_3d_rotation(
                clicked_atom_id=clicked_atom_id,
                clicked_bond_id=clicked_bond_id,
            )
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
        if self.current_tool in arrow_tools:
            arrow_kind = arrow_tools[self.current_tool]
            curved_tools = {"tool_arrow_curved", "tool_arrow_curved_fishhook"}
            is_curved_tool = self.current_tool in curved_tools
            if self._arrow_start_pos is None:
                self._arrow_start_pos = scene_pos
                self._arrow_end_pos = None
                self._arrow_curve_adjust_mode = False
                self._arrow_curve_factor = ArrowItem.CURVE_FACTOR_DEFAULT
                return
            if is_curved_tool:
                if self._arrow_end_pos is None:
                    if (scene_pos - self._arrow_start_pos).manhattanLength() < 2.0:
                        self._arrow_start_pos = None
                        self._arrow_end_pos = None
                        self._arrow_curve_adjust_mode = False
                        self._arrow_curve_factor = ArrowItem.CURVE_FACTOR_DEFAULT
                        self._preview_arrow_item.hide_preview()
                        return
                    self._arrow_end_pos = scene_pos
                    self._arrow_curve_adjust_mode = True
                    self._preview_arrow_item.update_preview(
                        self._arrow_start_pos,
                        self._arrow_end_pos,
                        arrow_kind,
                        curve_factor=self._arrow_curve_factor,
                    )
                    return
                cmd = AddArrowCommand(
                    self,
                    self._arrow_start_pos,
                    self._arrow_end_pos,
                    arrow_kind,
                    curve_factor=self._arrow_curve_factor,
                )
                self.undo_stack.push(cmd)
                self._arrow_start_pos = None
                self._arrow_end_pos = None
                self._arrow_curve_adjust_mode = False
                self._arrow_curve_factor = ArrowItem.CURVE_FACTOR_DEFAULT
                self._preview_arrow_item.hide_preview()
                return
            if (scene_pos - self._arrow_start_pos).manhattanLength() < 2.0:
                self._arrow_start_pos = None
                self._arrow_end_pos = None
                self._arrow_curve_adjust_mode = False
                self._arrow_curve_factor = ArrowItem.CURVE_FACTOR_DEFAULT
                self._preview_arrow_item.hide_preview()
                return
            cmd = AddArrowCommand(self, self._arrow_start_pos, scene_pos, arrow_kind)
            self.undo_stack.push(cmd)
            self._arrow_start_pos = None
            self._arrow_end_pos = None
            self._arrow_curve_adjust_mode = False
            self._arrow_curve_factor = ArrowItem.CURVE_FACTOR_DEFAULT
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
                if self.current_tool == "tool_charge":
                    self._insert_symbol_item("±", scene_pos, None, 1.0, False)
                return
            atom = self.model.get_atom(clicked_atom_id)
            if self.current_tool == "tool_charge_plus":
                new_charge = int(atom.charge) + 1
            elif self.current_tool == "tool_charge_minus":
                new_charge = int(atom.charge) - 1
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

        if self.current_tool == "tool_symbol_wavy_anchor":
            anchor_atom_id = clicked_atom_id
            if anchor_atom_id is None:
                item = self._get_item_at(scene_pos)
                if isinstance(item, AtomItem):
                    anchor_atom_id = item.atom_id
            if anchor_atom_id is None:
                return
            self._insert_wavy_anchor(scene_pos, anchor_atom_id)
            return

        symbol_spec = SYMBOL_TEXT_TOOLS.get(self.current_tool)
        if symbol_spec is not None:
            anchor_atom_id = clicked_atom_id
            if anchor_atom_id is None:
                item = self._get_item_at(scene_pos)
                if isinstance(item, AtomItem):
                    anchor_atom_id = item.atom_id
            if anchor_atom_id is None and self.model.atoms:
                atoms = [(atom.id, atom.x, atom.y) for atom in self.model.atoms.values()]
                pick_radius = max(ATOM_HIT_RADIUS * 3, self.state.bond_length * 1.5)
                anchor_atom_id = closest_atom(scene_pos, atoms, pick_radius)
            if symbol_spec.get("electrons"):
                self._insert_electron_dots(
                    scene_pos,
                    anchor_atom_id,
                    symbol_spec["electrons"],
                    symbol_spec.get("scale", 1.0),
                    spread=None,
                    mode="lone_pair" if self.current_tool == "tool_symbol_lone_pair" else "default",
                )
                return
            self._insert_symbol_item(
                symbol_spec["text"],
                scene_pos,
                anchor_atom_id,
                symbol_spec.get("scale", 1.0),
                symbol_spec.get("anchor", True),
                symbol_spec.get("rotate", False),
            )
            return

        if self.current_tool == "tool_atom":
            element = self.state.default_element
            if clicked_atom_id is None:
                is_explicit = None
                if element == "C" and not self.state.show_implicit_carbons:
                    is_explicit = True
                elif element == "H" and not self.state.show_implicit_hydrogens:
                    is_explicit = True
                cmd = AddAtomCommand(
                    self.model,
                    self,
                    element,
                    scene_pos.x(),
                    scene_pos.y(),
                    is_explicit=is_explicit,
                )
                self.undo_stack.push(cmd)
            else:
                cmd = ChangeAtomCommand(self.model, self, clicked_atom_id, element)
                self.undo_stack.push(cmd)
            return

        if self.current_tool == "tool_coordination_center":
            element = self.state.default_element or "Pd"
            if element == "C":
                element = "Pd"
            if clicked_atom_id is None:
                cmd = AddAtomCommand(
                    self.model,
                    self,
                    element,
                    scene_pos.x(),
                    scene_pos.y(),
                    is_explicit=True,
                    is_coordination_center=True,
                    sphere_radius=16.0,
                    sphere_filled=True,
                )
                self.undo_stack.push(cmd)
            else:
                atom = self.model.get_atom(clicked_atom_id)
                needs_center_flag = not getattr(atom, "is_coordination_center", False)
                if not needs_center_flag:
                    return
                self.undo_stack.push(
                    SetCoordinationCenterCommand(
                        self.model,
                        self,
                        clicked_atom_id,
                        True,
                    )
                )
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
            if self.state.active_ring_template:
                self._insert_ring_template(scene_pos)
                return
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
        """Método auxiliar para mouseMoveEvent.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._panning and self._pan_last_pos is not None:
            delta = event.pos() - self._pan_last_pos
            self._pan_last_pos = event.pos()
            hbar = self.horizontalScrollBar()
            vbar = self.verticalScrollBar()
            hbar.setValue(hbar.value() - delta.x())
            vbar.setValue(vbar.value() - delta.y())
            event.accept()
            return
        scene_pos = self.mapToScene(event.pos())
        self._last_scene_pos = scene_pos

        if self._is_rotating_3d:
            self._update_3d_rotation_drag(scene_pos)
            return

        if self._scale_dragging:
            self._update_scale_drag(scene_pos)
            return
        if self._rotation_dragging:
            self._update_rotation_drag(scene_pos)
            return

        if self._select_drag_mode is not None:
            self._update_selection_drag(scene_pos)
            return

        if self._drag_mode == "flex_adjust":
            self._update_flex_bond_preview(scene_pos, event.modifiers())
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
            curved_tools = {"tool_arrow_curved", "tool_arrow_curved_fishhook"}
            is_curved_tool = self.current_tool in curved_tools
            if is_curved_tool and self._arrow_end_pos is not None and self._arrow_curve_adjust_mode:
                curve_factor = self._curve_factor_from_pointer(
                    self._arrow_start_pos,
                    self._arrow_end_pos,
                    scene_pos,
                )
                if event.modifiers() & Qt.KeyboardModifier.ShiftModifier:
                    curve_factor = -curve_factor
                self._arrow_curve_factor = curve_factor
                self._preview_arrow_item.update_preview(
                    self._arrow_start_pos,
                    self._arrow_end_pos,
                    arrow_kind,
                    curve_factor=self._arrow_curve_factor,
                )
            else:
                preview_end = self._arrow_end_pos if (is_curved_tool and self._arrow_end_pos is not None) else scene_pos
                preview_curve = self._arrow_curve_factor if is_curved_tool else None
                self._preview_arrow_item.update_preview(
                    self._arrow_start_pos,
                    preview_end,
                    arrow_kind,
                    curve_factor=preview_curve,
                )
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

                # Move arrows
                if hasattr(self, "_drag_start_arrow_positions"):
                    for item, (start, end) in self._drag_start_arrow_positions.items():
                        item.update_positions(start + delta, end + delta)

                # Move brackets
                if hasattr(self, "_drag_start_bracket_rects"):
                    for item, rect in self._drag_start_bracket_rects.items():
                        item.set_rect(rect.translated(delta))
                
                self._update_selection_overlay()
            return

        self._update_hover(scene_pos)

        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event) -> None:
        """Método auxiliar para mouseReleaseEvent.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._panning:
            self._panning = False
            self._pan_last_pos = None
            self.setCursor(Qt.CursorShape.ArrowCursor)
            event.accept()
            return
        if self._is_rotating_3d:
            self._finalize_3d_rotation_drag()
            return
        if self._scale_dragging:
            self._finalize_scale_drag()
            return
        if self._rotation_dragging:
            self._finalize_rotation_drag()
            return
        if self._select_drag_mode is not None:
            self._finalize_selection_drag()
            return
        if self._drag_mode == "flex_adjust":
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
            kind = self.state.active_bracket_type
            pair = self._split_bracket_kind(kind)
            if pair:
                self.undo_stack.beginMacro("Add brackets")
                for side in pair:
                    self.undo_stack.push(AddBracketCommand(self, bbox, side))
                self.undo_stack.endMacro()
            else:
                cmd = AddBracketCommand(self, bbox, kind)
                self.undo_stack.push(cmd)
            return

        if self._dragging_selection:
            if self._selection_move_handle is not None:
                self._selection_move_handle.setCursor(Qt.CursorShape.OpenHandCursor)
            if self._drag_has_moved:
                # Capture end states
                move_atoms = bool(self._drag_start_positions)
                move_text = bool(getattr(self, "_drag_start_text_positions", {}))
                move_arrows = bool(getattr(self, "_drag_start_arrow_positions", {}))
                move_brackets = bool(getattr(self, "_drag_start_bracket_rects", {}))
                
                if sum([move_atoms, move_text, move_arrows, move_brackets]) > 1:
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

                if move_arrows:
                    before = self._drag_start_arrow_positions
                    after = {
                        item: (item.start_point(), item.end_point())
                        for item in before.keys()
                    }
                    self.undo_stack.push(MoveArrowItemsCommand(self, before, after))

                if move_brackets:
                    before = self._drag_start_bracket_rects
                    after = {item: item.base_rect() for item in before.keys()}
                    self.undo_stack.push(MoveBracketItemsCommand(self, before, after))
                    
                if sum([move_atoms, move_text, move_arrows, move_brackets]) > 1:
                    self.undo_stack.endMacro()

            self._dragging_selection = False
            self._drag_start_pos = None
            self._drag_start_positions = {}
            self._drag_start_text_positions = {}
            self._drag_start_arrow_positions = {}
            self._drag_start_bracket_rects = {}
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
        wavy_items: Iterable = (),
    ) -> None:
        """Método auxiliar para  delete selection.

        Args:
            atom_ids: Descripción del parámetro.
            bond_ids: Descripción del parámetro.
            arrow_items: Descripción del parámetro.
            bracket_items: Descripción del parámetro.
            text_items: Descripción del parámetro.
            wavy_items: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        extra_text_items = list(text_items)
        extra_wavy_items = list(wavy_items)
        if atom_ids and self._electron_dots:
            for item in list(self._electron_dots):
                anchor_id = item.data(ELECTRON_ANCHOR_ROLE)
                if anchor_id in atom_ids:
                    extra_text_items.append(item)
        if atom_ids and self._wavy_anchors:
            for item in list(self._wavy_anchors):
                anchor_id = item.data(WAVY_ANCHOR_ROLE)
                if anchor_id in atom_ids:
                    extra_wavy_items.append(item)
        # Prioritize delete logic
        cmd = DeleteSelectionCommand(
            self.model,
            self,
            atom_ids,
            bond_ids,
            arrow_items=arrow_items,
            bracket_items=bracket_items,
            text_items=extra_text_items,
            wavy_items=extra_wavy_items,
        )
        self.undo_stack.push(cmd)
        self.scene.clearSelection()

    def _get_item_at(self, scene_pos: QPointF):
        """Método auxiliar para  get item at.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        for item in self.scene.items(scene_pos):
             if isinstance(item, (AtomItem, BondItem, ArrowItem, BracketItem, TextAnnotationItem, WavyAnchorItem)):
                return item
             # If we clicked a label/text child of an atom, return the atom.
             if isinstance(item, QGraphicsTextItem):
                parent = item.parentItem()
                if isinstance(parent, AtomItem):
                    return parent
        return None

    def _symbol_insert_position(self, scene_pos: QPointF, atom_id: Optional[int]) -> QPointF:
        """Método auxiliar para  symbol insert position.

        Args:
            scene_pos: Descripción del parámetro.
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if atom_id is None:
            return scene_pos
        atom = self.model.get_atom(atom_id)
        center = QPointF(atom.x, atom.y)
        direction = scene_pos - center
        if direction.manhattanLength() < 2.0:
            direction = self._label_open_direction(atom_id)
        if direction.manhattanLength() < 1e-3:
            direction = QPointF(1.0, -1.0)
        length = math.hypot(direction.x(), direction.y())
        if length <= 1e-6:
            return scene_pos
        nx = direction.x() / length
        ny = direction.y() / length
        offset = max(12.0, self.state.bond_length * 0.3)
        return QPointF(center.x() + nx * offset, center.y() + ny * offset)

    def _symbol_anchor_data(
        self,
        atom_id: int,
        scene_pos: QPointF,
        offset: Optional[float] = None,
    ) -> tuple[QPointF, float]:
        """Método auxiliar para  symbol anchor data.

        Args:
            atom_id: Descripción del parámetro.
            scene_pos: Descripción del parámetro.
            offset: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atom = self.model.get_atom(atom_id)
        center = QPointF(atom.x, atom.y)
        direction = scene_pos - center
        if direction.manhattanLength() < 2.0:
            direction = self._label_open_direction(atom_id)
        if direction.manhattanLength() < 1e-3:
            direction = QPointF(1.0, -1.0)
        length = math.hypot(direction.x(), direction.y())
        if length <= 1e-6:
            return center, 0.0
        nx = direction.x() / length
        ny = direction.y() / length
        use_offset = max(12.0, self.state.bond_length * 0.3) if offset is None else offset
        pos = QPointF(center.x() + nx * use_offset, center.y() + ny * use_offset)
        angle = math.degrees(math.atan2(ny, nx))
        return pos, angle

    def _electron_dot_font_size(self, atom_id: int, scale: float) -> float:
        """Método auxiliar para  electron dot font size.

        Args:
            atom_id: Descripción del parámetro.
            scale: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = self.atom_items.get(atom_id)
        if item is not None:
            base_font = item.label.font()
        else:
            base_font = self._label_font()
        base_size = base_font.pointSizeF() or float(base_font.pointSize())
        return max(5.0, base_size * 0.55 * scale)

    def _electron_radial_offset(self, atom_id: int, scale: float) -> float:
        """Método auxiliar para  electron radial offset.

        Args:
            atom_id: Descripción del parámetro.
            scale: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = self.atom_items.get(atom_id)
        dot_size = self._electron_dot_font_size(atom_id, scale)
        if item is not None and item.label.isVisible():
            rect = item.label.boundingRect()
            label_radius = max(rect.width(), rect.height()) * 0.5
            return label_radius + dot_size * 0.45 + 2.0
        return max(14.0, self.state.bond_length * 0.25 + dot_size * 0.4)

    def _electron_pair_spread(self, atom_id: int, scale: float) -> float:
        """Método auxiliar para  electron pair spread.

        Args:
            atom_id: Descripción del parámetro.
            scale: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        dot_size = self._electron_dot_font_size(atom_id, scale)
        return max(3.0, dot_size * 0.55)

    def _wavy_anchor_length(self) -> float:
        """Método auxiliar para  wavy anchor length.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        return max(18.0, self.state.bond_length * 1.5)

    def _wavy_anchor_bond_angle(self, atom_id: int, scene_pos: QPointF) -> tuple[Optional[int], float]:
        """Método auxiliar para  wavy anchor bond angle.

        Args:
            atom_id: Descripción del parámetro.
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atom = self.model.get_atom(atom_id)
        center = QPointF(atom.x, atom.y)
        direction = scene_pos - center
        if direction.manhattanLength() < 2.0:
            direction = self._label_open_direction(atom_id)
        target_angle = None
        if direction.manhattanLength() >= 1e-3:
            target_angle = (
                math.degrees(math.atan2(direction.y(), direction.x())) + 360.0
            ) % 360.0

        bonds: list[tuple[int, float]] = []
        for bond in self.model.bonds.values():
            if bond.a1_id == atom_id:
                other = self.model.get_atom(bond.a2_id)
            elif bond.a2_id == atom_id:
                other = self.model.get_atom(bond.a1_id)
            else:
                continue
            dx = other.x - atom.x
            dy = other.y - atom.y
            angle = (math.degrees(math.atan2(dy, dx)) + 360.0) % 360.0
            bonds.append((bond.id, angle))

        if bonds:
            if target_angle is None:
                return bonds[0][0], bonds[0][1]
            best = min(bonds, key=lambda item: angle_distance_deg(item[1], target_angle))
            return best[0], best[1]

        if target_angle is None:
            return None, 0.0
        return None, target_angle

    def _position_wavy_anchor(self, item: WavyAnchorItem) -> None:
        """Método auxiliar para  position wavy anchor.

        Args:
            item: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        anchor_id = item.data(WAVY_ANCHOR_ROLE)
        base_angle = item.data(WAVY_ANCHOR_ANGLE_ROLE)
        bond_id = item.data(WAVY_ANCHOR_BOND_ROLE)
        length = item.data(WAVY_ANCHOR_LENGTH_ROLE)
        if anchor_id is None:
            return
        if anchor_id not in self.model.atoms:
            return
        atom = self.model.get_atom(anchor_id)
        use_length = float(length) if length is not None else self._wavy_anchor_length()

        angle_deg_value: Optional[float] = None
        if bond_id is not None and bond_id in self.model.bonds:
            bond = self.model.get_bond(int(bond_id))
            if bond.a1_id == anchor_id:
                other = self.model.get_atom(bond.a2_id)
            elif bond.a2_id == anchor_id:
                other = self.model.get_atom(bond.a1_id)
            else:
                other = None
            if other is not None:
                dx = other.x - atom.x
                dy = other.y - atom.y
                angle_deg_value = (math.degrees(math.atan2(dy, dx)) + 360.0) % 360.0

        if angle_deg_value is None and base_angle is not None:
            angle_deg_value = float(base_angle)
        if angle_deg_value is None:
            return

        perp_angle = math.radians(angle_deg_value + 90.0)
        nx = math.cos(perp_angle)
        ny = math.sin(perp_angle)
        center = QPointF(atom.x, atom.y)
        half = use_length * 0.5
        start = QPointF(center.x() - nx * half, center.y() - ny * half)
        end = QPointF(center.x() + nx * half, center.y() + ny * half)
        try:
            item.update_positions(start, end)
        except RuntimeError:
            self._wavy_anchors.discard(item)

    def _update_wavy_anchors_for_atom(self, atom_id: int) -> None:
        """Método auxiliar para  update wavy anchors for atom.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        for item in list(self._wavy_anchors):
            if item.data(WAVY_ANCHOR_ROLE) != atom_id:
                continue
            self._position_wavy_anchor(item)

    def _bond_angles_for_atom(self, atom_id: int) -> list[float]:
        """Método auxiliar para  bond angles for atom.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atom = self.model.get_atom(atom_id)
        angles: list[float] = []
        for bond in self.model.bonds.values():
            if bond.a1_id == atom_id:
                other = self.model.get_atom(bond.a2_id)
            elif bond.a2_id == atom_id:
                other = self.model.get_atom(bond.a1_id)
            else:
                continue
            dx = other.x - atom.x
            dy = other.y - atom.y
            angle = (math.degrees(math.atan2(dy, dx)) + 360.0) % 360.0
            angles.append(angle)
        return angles

    def _occupied_electron_slots(self, atom_id: int) -> set[int]:
        """Método auxiliar para  occupied electron slots.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        occupied: set[int] = set()
        for item in self._electron_dots:
            if item.data(ELECTRON_ANCHOR_ROLE) != atom_id:
                continue
            slot_idx = item.data(ELECTRON_SLOT_ROLE)
            if slot_idx is None:
                continue
            occupied.add(int(slot_idx))
        return occupied

    def _candidate_electron_slots(self, atom_id: int, mode: str) -> list[int]:
        """Método auxiliar para  candidate electron slots.

        Args:
            atom_id: Descripción del parámetro.
            mode: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        bond_angles = self._bond_angles_for_atom(atom_id)
        bond_count = len(bond_angles)
        if mode != "lone_pair":
            return list(range(len(ELECTRON_SLOT_ANGLES)))

        # Lone pair tool: 90-degree spacing when free, tetrahedral when bonded.
        if bond_count == 0:
            return [0, 2, 4, 6]

        ideal_angles: list[float] = []
        if bond_count == 1:
            b = bond_angles[0]
            ideal_angles = [b + 109.5, b - 109.5, b + 180.0]
        elif bond_count == 2:
            b1, b2 = bond_angles[0], bond_angles[1]
            ideal_angles = [
                b1 + 109.5,
                b1 - 109.5,
                b2 + 109.5,
                b2 - 109.5,
            ]
        else:
            ideal_angles = []

        candidates: list[int] = []
        if ideal_angles:
            for ang in ideal_angles:
                ang = (ang + 360.0) % 360.0
                best_idx = None
                best_dist = 1e9
                for idx, slot_angle in enumerate(ELECTRON_SLOT_ANGLES):
                    dist = angle_distance_deg(slot_angle, ang)
                    if dist < best_dist:
                        best_dist = dist
                        best_idx = idx
                if best_idx is not None and best_idx not in candidates:
                    candidates.append(best_idx)

        if not candidates:
            candidates = list(range(len(ELECTRON_SLOT_ANGLES)))

        # Remove slots blocked by bonds
        filtered: list[int] = []
        for idx in candidates:
            slot_angle = ELECTRON_SLOT_ANGLES[idx]
            blocked = any(
                angle_distance_deg(slot_angle, b) <= ELECTRON_SLOT_TOLERANCE_DEG
                for b in bond_angles
            )
            if not blocked:
                filtered.append(idx)
        return filtered if filtered else candidates

    def _select_electron_slot(
        self,
        atom_id: int,
        scene_pos: QPointF,
        candidate_slots: Optional[list[int]] = None,
    ) -> Optional[int]:
        """Método auxiliar para  select electron slot.

        Args:
            atom_id: Descripción del parámetro.
            scene_pos: Descripción del parámetro.
            candidate_slots: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atom = self.model.get_atom(atom_id)
        center = QPointF(atom.x, atom.y)
        direction = scene_pos - center
        if direction.manhattanLength() < 2.0:
            direction = self._label_open_direction(atom_id)
        if direction.manhattanLength() < 1e-3:
            direction = QPointF(1.0, 0.0)
        target_angle = (math.degrees(math.atan2(direction.y(), direction.x())) + 360.0) % 360.0

        bond_angles = self._bond_angles_for_atom(atom_id)
        occupied = self._occupied_electron_slots(atom_id)

        if candidate_slots is None:
            candidate_slots = list(range(len(ELECTRON_SLOT_ANGLES)))

        best_idx: Optional[int] = None
        best_dist = 1e9
        for idx in candidate_slots:
            slot_angle = ELECTRON_SLOT_ANGLES[idx]
            if idx in occupied:
                continue
            blocked = any(angle_distance_deg(slot_angle, b) <= ELECTRON_SLOT_TOLERANCE_DEG for b in bond_angles)
            if blocked:
                continue
            dist = angle_distance_deg(slot_angle, target_angle)
            if dist < best_dist:
                best_dist = dist
                best_idx = idx
        return best_idx

    def _electron_slot_position(self, atom_id: int, slot_idx: int, scale: float) -> tuple[QPointF, QPointF]:
        """Método auxiliar para  electron slot position.

        Args:
            atom_id: Descripción del parámetro.
            slot_idx: Descripción del parámetro.
            scale: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atom = self.model.get_atom(atom_id)
        angle = math.radians(ELECTRON_SLOT_ANGLES[slot_idx])
        nx = math.cos(angle)
        ny = math.sin(angle)
        px, py = -ny, nx
        radial_offset = self._electron_radial_offset(atom_id, scale)
        center = QPointF(atom.x, atom.y)
        pos = QPointF(center.x() + nx * radial_offset, center.y() + ny * radial_offset)
        tangent = QPointF(px, py)
        return pos, tangent

    def _position_electron_dot(self, item: TextAnnotationItem) -> None:
        """Método auxiliar para  position electron dot.

        Args:
            item: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if not item.data(ELECTRON_DOT_ROLE):
            return
        anchor_id = item.data(ELECTRON_ANCHOR_ROLE)
        slot_idx = item.data(ELECTRON_SLOT_ROLE)
        side = item.data(ELECTRON_SIDE_ROLE) or 0
        if anchor_id is None or slot_idx is None:
            return
        scale = float(item.data(ELECTRON_SCALE_ROLE) or 1.0)
        pos, tangent = self._electron_slot_position(int(anchor_id), int(slot_idx), scale)
        offset = self._electron_pair_spread(int(anchor_id), scale) * int(side)
        pos = QPointF(pos.x() + tangent.x() * offset, pos.y() + tangent.y() * offset)
        # Update font size to match atom label scale
        item_font = item.font()
        item_font.setPointSizeF(self._electron_dot_font_size(int(anchor_id), scale))
        item_font.setBold(False)
        item.setFont(item_font)
        rect = item.boundingRect()
        item.setPos(pos.x() - rect.width() / 2, pos.y() - rect.height() / 2)

    def _update_electron_dots_for_atom(self, atom_id: int) -> None:
        """Método auxiliar para  update electron dots for atom.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        for item in list(self._electron_dots):
            if item.data(ELECTRON_ANCHOR_ROLE) != atom_id:
                continue
            self._position_electron_dot(item)

    def _create_dot_item(
        self,
        pos: QPointF,
        atom_item: Optional[AtomItem],
        scale: float,
        anchor_atom_id: Optional[int] = None,
        slot_idx: Optional[int] = None,
        side: int = 0,
    ) -> None:
        """Método auxiliar para  create dot item.

        Args:
            pos: Descripción del parámetro.
            atom_item: Descripción del parámetro.
            scale: Descripción del parámetro.
            anchor_atom_id: Descripción del parámetro.
            slot_idx: Descripción del parámetro.
            side: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = TextAnnotationItem("•", 0.0, 0.0)
        self._apply_text_settings(item)
        item.setDefaultTextColor(QColor("#222222"))
        item.setTextInteractionFlags(Qt.TextInteractionFlag.NoTextInteraction)
        try:
            item.document().setDocumentMargin(0)
            item.document().setDefaultStyleSheet("body { background: transparent; }")
        except Exception:
            pass
        cursor = item.textCursor()
        cursor.clearSelection()
        item.setTextCursor(cursor)
        fmt = QTextCharFormat()
        fmt.setBackground(QBrush(Qt.BrushStyle.NoBrush))
        cursor.select(cursor.SelectionType.Document)
        cursor.mergeCharFormat(fmt)

        if scale and abs(scale - 1.0) > 1e-3:
            font = item.font()
            font.setPointSizeF(font.pointSizeF() * scale)
            item.setFont(font)
        if anchor_atom_id is not None:
            font = item.font()
            font.setPointSizeF(self._electron_dot_font_size(anchor_atom_id, scale))
            font.setBold(False)
            item.setFont(font)

        rect = item.boundingRect()
        item.setPos(pos.x() - rect.width() / 2, pos.y() - rect.height() / 2)
        item.setData(ELECTRON_DOT_ROLE, True)
        if anchor_atom_id is not None:
            item.setData(ELECTRON_ANCHOR_ROLE, anchor_atom_id)
        if slot_idx is not None:
            item.setData(ELECTRON_SLOT_ROLE, slot_idx)
            item.setData(ELECTRON_SIDE_ROLE, int(side))
        item.setData(ELECTRON_SCALE_ROLE, float(scale))
        self._electron_dots.add(item)

        self.undo_stack.push(AddTextItemCommand(self, item))

    def _insert_electron_dots(
        self,
        scene_pos: QPointF,
        atom_id: Optional[int],
        count: int,
        scale: float,
        spread: Optional[float] = None,
        mode: str = "default",
    ) -> None:
        """Método auxiliar para  insert electron dots.

        Args:
            scene_pos: Descripción del parámetro.
            atom_id: Descripción del parámetro.
            count: Descripción del parámetro.
            scale: Descripción del parámetro.
            spread: Descripción del parámetro.
            mode: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if atom_id is None:
            return
        atom = self.model.get_atom(atom_id)
        if atom.element not in HETERO_ELECTRON_ATOMS:
            return

        slot_candidates = self._candidate_electron_slots(atom_id, mode)
        slot_idx = self._select_electron_slot(atom_id, scene_pos, slot_candidates)
        if slot_idx is None:
            return

        pos, tangent = self._electron_slot_position(atom_id, slot_idx, scale)
        if count <= 1:
            self._create_dot_item(
                pos,
                self.atom_items.get(atom_id),
                scale,
                anchor_atom_id=atom_id,
                slot_idx=slot_idx,
                side=0,
            )
            return

        base_spread = self._electron_pair_spread(atom_id, scale)
        if mode == "lone_pair":
            base_spread *= 1.25
        use_spread = base_spread if spread is None else spread
        dx = tangent.x() * use_spread
        dy = tangent.y() * use_spread
        p1 = QPointF(pos.x() + dx, pos.y() + dy)
        p2 = QPointF(pos.x() - dx, pos.y() - dy)
        self._create_dot_item(
            p1,
            self.atom_items.get(atom_id),
            scale,
            anchor_atom_id=atom_id,
            slot_idx=slot_idx,
            side=1,
        )
        self._create_dot_item(
            p2,
            self.atom_items.get(atom_id),
            scale,
            anchor_atom_id=atom_id,
            slot_idx=slot_idx,
            side=-1,
        )

    def _insert_wavy_anchor(self, scene_pos: QPointF, atom_id: int) -> None:
        """Método auxiliar para  insert wavy anchor.

        Args:
            scene_pos: Descripción del parámetro.
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if atom_id not in self.model.atoms:
            return
        bond_id, base_angle = self._wavy_anchor_bond_angle(atom_id, scene_pos)
        length = self._wavy_anchor_length()

        atom = self.model.get_atom(atom_id)
        perp_angle = math.radians(base_angle + 90.0)
        nx = math.cos(perp_angle)
        ny = math.sin(perp_angle)
        center = QPointF(atom.x, atom.y)
        half = length * 0.5
        start = QPointF(center.x() - nx * half, center.y() - ny * half)
        end = QPointF(center.x() + nx * half, center.y() + ny * half)

        item = WavyAnchorItem(start, end, style=self.drawing_style)
        item.setData(WAVY_ANCHOR_ROLE, atom_id)
        item.setData(WAVY_ANCHOR_ANGLE_ROLE, base_angle)
        item.setData(WAVY_ANCHOR_LENGTH_ROLE, length)
        if bond_id is not None:
            item.setData(WAVY_ANCHOR_BOND_ROLE, int(bond_id))
        self._wavy_anchors.add(item)
        self.undo_stack.push(AddWavyAnchorCommand(self, item))

    def _insert_symbol_item(
        self,
        text: str,
        scene_pos: QPointF,
        atom_id: Optional[int],
        scale: float,
        anchor_to_atom: bool,
        rotate: bool = False,
    ) -> None:
        """Método auxiliar para  insert symbol item.

        Args:
            text: Descripción del parámetro.
            scene_pos: Descripción del parámetro.
            atom_id: Descripción del parámetro.
            scale: Descripción del parámetro.
            anchor_to_atom: Descripción del parámetro.
            rotate: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        pos = scene_pos
        angle = 0.0
        if anchor_to_atom and atom_id is not None:
            pos, angle = self._symbol_anchor_data(atom_id, scene_pos)
        item = TextAnnotationItem(text, 0.0, 0.0)
        self._apply_text_settings(item)
        item.setDefaultTextColor(QColor("#222222"))
        item.setTextInteractionFlags(Qt.TextInteractionFlag.NoTextInteraction)
        try:
            item.document().setDocumentMargin(0)
            item.document().setDefaultStyleSheet("body { background: transparent; }")
        except Exception:
            pass
        cursor = item.textCursor()
        cursor.clearSelection()
        item.setTextCursor(cursor)
        fmt = QTextCharFormat()
        fmt.setBackground(QBrush(Qt.BrushStyle.NoBrush))
        cursor.select(cursor.SelectionType.Document)
        cursor.mergeCharFormat(fmt)
        if scale and abs(scale - 1.0) > 1e-3:
            font = item.font()
            font.setPointSizeF(font.pointSizeF() * scale)
            item.setFont(font)
        rect = item.boundingRect()
        if anchor_to_atom and atom_id is not None and atom_id in self.atom_items:
            atom_item = self.atom_items[atom_id]
            local = atom_item.mapFromScene(QPointF(pos.x(), pos.y()))
            item.setParentItem(atom_item)
            item.setPos(local.x() - rect.width() / 2, local.y() - rect.height() / 2)
        else:
            item.setPos(pos.x() - rect.width() / 2, pos.y() - rect.height() / 2)
        if rotate:
            item.setRotation(angle)
        self.undo_stack.push(AddTextItemCommand(self, item))

    def mouseDoubleClickEvent(self, event) -> None:
        """Método auxiliar para mouseDoubleClickEvent.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para wheelEvent.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if event.modifiers() & Qt.KeyboardModifier.ShiftModifier:
            hbar = self.horizontalScrollBar()
            delta = event.angleDelta().y()
            if delta == 0:
                delta = event.angleDelta().x()
            if delta != 0:
                hbar.setValue(hbar.value() - delta)
                event.accept()
                return
        if event.angleDelta().y() > 0:
            self.zoom_in()
        else:
            self.zoom_out()

    def keyPressEvent(self, event) -> None:
        """Método auxiliar para keyPressEvent.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        focus_item = self.scene.focusItem()
        if isinstance(focus_item, TextAnnotationItem) and (
            focus_item.textInteractionFlags()
            & Qt.TextInteractionFlag.TextEditorInteraction
        ):
            # Let the text editor handle Delete/Backspace while editing text.
            super().keyPressEvent(event)
            return
        if event.key() == Qt.Key.Key_Space and not self._space_panning:
            self._space_panning = True
            self.setCursor(Qt.CursorShape.OpenHandCursor)
            event.accept()
            return
        if event.key() in (Qt.Key.Key_Delete, Qt.Key.Key_Backspace):
            has_selected_arrows = any(
                isinstance(item, ArrowItem) for item in self.scene.selectedItems()
            )
            has_selected_brackets = any(
                isinstance(item, BracketItem) for item in self.scene.selectedItems()
            )
            has_selected_text = any(
                isinstance(item, TextAnnotationItem) for item in self.scene.selectedItems()
            )
            if (
                self.state.selected_atoms
                or self.state.selected_bonds
                or has_selected_arrows
                or has_selected_brackets
                or has_selected_text
            ):
                self.delete_selection()
                return
            if self._delete_hovered():
                return
        if event.key() == Qt.Key.Key_Escape and self._pending_template_graph is not None:
            self.cancel_template_insert_mode()
            event.accept()
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

    def remember_text_edit_item(self, item: TextAnnotationItem | None) -> None:
        """Método auxiliar para remember text edit item.

        Args:
            item: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._last_text_edit_item = item

    def restore_text_edit_focus(self) -> None:
        """Método auxiliar para restore text edit focus.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = getattr(self, "_last_text_edit_item", None)
        if item is None:
            return
        if item.scene() is not self.scene:
            self._last_text_edit_item = None
            return
        item.setTextInteractionFlags(Qt.TextInteractionFlag.TextEditorInteraction)
        item.setFocus()

    def keyReleaseEvent(self, event) -> None:
        """Método auxiliar para keyReleaseEvent.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if event.key() == Qt.Key.Key_Space and self._space_panning:
            self._space_panning = False
            if not self._panning:
                self.setCursor(Qt.CursorShape.ArrowCursor)
            event.accept()
            return
        super().keyReleaseEvent(event)

    def _handle_atom_text_entry(self, event) -> bool:
        """Método auxiliar para  handle atom text entry.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  handle nudge.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        selected_arrows = self._selected_arrow_items()
        selected_brackets = self._selected_bracket_items()
        if not selected_atom_ids and not selected_text_items and not selected_arrows and not selected_brackets:
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

        if selected_arrows:
            before = {item: (item.start_point(), item.end_point()) for item in selected_arrows}
            after = {
                item: (start + QPointF(dx, dy), end + QPointF(dx, dy))
                for item, (start, end) in before.items()
            }
            self.undo_stack.push(MoveArrowItemsCommand(self, before, after))

        if selected_brackets:
            before = {item: item.base_rect() for item in selected_brackets}
            after = {item: rect.translated(dx, dy) for item, rect in before.items()}
            self.undo_stack.push(MoveBracketItemsCommand(self, before, after))

        self._update_selection_overlay()
        return True

    def _select_all_items(self) -> None:
        """Método auxiliar para  select all items.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self.scene.clearSelection()
        for item in self.scene.items():
            if isinstance(
                item,
                (AtomItem, BondItem, ArrowItem, BracketItem, TextAnnotationItem, WavyAnchorItem),
            ) and item.isVisible():
                item.setSelected(True)
        self._sync_selection_from_scene()

    def _handle_select_all(self, event) -> bool:
        """Método auxiliar para  handle select all.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  normalize atom label.

        Args:
            text: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  atom label items.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        elements = ["C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "H"]
        return elements + FUNCTIONAL_GROUP_LABELS

    def _label_anchor_candidates(self, label: str) -> list[str]:
        """Método auxiliar para  label anchor candidates.

        Args:
            label: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  prompt anchor for atom.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  cycle anchor override.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  prompt atom label.

        Args:
            current: Descripción del parámetro.
            atom_id: Descripción del parámetro.
            initial: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
                "text_items": [],
                "wavy_anchors": []
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
                    "end": {"x": item.end_point().x(), "y": item.end_point().y()},
                    "curve_factor": item.curve_factor(),
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
                if bool(item.data(NUMBERING_TEXT_ROLE)):
                    continue
                data["annotations"]["text_items"].append({
                    "text": item.toPlainText(),
                    "html": item.toHtml(), # Store HTML for rich text (formatting, colors)
                    "x": item.pos().x(),
                    "y": item.pos().y(),
                    "rotation": item.rotation(),
                    "font": item.font().toString(),
                    "color": item.defaultTextColor().name()
                })
            elif isinstance(item, WavyAnchorItem):
                anchor_id = item.data(WAVY_ANCHOR_ROLE)
                angle = item.data(WAVY_ANCHOR_ANGLE_ROLE)
                length = item.data(WAVY_ANCHOR_LENGTH_ROLE)
                bond_id = item.data(WAVY_ANCHOR_BOND_ROLE)
                if anchor_id is None:
                    continue
                data["annotations"]["wavy_anchors"].append({
                    "anchor_id": int(anchor_id),
                    "angle": float(angle) if angle is not None else 0.0,
                    "length": float(length) if length is not None else self._wavy_anchor_length(),
                    "bond_id": int(bond_id) if bond_id is not None else None,
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
            curve_factor = arrow_d.get("curve_factor")
            try:
                curve_factor_value = float(curve_factor) if curve_factor is not None else None
            except Exception:
                curve_factor_value = None
            self.add_arrow_item(start, end, kind=arrow_d["kind"], curve_factor=curve_factor_value)

        for br_d in annotations.get("brackets", []):
            rect_d = br_d["rect"]
            rect = QRectF(rect_d["x"], rect_d["y"], rect_d["w"], rect_d["h"])
            kind = br_d.get("kind", "[]")
            padding = br_d.get("padding", 8.0)
            pair = self._split_bracket_kind(kind)
            if pair:
                for side in pair:
                    bracket = BracketItem(rect, kind=side, padding=padding, style=self.drawing_style)
                    self.readd_bracket_item(bracket, rect, side, padding=padding)
            else:
                bracket = BracketItem(rect, kind=kind, padding=padding, style=self.drawing_style)
                self.readd_bracket_item(bracket, rect, kind, padding=padding)

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

        for anchor_d in annotations.get("wavy_anchors", []):
            anchor_id = anchor_d.get("anchor_id")
            if anchor_id is None or anchor_id not in self.model.atoms:
                continue
            angle = float(anchor_d.get("angle", 0.0))
            length = float(anchor_d.get("length", self._wavy_anchor_length()))
            bond_id = anchor_d.get("bond_id")
            item = WavyAnchorItem(QPointF(0.0, 0.0), QPointF(1.0, 0.0), style=self.drawing_style)
            item.setData(WAVY_ANCHOR_ROLE, int(anchor_id))
            item.setData(WAVY_ANCHOR_ANGLE_ROLE, angle)
            item.setData(WAVY_ANCHOR_LENGTH_ROLE, length)
            if bond_id is not None:
                item.setData(WAVY_ANCHOR_BOND_ROLE, int(bond_id))
            self.readd_wavy_anchor_item(item)

        # Full refresh to update atom visibility and circles
        self.refresh_atom_visibility()
        self.refresh_aromatic_circles()

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

    def clear_canvas(self) -> None:
        """Método auxiliar para clear canvas.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        self._electron_dots.clear()
        self._wavy_anchors.clear()
        self._numbering_overlay_items.clear()
        self._numbering_atom_items.clear()
        self._numbering_structure_items.clear()
        self._numbering_atom_base_pos.clear()
        self._numbering_structure_base_pos.clear()
        self._atom_numbering.clear()
        self._structure_numbering.clear()

        self.scene.clear()
        self.model.clear()
        self.undo_stack.clear()
        self.undo_stack.blockSignals(False)
        self.atom_items.clear()
        self.bond_items.clear()
        self.aromatic_circles.clear()
        self.arrow_items.clear()
        self.bracket_items.clear()
        self._electron_dots.clear()
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
        self.recompute_numbering()

    def _rebuild_items_from_model(self) -> None:
        """Create AtomItems and BondItems for everything currently in self.model."""
        self.atom_items.clear()
        self.bond_items.clear()
        self._suspend_numbering_refresh = True
        try:
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
            ring_pairs = self._compute_ring_bond_pairs()
            for bond in self.model.bonds.values():
                a1 = self.model.atoms.get(bond.a1_id)
                a2 = self.model.atoms.get(bond.a2_id)
                if a1 and a2:
                    item = BondItem(
                        bond, a1, a2,
                        render_aromatic_as_circle=self.state.use_aromatic_circles,
                        style=self.drawing_style
                    )
                    ring_center = self._aromatic_ring_center_for_bond(bond) if bond.is_aromatic else None
                    if ring_center is None and bond.ring_id is not None:
                        ring_center = self._ring_centers.get(bond.ring_id)
                    item.set_ring_context(ring_center)
                    item.set_bond_in_ring(self._bond_in_ring_for_pairs(bond, ring_pairs))
                    item.set_offset_sign(self._bond_offset_sign(bond))
                    self._configure_bond_rendering(bond, item)
                    self._set_bond_item_join_context(bond, item)
                    trim_start = self._bond_endpoint_trim(bond, bond.a1_id)
                    trim_end = self._bond_endpoint_trim(bond, bond.a2_id)
                    item.set_endpoint_trim(trim_start, trim_end)
                    extend_start = self._bond_endpoint_extend(bond, bond.a1_id)
                    extend_end = self._bond_endpoint_extend(bond, bond.a2_id)
                    item.set_endpoint_extend(extend_start, extend_end)
                    # Explicitly update positions now that context and sign are set
                    item.update_positions(a1, a2)
                    self.scene.addItem(item)
                    self.bond_items[bond.id] = item

            # 4. Global refresh: labels, visibility, shrinks, and implicit H overlays
            self.refresh_atom_visibility()
            self.refresh_aromatic_circles()
            self.validate_structure()
        finally:
            self._suspend_numbering_refresh = False
        self.recompute_numbering()

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

    def _compute_ring_bond_pairs(self) -> set[frozenset[int]]:
        """Método auxiliar para  compute ring bond pairs.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        view = MolView(self.model)
        rings = find_rings_simple(view)
        if not rings:
            return set()
        pairs: set[frozenset[int]] = set()
        for ring in rings:
            pairs.update(ring_bonds(view, ring))
        return pairs

    @staticmethod
    def _bond_in_ring_for_pairs(bond: Bond, ring_pairs: set[frozenset[int]]) -> bool:
        """Método auxiliar para  bond in ring for pairs.

        Args:
            bond: Descripción del parámetro.
            ring_pairs: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if bond.ring_id is not None:
            return True
        return frozenset({bond.a1_id, bond.a2_id}) in ring_pairs

    def _refresh_bond_ring_flags(
        self,
        ring_pairs: Optional[set[frozenset[int]]] = None,
    ) -> None:
        """Método auxiliar para  refresh bond ring flags.

        Args:
            ring_pairs: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if ring_pairs is None:
            ring_pairs = self._compute_ring_bond_pairs()
        for bond in self.model.bonds.values():
            item = self.bond_items.get(bond.id)
            if item is None:
                continue
            item.set_bond_in_ring(self._bond_in_ring_for_pairs(bond, ring_pairs))
            if item.style in (BondStyle.WEDGE, BondStyle.HASHED):
                atom1 = self.model.get_atom(bond.a1_id)
                atom2 = self.model.get_atom(bond.a2_id)
                if atom1 and atom2:
                    item.update_positions(atom1, atom2)

    def add_text_item(self, item: TextAnnotationItem) -> None:
        """Añade text item.

        Args:
            item: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if item.scene() is not self.scene:
            self.scene.addItem(item)

    def add_wavy_anchor_item(self, item: WavyAnchorItem) -> None:
        """Añade wavy anchor item.

        Args:
            item: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if item.scene() is not self.scene:
            self.scene.addItem(item)
        self._wavy_anchors.add(item)
        self._position_wavy_anchor(item)

    def remove_text_item(self, item: TextAnnotationItem) -> None:
        """Elimina text item.

        Args:
            item: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if item.scene() is self.scene:
            self.scene.removeItem(item)
        if item in self._electron_dots:
            self._electron_dots.discard(item)

    def remove_wavy_anchor_item(self, item: WavyAnchorItem) -> None:
        """Elimina wavy anchor item.

        Args:
            item: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if item.scene() is self.scene:
            self.scene.removeItem(item)
        if item in self._wavy_anchors:
            self._wavy_anchors.discard(item)

    def readd_text_item(self, item: TextAnnotationItem) -> None:
        """Método auxiliar para readd text item.

        Args:
            item: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if item.scene() is not self.scene:
            self.scene.addItem(item)
        if item.data(ELECTRON_DOT_ROLE):
            self._electron_dots.add(item)
            self._position_electron_dot(item)

    def readd_wavy_anchor_item(self, item: WavyAnchorItem) -> None:
        """Método auxiliar para readd wavy anchor item.

        Args:
            item: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if item.scene() is not self.scene:
            self.scene.addItem(item)
        self._wavy_anchors.add(item)
        self._position_wavy_anchor(item)

    def delete_selection(self) -> None:
        """Método auxiliar para delete selection.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        selected_arrows = [
            item for item in self.scene.selectedItems() if isinstance(item, ArrowItem)
        ]
        selected_brackets = [
            item for item in self.scene.selectedItems() if isinstance(item, BracketItem)
        ]
        selected_text_items = [
            item for item in self.scene.selectedItems() if isinstance(item, TextAnnotationItem)
        ]
        selected_wavy = [
            item for item in self.scene.selectedItems() if isinstance(item, WavyAnchorItem)
        ]
        if (
            not self.state.selected_atoms
            and not self.state.selected_bonds
            and not selected_arrows
            and not selected_brackets
            and not selected_text_items
            and not selected_wavy
        ):
            return
        self._delete_selection(
            set(self.state.selected_atoms),
            set(self.state.selected_bonds),
            arrow_items=selected_arrows,
            bracket_items=selected_brackets,
            text_items=selected_text_items,
            wavy_items=selected_wavy,
        )

    def _delete_hovered(self) -> bool:
        """Método auxiliar para  delete hovered.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self.hovered_atom_id is not None:
            self._delete_selection({self.hovered_atom_id}, set())
            return True
        if self.hovered_bond_id is not None:
            self._delete_selection(set(), {self.hovered_bond_id})
            return True
        return False

    def _selected_structure_ids(self) -> tuple[set[int], list[Bond]]:
        """Método auxiliar para  selected structure ids.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atom_ids = self._selected_atom_ids_for_transform()
        if not atom_ids:
            return set(), []
        bonds: list[Bond] = []
        for bond in self.model.bonds.values():
            if bond.a1_id in atom_ids and bond.a2_id in atom_ids:
                bonds.append(bond)
        return atom_ids, bonds

    def _build_selection_graph(self, atom_ids: set[int], bonds: list[Bond]) -> MolGraph:
        """Método auxiliar para  build selection graph.

        Args:
            atom_ids: Descripción del parámetro.
            bonds: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        graph = MolGraph()
        for atom_id in atom_ids:
            atom = self.model.get_atom(atom_id)
            graph.add_atom(
                atom.element,
                atom.x,
                atom.y,
                atom_id=atom.id,
                charge=atom.charge,
                isotope=atom.isotope,
                explicit_h=atom.explicit_h,
                mapping=atom.mapping,
                is_query=atom.is_query,
                is_explicit=atom.is_explicit,
                no_implicit=bool(getattr(atom, "no_implicit", False)),
                is_coordination_center=getattr(atom, "is_coordination_center", False),
                sphere_radius=getattr(atom, "sphere_radius", None),
                sphere_color=getattr(atom, "sphere_color", None),
                sphere_filled=bool(getattr(atom, "sphere_filled", True)),
                sphere_transparent=bool(getattr(atom, "sphere_transparent", False)),
            )
        for bond in bonds:
            graph.add_bond(
                bond.a1_id,
                bond.a2_id,
                bond.order,
                bond_id=bond.id,
                style=bond.style,
                stereo=bond.stereo,
                is_aromatic=bond.is_aromatic,
                display_order=bond.display_order,
                is_query=bond.is_query,
                ring_id=bond.ring_id,
                length_px=bond.length_px,
                stroke_px=bond.stroke_px,
                color=bond.color,
                donor_atom_id=getattr(bond, "donor_atom_id", None),
            )
        return graph

    def _build_selection_payload(self) -> Optional[dict]:
        """Método auxiliar para  build selection payload.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atom_ids, bonds = self._selected_structure_ids()
        arrows = self._selected_arrow_items()
        brackets = self._selected_bracket_items()
        texts = self._selected_text_items()

        if not atom_ids and not bonds and not arrows and not brackets and not texts:
            return None

        bbox = self._selected_items_bbox()
        if bbox is None:
            return None
        left = bbox.left()
        top = bbox.top()

        atoms_payload = []
        for atom_id in atom_ids:
            atom = self.model.get_atom(atom_id)
            atoms_payload.append(
                {
                    "id": atom_id,
                    "element": atom.element,
                    "x": atom.x - left,
                    "y": atom.y - top,
                    "charge": atom.charge,
                    "formal_charge": int(getattr(atom, "formal_charge", atom.charge)),
                    "isotope": atom.isotope,
                    "explicit_h": atom.explicit_h,
                    "mapping": atom.mapping,
                    "is_query": atom.is_query,
                    "is_explicit": atom.is_explicit,
                    "no_implicit": bool(getattr(atom, "no_implicit", False)),
                    "is_coordination_center": getattr(atom, "is_coordination_center", False),
                    "sphere_radius": getattr(atom, "sphere_radius", None),
                    "sphere_color": getattr(atom, "sphere_color", None),
                    "sphere_filled": bool(getattr(atom, "sphere_filled", True)),
                    "sphere_transparent": bool(getattr(atom, "sphere_transparent", False)),
                    "anchor": self._group_anchor_overrides.get(atom_id),
                }
            )

        bonds_payload = []
        for bond in bonds:
                bonds_payload.append(
                    {
                        "a1": bond.a1_id,
                        "a2": bond.a2_id,
                        "order": bond.order,
                        "style": bond.style.value if bond.style is not None else BondStyle.PLAIN.value,
                        "type": bond.style.value if bond.style is not None else BondStyle.PLAIN.value,
                        "stereo": bond.stereo.value if bond.stereo is not None else BondStereo.NONE.value,
                        "is_aromatic": bond.is_aromatic,
                        "display_order": bond.display_order,
                        "is_query": bond.is_query,
                        "ring_id": bond.ring_id,
                        "length_px": bond.length_px,
                        "stroke_px": bond.stroke_px,
                        "color": bond.color,
                        "donor_atom_id": getattr(bond, "donor_atom_id", None),
                        "flex_curve_1": getattr(bond, "flex_curve_1", None),
                        "flex_curve_2": getattr(bond, "flex_curve_2", None),
                    }
                )

        arrows_payload = []
        for item in arrows:
            start = item.start_point()
            end = item.end_point()
            arrows_payload.append(
                {
                    "start": [start.x() - left, start.y() - top],
                    "end": [end.x() - left, end.y() - top],
                    "kind": item.kind(),
                    "curve_factor": item.curve_factor(),
                }
            )

        brackets_payload = []
        for item in brackets:
            rect = item.base_rect()
            brackets_payload.append(
                {
                    "rect": [rect.x() - left, rect.y() - top, rect.width(), rect.height()],
                    "kind": getattr(item, "_kind", "[]"),
                    "padding": getattr(item, "_padding", None),
                }
            )

        texts_payload = []
        for item in texts:
            texts_payload.append(
                {
                    "text": item.toPlainText(),
                    "html": item.toHtml(),
                    "x": item.pos().x() - left,
                    "y": item.pos().y() - top,
                    "rotation": item.rotation(),
                    "font": item.font().toString(),
                    "color": item.defaultTextColor().name(),
                }
            )

        return {
            "atoms": atoms_payload,
            "bonds": bonds_payload,
            "arrows": arrows_payload,
            "brackets": brackets_payload,
            "texts": texts_payload,
        }

    def _paste_selection_payload(self, payload: dict) -> None:
        """Método auxiliar para  paste selection payload.

        Args:
            payload: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atoms = payload.get("atoms", [])
        bonds = payload.get("bonds", [])
        arrows = payload.get("arrows", [])
        brackets = payload.get("brackets", [])
        texts = payload.get("texts", [])
        if not atoms and not bonds and not arrows and not brackets and not texts:
            return

        target = self._last_scene_pos
        if target is None:
            target = self.mapToScene(self.viewport().rect().center())
        dx = target.x()
        dy = target.y()

        ring_map: dict[int, int] = {}

        def map_ring(ring_id: Optional[int]) -> Optional[int]:
            """Método auxiliar para map ring.

            Args:
                ring_id: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
            if ring_id is None:
                return None
            if ring_id not in ring_map:
                ring_map[ring_id] = self.allocate_ring_id()
            return ring_map[ring_id]

        has_undo_items = bool(atoms or bonds or arrows or brackets)
        if has_undo_items:
            self.undo_stack.beginMacro("Paste selection")
        id_map: Dict[int, int] = {}
        for atom_d in atoms:
            cmd = AddAtomCommand(
                self.model,
                self,
                atom_d.get("element", "C"),
                float(atom_d.get("x", 0.0)) + dx,
                float(atom_d.get("y", 0.0)) + dy,
                is_explicit=atom_d.get("is_explicit"),
                charge=atom_d.get("formal_charge", atom_d.get("charge")),
                isotope=atom_d.get("isotope"),
                explicit_h=atom_d.get("explicit_h"),
                mapping=atom_d.get("mapping"),
                is_query=bool(atom_d.get("is_query", False)),
                anchor_override=atom_d.get("anchor"),
                auto_hydrogens=False,
                no_implicit=bool(atom_d.get("no_implicit", False)),
                is_coordination_center=bool(atom_d.get("is_coordination_center", False)),
                sphere_radius=atom_d.get("sphere_radius"),
                sphere_color=atom_d.get("sphere_color"),
                sphere_filled=bool(atom_d.get("sphere_filled", True)),
                sphere_transparent=bool(atom_d.get("sphere_transparent", False)),
            )
            if has_undo_items:
                self.undo_stack.push(cmd)
            else:
                cmd.redo()
            if cmd.atom_id is not None:
                id_map[int(atom_d.get("id"))] = cmd.atom_id

        for bond_d in bonds:
            a1 = id_map.get(int(bond_d.get("a1")))
            a2 = id_map.get(int(bond_d.get("a2")))
            if a1 is None or a2 is None:
                continue
            style = self._parse_bond_style_payload(bond_d)
            stereo = self._parse_bond_stereo_payload(bond_d)
            donor_atom_id = None
            raw_donor = bond_d.get("donor_atom_id")
            if raw_donor is not None:
                try:
                    donor_atom_id = id_map.get(int(raw_donor))
                except Exception:
                    donor_atom_id = None
            if style == BondStyle.COORDINATION:
                donor_atom_id = self._infer_coordination_donor_atom(
                    a1,
                    a2,
                    preferred=donor_atom_id,
                )
            cmd = AddBondCommand(
                self.model,
                self,
                a1,
                a2,
                int(bond_d.get("order", 1)),
                style,
                stereo,
                bool(bond_d.get("is_aromatic", False)),
                display_order=bond_d.get("display_order"),
                length_px=bond_d.get("length_px"),
                stroke_px=bond_d.get("stroke_px"),
                color=bond_d.get("color"),
                ring_id=map_ring(bond_d.get("ring_id")),
                donor_atom_id=donor_atom_id,
                flex_curve_1=bond_d.get("flex_curve_1"),
                flex_curve_2=bond_d.get("flex_curve_2"),
            )
            if has_undo_items:
                self.undo_stack.push(cmd)
            else:
                cmd.redo()

        for arrow_d in arrows:
            start = arrow_d.get("start", [0.0, 0.0])
            end = arrow_d.get("end", [10.0, 0.0])
            kind = arrow_d.get("kind", "forward")
            curve_factor = arrow_d.get("curve_factor")
            try:
                curve_factor_value = float(curve_factor) if curve_factor is not None else None
            except Exception:
                curve_factor_value = None
            start_pt = QPointF(float(start[0]) + dx, float(start[1]) + dy)
            end_pt = QPointF(float(end[0]) + dx, float(end[1]) + dy)
            cmd = AddArrowCommand(self, start_pt, end_pt, kind, curve_factor=curve_factor_value)
            if has_undo_items:
                self.undo_stack.push(cmd)
            else:
                cmd.redo()

        for bracket_d in brackets:
            rect_vals = bracket_d.get("rect", [0.0, 0.0, 10.0, 10.0])
            rect = QRectF(
                float(rect_vals[0]) + dx,
                float(rect_vals[1]) + dy,
                float(rect_vals[2]),
                float(rect_vals[3]),
            )
            kind = bracket_d.get("kind", "[]")
            pair = self._split_bracket_kind(kind)
            if pair:
                for side in pair:
                    cmd = AddBracketCommand(self, rect, side)
                    if has_undo_items:
                        self.undo_stack.push(cmd)
                    else:
                        cmd.redo()
            else:
                cmd = AddBracketCommand(self, rect, kind)
                if has_undo_items:
                    self.undo_stack.push(cmd)
                else:
                    cmd.redo()

        for txt_d in texts:
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
            text_item.setPos(
                float(txt_d.get("x", 0.0)) + dx, float(txt_d.get("y", 0.0)) + dy
            )
            self.scene.addItem(text_item)

        if has_undo_items:
            self.undo_stack.endMacro()
        if ring_map:
            self.refresh_ring_centers()



    def copy_to_clipboard(self) -> None:
        """Método auxiliar para copy to clipboard.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        selected_text_items = self._selected_text_items()
        selected_payload = self._build_selection_payload()
        atom_ids, bonds = self._selected_structure_ids()
        selected_arrows = bool(self._selected_arrow_items())
        selected_brackets = bool(self._selected_bracket_items())
        has_structure_selection = bool(atom_ids or bonds)
        only_text_selection = bool(selected_text_items) and not (
            has_structure_selection or selected_arrows or selected_brackets
        )
        if not self.model.atoms and not selected_text_items:
            return
        mime = QMimeData()
        if selected_payload is not None:
            mime.setData(
                "application/x-chemuson-selection",
                json.dumps(selected_payload).encode("utf-8"),
            )
        if has_structure_selection:
            graph = self._build_selection_graph(atom_ids, bonds)
        elif selected_payload is None:
            graph = self.model if self.model.atoms and not only_text_selection else None
        else:
            graph = None
        smiles = ""
        if graph is not None and graph.atoms:
            try:
                molfile = molgraph_to_molfile(graph)
                smiles = molgraph_to_smiles(graph)
                mime.setData("chemical/x-mdl-molfile", molfile.encode("utf-8"))
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
            background=Qt.GlobalColor.white,
        )
        if image is not None:
            buffer = QBuffer()
            buffer.open(QBuffer.OpenModeFlag.WriteOnly)
            image.save(buffer, "PNG")
            mime.setData("image/png", buffer.data())
            mime.setImageData(image)
            try:
                png_b64 = base64.b64encode(bytes(buffer.data())).decode("ascii")
                html = f'<img src="data:image/png;base64,{png_b64}" alt="{smiles}">'
                mime.setHtml(html)
            except Exception:
                pass
        elif smiles:
            mime.setText(smiles)
        svg_data = self._render_scene_svg(selected_only=has_selection)
        if svg_data:
            mime.setData("image/svg+xml", svg_data)

        QApplication.clipboard().setMimeData(mime)

    def paste_from_clipboard(self) -> None:
        """Método auxiliar para paste from clipboard.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        clipboard = QApplication.clipboard()
        mime = clipboard.mimeData()
        if mime is None:
            return

        if mime.hasFormat("application/x-chemuson-selection"):
            try:
                payload = json.loads(bytes(mime.data("application/x-chemuson-selection")).decode("utf-8"))
                self._paste_selection_payload(payload)
                return
            except Exception:
                pass

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
        """Método auxiliar para  paste text items.

        Args:
            payload: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  render scene bounds.

        Args:
            selected_only: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        rect: Optional[QRectF] = None

        def extend(candidate: QRectF) -> None:
            """Método auxiliar para extend.

            Args:
                candidate: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
            nonlocal rect
            if not candidate.isValid() or candidate.isNull():
                return
            if rect is None:
                rect = candidate
            else:
                rect = rect.united(candidate)

        def extend_atom_bounds(atom_id: int) -> None:
            """Método auxiliar para extend atom bounds.

            Args:
                atom_id: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
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
                if (
                    isinstance(item, TextAnnotationItem)
                    and item.isVisible()
                    and not bool(item.data(NUMBERING_TEXT_ROLE))
                ):
                    extend(item.sceneBoundingRect())
            if self.state.numbering_enabled and self.state.numbering_include_in_export:
                for item in self._numbering_overlay_items:
                    if item is None or item.scene() is not self.scene:
                        continue
                    if item.isVisible():
                        extend(item.sceneBoundingRect())

        if rect is None:
            return None
        pad = max(1.0, self.drawing_style.stroke_px)
        return rect.adjusted(-pad, -pad, pad, pad)

    def _hidden_render_items(self) -> list:
        """Método auxiliar para  hidden render items.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        if not self.state.numbering_include_in_export:
            items.extend(self._numbering_overlay_items)
        return [
            item
            for item in items
            if item is not None and item.scene() is self.scene
        ]

    def _with_hidden_render_items(self, render_fn):
        """Método auxiliar para  with hidden render items.

        Args:
            render_fn: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  with hidden unselected.

        Args:
            render_fn: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        self,
        scale: float = 1.0,
        selected_only: bool = False,
        background: Optional[QColor] = None,
    ) -> Optional[QImage]:
        """Método auxiliar para  render scene image.

        Args:
            scale: Descripción del parámetro.
            selected_only: Descripción del parámetro.
            background: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        rect = self._render_scene_bounds(selected_only=selected_only)
        if rect is None:
            return None
        scale = max(1.0, float(scale))
        width = max(1, math.ceil(rect.width() * scale))
        height = max(1, math.ceil(rect.height() * scale))

        def render():
            """Método auxiliar para render.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
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
            if background is not None:
                flattened = QImage(trimmed.width(), trimmed.height(), QImage.Format.Format_RGB32)
                flattened.fill(background)
                painter = QPainter(flattened)
                painter.setRenderHint(QPainter.RenderHint.Antialiasing)
                painter.setRenderHint(QPainter.RenderHint.TextAntialiasing)
                painter.drawImage(0, 0, trimmed)
                painter.end()
                self._apply_image_dpi(flattened, scale)
                return flattened
            self._apply_image_dpi(trimmed, scale)
            return trimmed

        if selected_only:
            return self._with_hidden_unselected(
                lambda: self._with_hidden_render_items(render)
            )
        return self._with_hidden_render_items(render)

    def _apply_image_dpi(self, image: QImage, scale: float) -> None:
        """Método auxiliar para  apply image dpi.

        Args:
            image: Descripción del parámetro.
            scale: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        base_dpi_x = self.logicalDpiX() or 96.0
        base_dpi_y = self.logicalDpiY() or 96.0
        dpi_x = base_dpi_x * scale
        dpi_y = base_dpi_y * scale
        dpm_x = max(1, round(dpi_x * 1000.0 / 25.4))
        dpm_y = max(1, round(dpi_y * 1000.0 / 25.4))
        image.setDotsPerMeterX(dpm_x)
        image.setDotsPerMeterY(dpm_y)

    def _trim_transparent_image(self, image: QImage) -> Optional[QImage]:
        """Método auxiliar para  trim transparent image.

        Args:
            image: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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

    def _render_scene_svg(self, selected_only: bool = False) -> Optional[bytes]:
        """Método auxiliar para  render scene svg.

        Args:
            selected_only: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if QSvgGenerator is None:
            return None
        rect = self._render_scene_bounds(selected_only=selected_only)
        if rect is None:
            return None
        width = max(1, math.ceil(rect.width()))
        height = max(1, math.ceil(rect.height()))

        def render():
            """Método auxiliar para render.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
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

        if selected_only:
            return self._with_hidden_unselected(
                lambda: self._with_hidden_render_items(render)
            )
        return self._with_hidden_render_items(render)

    def _insert_molgraph(self, graph: MolGraph) -> None:
        """Método auxiliar para  insert molgraph.

        Args:
            graph: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
                charge=atom.charge,
                no_implicit=bool(getattr(atom, "no_implicit", False)),
                auto_hydrogens=False,
                is_coordination_center=getattr(atom, "is_coordination_center", False),
                sphere_radius=getattr(atom, "sphere_radius", None),
                sphere_color=getattr(atom, "sphere_color", None),
                sphere_filled=bool(getattr(atom, "sphere_filled", True)),
                sphere_transparent=bool(getattr(atom, "sphere_transparent", False)),
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
                stroke_px=bond.stroke_px,
                color=bond.color,
                donor_atom_id=id_map.get(getattr(bond, "donor_atom_id", -1)),
                flex_curve_1=getattr(bond, "flex_curve_1", None),
                flex_curve_2=getattr(bond, "flex_curve_2", None),
            )
            self.undo_stack.push(cmd)
        self.undo_stack.endMacro()
        if any(bond.is_aromatic for bond in self.model.bonds.values()):
            self._kekulize_aromatic_bonds()

    def _pick_template_attachment_atom_id(self, graph: MolGraph) -> Optional[int]:
        """Elige un átomo de la plantilla para conectar al átomo ancla.

        Se priorizan vértices terminales (grado mínimo) y, entre ellos, el más
        alejado del centro para favorecer una orientación externa.
        """
        if not graph.atoms:
            return None
        degrees: Dict[int, int] = {atom_id: 0 for atom_id in graph.atoms.keys()}
        for bond in graph.bonds.values():
            if bond.a1_id in degrees:
                degrees[bond.a1_id] += 1
            if bond.a2_id in degrees:
                degrees[bond.a2_id] += 1
        min_degree = min(degrees.values()) if degrees else 0
        candidate_ids = [atom_id for atom_id, degree in degrees.items() if degree == min_degree]
        if len(candidate_ids) == 1:
            return candidate_ids[0]
        xs = [atom.x for atom in graph.atoms.values()]
        ys = [atom.y for atom in graph.atoms.values()]
        center_x = (min(xs) + max(xs)) / 2.0
        center_y = (min(ys) + max(ys)) / 2.0
        return max(
            candidate_ids,
            key=lambda atom_id: (
                (graph.get_atom(atom_id).x - center_x) ** 2
                + (graph.get_atom(atom_id).y - center_y) ** 2
            ),
        )

    def _insert_molgraph_at(
        self,
        graph: MolGraph,
        target: QPointF,
        attach_to_atom_id: Optional[int] = None,
    ) -> None:
        """Método auxiliar para  insert molgraph at.

        Args:
            graph: Descripción del parámetro.
            target: Descripción del parámetro.
            attach_to_atom_id: Átomo del lienzo al que se conectará la plantilla.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if not graph.atoms:
            return
        xs = [atom.x for atom in graph.atoms.values()]
        ys = [atom.y for atom in graph.atoms.values()]
        min_x, max_x = min(xs), max(xs)
        min_y, max_y = min(ys), max(ys)
        center_x = (min_x + max_x) / 2
        center_y = (min_y + max_y) / 2
        attach_template_id = None
        if attach_to_atom_id is not None and attach_to_atom_id in self.model.atoms:
            attach_template_id = self._pick_template_attachment_atom_id(graph)
        if attach_template_id is not None:
            attach_atom = graph.get_atom(attach_template_id)
            anchor = self.model.get_atom(attach_to_atom_id)
            anchor_pos = QPointF(anchor.x, anchor.y)
            attach_target = self._compute_default_bond_endpoint(anchor_pos, attach_to_atom_id)
            dx = attach_target.x() - attach_atom.x
            dy = attach_target.y() - attach_atom.y
        else:
            dx = target.x() - center_x
            dy = target.y() - center_y

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
                charge=atom.charge,
                no_implicit=bool(getattr(atom, "no_implicit", False)),
                auto_hydrogens=False,
                is_coordination_center=getattr(atom, "is_coordination_center", False),
                sphere_radius=getattr(atom, "sphere_radius", None),
                sphere_color=getattr(atom, "sphere_color", None),
                sphere_filled=bool(getattr(atom, "sphere_filled", True)),
                sphere_transparent=bool(getattr(atom, "sphere_transparent", False)),
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
                stroke_px=bond.stroke_px,
                color=bond.color,
                donor_atom_id=id_map.get(getattr(bond, "donor_atom_id", -1)),
                flex_curve_1=getattr(bond, "flex_curve_1", None),
                flex_curve_2=getattr(bond, "flex_curve_2", None),
            )
            self.undo_stack.push(cmd)
        if attach_to_atom_id is not None and attach_template_id is not None:
            template_new_id = id_map.get(attach_template_id)
            if template_new_id is not None:
                cmd = AddBondCommand(
                    self.model,
                    self,
                    attach_to_atom_id,
                    template_new_id,
                    1,
                    BondStyle.PLAIN,
                    BondStereo.NONE,
                    is_aromatic=False,
                )
                self.undo_stack.push(cmd)
        self.undo_stack.endMacro()
        if any(bond.is_aromatic for bond in self.model.bonds.values()):
            self._kekulize_aromatic_bonds()

    def begin_template_insert_mode(self, graph: MolGraph, label: str) -> None:
        """Activa la inserción por clic para una plantilla."""
        self._pending_template_graph = graph
        self._pending_template_label = label
        self.setCursor(Qt.CursorShape.CrossCursor)

    def cancel_template_insert_mode(self) -> None:
        """Cancela la inserción pendiente de plantillas."""
        self._pending_template_graph = None
        self._pending_template_label = None
        if self._space_panning:
            self.setCursor(Qt.CursorShape.OpenHandCursor)
        else:
            self.unsetCursor()

    def _insert_ring_template(self, scene_pos: QPointF) -> None:
        """Método auxiliar para  insert ring template.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        template = self.state.active_ring_template
        if not template:
            return
        anomeric = (self.state.active_ring_anomeric or "beta").lower()
        anomeric_up = anomeric != "alpha"
        graph = None
        if template == "haworth":
            graph = build_haworth_template(
                self.state.bond_length, anomeric_up=anomeric_up, bold_front=True
            )
        elif template == "chair":
            graph = build_chair_template(
                self.state.bond_length, anomeric_up=anomeric_up, bold_front=True
            )
        if graph is None:
            return
        self._insert_molgraph_at(graph, scene_pos)

    def _insert_image_from_clipboard(self, mime: QMimeData) -> None:
        """Método auxiliar para  insert image from clipboard.

        Args:
            mime: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Añade atom item.

        Args:
            atom: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        self.recompute_numbering()

    def update_atom_item(self, atom_id: int, x: float, y: float) -> None:
        """Actualiza atom item.

        Args:
            atom_id: Descripción del parámetro.
            x: Descripción del parámetro.
            y: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = self.atom_items.get(atom_id)
        if item is not None:
            atom_ref = self.model.atoms.get(atom_id)
            if atom_ref is not None and hasattr(item, "atom"):
                item.atom = atom_ref
            item.setPos(x, y)
        self._update_electron_dots_for_atom(atom_id)
        self._update_wavy_anchors_for_atom(atom_id)

    def update_atom_item_element(
        self,
        atom_id: int,
        element: str,
        is_explicit: Optional[bool] = None,
    ) -> None:
        """Actualiza atom item element.

        Args:
            atom_id: Descripción del parámetro.
            element: Descripción del parámetro.
            is_explicit: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = self.atom_items.get(atom_id)
        if item is not None:
            atom_ref = self.model.atoms.get(atom_id)
            if atom_ref is not None and hasattr(item, "atom"):
                item.atom = atom_ref
            if is_explicit is None:
                is_explicit = self.model.get_atom(atom_id).is_explicit
            item.set_element(element, is_explicit=is_explicit)
            self._refresh_atom_label(atom_id)
        self._update_electron_dots_for_atom(atom_id)
        self._update_wavy_anchors_for_atom(atom_id)

    def update_atom_item_charge(self, atom_id: int, charge: int) -> None:
        """Actualiza atom item charge.

        Args:
            atom_id: Descripción del parámetro.
            charge: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = self.atom_items.get(atom_id)
        if item is not None:
            atom_ref = self.model.atoms.get(atom_id)
            if atom_ref is not None and hasattr(item, "atom"):
                item.atom = atom_ref
            item.set_charge(charge)
            self._refresh_atom_label(atom_id)
        self._sync_selection_from_scene()

    def _label_font(self) -> QFont:
        """Método auxiliar para  label font.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para label font.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        return self._label_font()

    def apply_label_font(self, font: QFont) -> None:
        """Aplica label font.

        Args:
            font: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para refresh label fonts.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        font = self._label_font()
        for item in self.atom_items.values():
            item.set_label_font(font)
        self.refresh_atom_labels()
        self.recompute_numbering()

    def set_use_element_colors(self, use_element_colors: bool) -> None:
        """Actualiza el estado de use element colors.

        Args:
            use_element_colors: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self.state.use_element_colors = use_element_colors
        for item in self.atom_items.values():
            item.set_use_element_colors(use_element_colors)
        self.refresh_atom_labels()

    def refresh_atom_labels(self, atom_ids: Optional[Iterable[int]] = None) -> None:
        """Método auxiliar para refresh atom labels.

        Args:
            atom_ids: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if atom_ids is None:
            atom_ids = list(self.atom_items.keys())
        for atom_id in atom_ids:
            self._refresh_atom_label(atom_id)

    def _refresh_atom_label(self, atom_id: int) -> None:
        """Método auxiliar para  refresh atom label.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = self.atom_items.get(atom_id)
        atom = self.model.atoms.get(atom_id)
        if item is None or atom is None:
            return
        label, anchor, offset = self._build_atom_label(atom)
        item.set_display_label(label, anchor)
        item.set_label_offset(offset)
        self._update_bond_label_shrinks({atom_id})

    def _build_atom_label(self, atom) -> tuple[str, Optional[str], QPointF]:
        """Método auxiliar para  build atom label.

        Args:
            atom: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        label = atom.element
        anchor: Optional[str] = None
        anchor_override = self._group_anchor_overrides.get(atom.id)
        is_coordination_center = bool(getattr(atom, "is_coordination_center", False))
        if atom.element in ELEMENT_SYMBOLS:
            if (
                not is_coordination_center
                and atom.element == "C"
                and not (
                self.state.show_implicit_carbons or atom.is_explicit
                )
            ):
                return label, None, QPointF(0.0, 0.0)
            implicit_h = self._implicit_hydrogen_count(atom.id, atom.element)
            if (
                implicit_h > 0
                and atom.element in {"O", "S"}
                and self._atom_degree(atom.id) >= 2
            ):
                implicit_h = 0
            if (
                implicit_h > 0
                and atom.element != "H"
                and not self.state.show_implicit_hydrogens
                and not is_coordination_center
            ):
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

    def _atom_degree(self, atom_id: int) -> int:
        """Método auxiliar para  atom degree.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        return sum(
            1
            for bond in self.model.bonds.values()
            if bond.a1_id == atom_id or bond.a2_id == atom_id
        )

    def _reflow_group_label(
        self, label: str, atom_id: int, anchor_override: Optional[str] = None
    ) -> tuple[str, Optional[str]]:
        """Método auxiliar para  reflow group label.

        Args:
            label: Descripción del parámetro.
            atom_id: Descripción del parámetro.
            anchor_override: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  move anchor in label.

        Args:
            label: Descripción del parámetro.
            anchor: Descripción del parámetro.
            anchor_first: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para get anchor override.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        return self._group_anchor_overrides.get(atom_id)

    def set_anchor_override(self, atom_id: int, anchor: Optional[str]) -> None:
        """Actualiza el estado de anchor override.

        Args:
            atom_id: Descripción del parámetro.
            anchor: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if anchor:
            self._group_anchor_overrides[atom_id] = anchor
        else:
            self._group_anchor_overrides.pop(atom_id, None)

    def _group_label_tokens(self, label: str) -> Optional[list[tuple[str, str]]]:
        """Método auxiliar para  group label tokens.

        Args:
            label: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  implicit hydrogen count.

        Args:
            atom_id: Descripción del parámetro.
            element: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if element not in IMPLICIT_H_ELEMENTS:
            return 0
        try:
            return int(max(0, self.model.implicit_h_count(atom_id)))
        except Exception:
            return 0

    def _label_open_direction(self, atom_id: int) -> QPointF:
        """Método auxiliar para  label open direction.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        def collect_vectors(skip_hydrogen: bool) -> list[tuple[float, float]]:
            """Método auxiliar para collect vectors.

            Args:
                skip_hydrogen: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
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
        """Método auxiliar para  label offset.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  neighbor angles deg.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  atom hybridization.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  angle of vector.

        Args:
            vx: Descripción del parámetro.
            vy: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        return math.degrees(math.atan2(-vy, vx))

    def _normalize_angle(self, angle: float) -> float:
        """Método auxiliar para  normalize angle.

        Args:
            angle: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
            """Método auxiliar para norm deg.

            Args:
                a: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
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
        """Método auxiliar para  clear implicit h overlays.

        Args:
            atom_ids: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  refresh implicit h overlays.

        Args:
            atom_ids: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  local bond length.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  update bond label shrinks.

        Args:
            atom_ids: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  label shrink for atom.

        Args:
            atom_id: Descripción del parámetro.
            ux: Descripción del parámetro.
            uy: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = self.atom_items.get(atom_id)
        atom = self.model.atoms.get(atom_id)
        if (
            atom is not None
            and bool(getattr(atom, "is_coordination_center", False))
            and not bool(getattr(atom, "sphere_transparent", False))
        ):
            if item is not None and hasattr(item, "_coordination_draw_radius"):
                try:
                    return max(0.0, float(item._coordination_draw_radius()) + 1.5)
                except Exception:
                    pass
            configured = getattr(atom, "sphere_radius", None)
            return max(0.0, float(configured) if configured is not None else 16.0)
        if item is None or not item.label.isVisible():
            return 0.0
        rect = item.label.mapRectToParent(item.label.boundingRect())
        if rect.isNull() or rect.width() <= 0 or rect.height() <= 0:
            return 0.0
        # Increased padding to clear label characters (especially "C")
        # 6.0px provides a comfortable margin for standard font sizes.
        pad = 6.0
        if atom is not None and atom.element in ELEMENT_SYMBOLS and atom.element not in {"C", "H"}:
            pad = 3.5
        rect = rect.adjusted(-pad, -pad, pad, pad)
        distance = self._ray_ellipse_distance(rect, ux, uy)
        return distance if distance is not None else 0.0

    def _configure_bond_rendering(self, bond: Bond, item: BondItem) -> None:
        """Método auxiliar para  configure bond rendering.

        Args:
            bond: Descripción del parámetro.
            item: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        prefer_full_length = self.state.show_implicit_carbons and self.state.show_implicit_hydrogens
        is_aromatic_ring = bond.is_aromatic and bond.ring_id is not None
        a1 = self.model.get_atom(bond.a1_id)
        a2 = self.model.get_atom(bond.a2_id)
        has_hetero = False
        if a1 is not None and a2 is not None:
            has_hetero = (a1.element != "C") or (a2.element != "C")
        # Estilo base: vértice tipo C=C (línea sigma central + línea pi desplazada).
        # En heteroátomos la línea pi debe ser completa (sin recorte), pero
        # manteniendo la misma orientación del doble enlace tipo C=C.
        effective_order = bond.display_order if bond.display_order is not None else bond.order
        is_plain_double = (
            bond.style == BondStyle.PLAIN and effective_order == 2 and not is_aromatic_ring
        )
        def _neighbor_bonds(atom_id: int) -> list[Bond]:
            neighbors: list[Bond] = []
            for other in self.model.bonds.values():
                if other.id == bond.id:
                    continue
                if other.style == BondStyle.COORDINATION:
                    continue
                if other.a1_id != atom_id and other.a2_id != atom_id:
                    continue
                other_id = other.a2_id if other.a1_id == atom_id else other.a1_id
                other_atom = self.model.atoms.get(other_id)
                if other_atom is not None and other_atom.element == "H":
                    continue
                neighbors.append(other)
            return neighbors

        def _is_terminal_like(atom_id: int) -> bool:
            total = 0
            for other in self.model.bonds.values():
                if other.style == BondStyle.COORDINATION:
                    continue
                if other.a1_id != atom_id and other.a2_id != atom_id:
                    continue
                other_id = other.a2_id if other.a1_id == atom_id else other.a1_id
                other_atom = self.model.atoms.get(other_id)
                if other_atom is not None and other_atom.element == "H":
                    continue
                total += 1
            return total <= 1

        def _needs_centered_double() -> bool:
            if not is_plain_double or a1 is None or a2 is None:
                return False
            candidates = ((a1.id, a2.id), (a2.id, a1.id))
            for center_id, terminal_id in candidates:
                neighbors = _neighbor_bonds(center_id)
                if len(neighbors) != 2:
                    continue
                if any((nb.display_order if nb.display_order is not None else nb.order) != 1 for nb in neighbors):
                    continue
                if _is_terminal_like(terminal_id):
                    return True
            return False

        if has_hetero and is_plain_double:
            prefer_full_length = True
        symmetric_double = _needs_centered_double()
        item.set_multibond_rendering(prefer_full_length, symmetric_double)

    @staticmethod
    def _ray_ellipse_distance(rect: QRectF, ux: float, uy: float) -> Optional[float]:
        # Ray P = t * D, where D = (ux, uy). Origin (0,0).
        # Ellipse defined by rect: Center C, radii a, b.
        # Equation: ((x - cx)/a)^2 + ((y - cy)/b)^2 = 1
        # Substitute x = t*ux, y = t*uy
        
        """Método auxiliar para  ray ellipse distance.

        Args:
            rect: Descripción del parámetro.
            ux: Descripción del parámetro.
            uy: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  prefer prefix h.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        direction = self._label_open_direction(atom_id)
        return direction.x() < -0.2

    def remove_atom_item(self, atom_id: int) -> None:
        """Elimina atom item.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = self.atom_items.pop(atom_id, None)
        if item is not None:
            self.scene.removeItem(item)
        self._clear_implicit_h_overlays([atom_id])
        self._group_anchor_overrides.pop(atom_id, None)
        if self.bond_anchor_id == atom_id:
            self.bond_anchor_id = None
        if self.hovered_atom_id == atom_id:
            self.hovered_atom_id = None
        self.recompute_numbering()

    def add_bond_item(self, bond) -> None:
        """Añade bond item.

        Args:
            bond: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        ring_pairs = self._compute_ring_bond_pairs()
        ring_center = self._aromatic_ring_center_for_bond(bond) if bond.is_aromatic else None
        if ring_center is None and bond.ring_id is not None:
            ring_center = self._ring_centers.get(bond.ring_id)
        item.set_ring_context(ring_center)
        item.set_bond_in_ring(self._bond_in_ring_for_pairs(bond, ring_pairs))
        item.set_offset_sign(self._bond_offset_sign(bond))
        self._configure_bond_rendering(bond, item)
        self._set_bond_item_join_context(bond, item)
        trim_start = self._bond_endpoint_trim(bond, bond.a1_id)
        trim_end = self._bond_endpoint_trim(bond, bond.a2_id)
        item.set_endpoint_trim(trim_start, trim_end)
        extend_start = self._bond_endpoint_extend(bond, bond.a1_id)
        extend_end = self._bond_endpoint_extend(bond, bond.a2_id)
        item.set_endpoint_extend(extend_start, extend_end)
        item.update_positions(atom1, atom2)
        self.scene.addItem(item)
        self.bond_items[bond.id] = item
        self._refresh_atom_label(bond.a1_id)
        self._refresh_atom_label(bond.a2_id)
        self._refresh_bond_ring_flags(ring_pairs)
        self.recompute_numbering()

    def add_arrow_item(
        self,
        start: QPointF,
        end: QPointF,
        kind: str,
        curve_factor: float | None = None,
    ) -> ArrowItem:
        """Añade arrow item.

        Args:
            start: Descripción del parámetro.
            end: Descripción del parámetro.
            kind: Descripción del parámetro.
            curve_factor: Curvatura normalizada para flechas curvas.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = ArrowItem(
            start,
            end,
            kind=kind,
            curve_factor=curve_factor,
            style=self.drawing_style,
        )
        self.scene.addItem(item)
        self.arrow_items.append(item)
        return item

    def readd_arrow_item(
        self,
        item: ArrowItem,
        start: QPointF,
        end: QPointF,
        kind: str,
        curve_factor: float | None = None,
    ) -> None:
        """Método auxiliar para readd arrow item.

        Args:
            item: Descripción del parámetro.
            start: Descripción del parámetro.
            end: Descripción del parámetro.
            kind: Descripción del parámetro.
            curve_factor: Curvatura normalizada para flechas curvas.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item.set_kind(kind)
        item.update_positions(start, end, curve_factor=curve_factor)
        if item.scene() is not self.scene:
            self.scene.addItem(item)
        if item not in self.arrow_items:
            self.arrow_items.append(item)

    def remove_arrow_item(self, item: ArrowItem) -> None:
        """Elimina arrow item.

        Args:
            item: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if item in self.arrow_items:
            self.arrow_items.remove(item)
        if item.scene() is self.scene:
            self.scene.removeItem(item)

    def add_bracket_item(self, rect: QRectF, kind: str) -> BracketItem:
        """Añade bracket item.

        Args:
            rect: Descripción del parámetro.
            kind: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para readd bracket item.

        Args:
            item: Descripción del parámetro.
            rect: Descripción del parámetro.
            kind: Descripción del parámetro.
            padding: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item.set_rect(rect, padding=padding)
        item._kind = kind
        if item.scene() is not self.scene:
            self.scene.addItem(item)
        if item not in self.bracket_items:
            self.bracket_items.append(item)

    def remove_bracket_item(self, item: BracketItem) -> None:
        """Elimina bracket item.

        Args:
            item: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if item in self.bracket_items:
            self.bracket_items.remove(item)
        if item.scene() is self.scene:
            self.scene.removeItem(item)

    @staticmethod
    def _split_bracket_kind(kind: str) -> Optional[tuple[str, str]]:
        """Método auxiliar para  split bracket kind.

        Args:
            kind: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if kind in {"()", "[]", "{}"} and len(kind) == 2:
            return kind[0], kind[1]
        return None

    def _ensure_bracket_preview(self) -> QGraphicsRectItem:
        """Método auxiliar para  ensure bracket preview.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  clear bracket preview.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._bracket_preview is not None:
            self.scene.removeItem(self._bracket_preview)
            self._bracket_preview = None
        self._bracket_drag_start = None

    def update_bond_item(self, bond_id: int) -> None:
        """Actualiza bond item.

        Args:
            bond_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        bond = self.model.get_bond(bond_id)
        atom1 = self.model.get_atom(bond.a1_id)
        atom2 = self.model.get_atom(bond.a2_id)
        item = self.bond_items.get(bond_id)
        if item is not None:
            ring_pairs = self._compute_ring_bond_pairs()
            ring_center = self._aromatic_ring_center_for_bond(bond) if bond.is_aromatic else None
            if ring_center is None and bond.ring_id is not None:
                ring_center = self._ring_centers.get(bond.ring_id)
            item.set_ring_context(ring_center)
            item.set_bond_in_ring(self._bond_in_ring_for_pairs(bond, ring_pairs))
            item.set_offset_sign(self._bond_offset_sign(bond))
            self._configure_bond_rendering(bond, item)
            self._set_bond_item_join_context(bond, item)
            trim_start = self._bond_endpoint_trim(bond, bond.a1_id)
            trim_end = self._bond_endpoint_trim(bond, bond.a2_id)
            item.set_endpoint_trim(trim_start, trim_end)
            extend_start = self._bond_endpoint_extend(bond, bond.a1_id)
            extend_end = self._bond_endpoint_extend(bond, bond.a2_id)
            item.set_endpoint_extend(extend_start, extend_end)
            item.set_bond(bond, atom1, atom2)
            item.set_style(self.drawing_style, atom1, atom2)
        self._refresh_bond_ring_flags(ring_pairs if item is not None else None)
        self._refresh_atom_label(bond.a1_id)
        self._refresh_atom_label(bond.a2_id)

    def _bond_render_width(self, bond: Bond) -> float:
        """Método auxiliar para  bond render width.

        Args:
            bond: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        base = bond.stroke_px if bond.stroke_px is not None else self.drawing_style.stroke_px
        if bond.style == BondStyle.BOLD:
            return max(base * 2.2, base + 1.0)
        if bond.style == BondStyle.WEDGE:
            scale = base / self.drawing_style.stroke_px if self.drawing_style.stroke_px > 1e-6 else 1.0
            width = self.drawing_style.wedge_width_px * (0.72 + 0.28 * math.sqrt(max(scale, 1e-6)))
            return max(width, base * 2.3)
        if bond.style == BondStyle.HASHED:
            scale = base / self.drawing_style.stroke_px if self.drawing_style.stroke_px > 1e-6 else 1.0
            return max(self.drawing_style.hash_stroke_px * scale, base * 0.85)
        return base

    def _bond_neighbor_vectors(
        self,
        bond: Bond,
        atom_id: int,
        *,
        simple_only: bool = False,
    ) -> list[tuple[float, float, float, float, float]]:
        """Devuelve vecinos con dirección, ancho y centro de extremo renderizado."""
        anchor = self.model.get_atom(atom_id)
        if anchor is None:
            return []
        neighbors: list[tuple[float, float, float, float, float]] = []
        for other in self.model.bonds.values():
            if other.id == bond.id:
                continue
            if simple_only and other.style not in {BondStyle.PLAIN, BondStyle.BOLD}:
                continue
            if other.a1_id == atom_id:
                other_atom = self.model.get_atom(other.a2_id)
            elif other.a2_id == atom_id:
                other_atom = self.model.get_atom(other.a1_id)
            else:
                continue
            if other_atom is None:
                continue
            dx = other_atom.x - anchor.x
            dy = other_atom.y - anchor.y
            length = math.hypot(dx, dy)
            if length <= 1e-6:
                continue
            nux = dx / length
            nuy = dy / length
            nwidth = self._bond_render_width(other)
            shrink = self._label_shrink_for_atom(atom_id, nux, nuy)
            trim = self._bond_endpoint_trim(other, atom_id)
            advance = max(0.0, shrink + trim)
            end_x = anchor.x + nux * advance
            end_y = anchor.y + nuy * advance

            cap_extra = 0.0
            if self.drawing_style.cap_style in (
                Qt.PenCapStyle.RoundCap,
                Qt.PenCapStyle.SquareCap,
            ):
                cap_extra = nwidth * 0.5

            extend = 0.0
            if self.drawing_style.cap_style == Qt.PenCapStyle.FlatCap and advance <= 1e-6:
                extend = self._bond_endpoint_extend(other, atom_id)

            edge_cx = end_x - nux * (cap_extra + extend)
            edge_cy = end_y - nuy * (cap_extra + extend)
            neighbors.append((nux, nuy, nwidth, edge_cx, edge_cy))
        return neighbors

    def _set_bond_item_join_context(self, bond: Bond, item: BondItem) -> None:
        """Inyecta contexto de enlaces vecinos para geometrías dependientes de unión."""
        if bond.style != BondStyle.WEDGE:
            item.set_wedge_join_neighbors([], [])
            return
        start_neighbors = self._bond_neighbor_vectors(bond, bond.a1_id)
        # El extremo ancho de la cuña se integra contra enlaces planos:
        # - Plain (incluye dobles/aromáticos por order/display_order)
        # - Bold
        # Esto permite vértices limpios en uniones cuña->bold.
        end_neighbors = self._bond_neighbor_vectors(bond, bond.a2_id, simple_only=True)
        item.set_wedge_join_neighbors(start_neighbors, end_neighbors)

    def _bond_endpoint_extend(self, bond: Bond, atom_id: int) -> float:
        """Método auxiliar para  bond endpoint extend.

        Args:
            bond: Descripción del parámetro.
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self.drawing_style.cap_style != Qt.PenCapStyle.FlatCap:
            return 0.0
        if self._atom_degree(atom_id) != 2:
            return 0.0
        atom = self.model.get_atom(atom_id)
        if self._atom_label_visible(atom):
            return 0.0
        other_bond = None
        for candidate in self.model.bonds.values():
            if candidate.id == bond.id:
                continue
            if candidate.a1_id == atom_id or candidate.a2_id == atom_id:
                other_bond = candidate
                break
        if other_bond is None:
            return 0.0

        def unit_vector(b: Bond) -> tuple[float, float]:
            """Método auxiliar para unit vector.

            Args:
                b: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
            other_id = b.a2_id if b.a1_id == atom_id else b.a1_id
            other_atom = self.model.get_atom(other_id)
            dx = other_atom.x - atom.x
            dy = other_atom.y - atom.y
            length = math.hypot(dx, dy) or 1.0
            return dx / length, dy / length

        v1x, v1y = unit_vector(bond)
        v2x, v2y = unit_vector(other_bond)
        dot = max(-1.0, min(1.0, v1x * v2x + v1y * v2y))
        theta = math.acos(dot)
        if theta <= 1e-3:
            return 0.0
        width = max(self._bond_render_width(bond), self._bond_render_width(other_bond))
        extension = (width * 0.5) / math.tan(theta / 2.0)
        return max(0.0, min(extension, width * 2.0))

    def _atom_label_visible(self, atom) -> bool:
        """Método auxiliar para  atom label visible.

        Args:
            atom: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if atom.is_explicit:
            return True
        if atom.element == "C":
            return self.state.show_implicit_carbons
        if atom.element == "H":
            return self.state.show_implicit_hydrogens
        return True

    def _bond_endpoint_trim(self, bond: Bond, atom_id: int) -> float:
        """Método auxiliar para  bond endpoint trim.

        Args:
            bond: Descripción del parámetro.
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        # Keep stereo wedges visually full-length at junctions (ChemDraw-like).
        # Label collision avoidance is handled separately via label_shrink.
        if bond.style in (BondStyle.WEDGE, BondStyle.HASHED, BondStyle.COORDINATION):
            return 0.0
        if self._atom_degree(atom_id) < 2:
            return 0.0
        width = self._bond_render_width(bond)
        max_other = 0.0
        for other in self.model.bonds.values():
            if other.id == bond.id:
                continue
            if other.a1_id == atom_id or other.a2_id == atom_id:
                max_other = max(max_other, self._bond_render_width(other))
        if max_other <= 0.0:
            return 0.0
        if width <= max_other:
            return 0.0
        return max(0.0, (width - max_other) * 0.5)

    def update_bond_items_for_atoms(self, atom_ids: set[int]) -> None:
        """Actualiza bond items for atoms.

        Args:
            atom_ids: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._clear_trackball_reference_if_desynced(atom_ids)
        for bond in self.model.bonds.values():
            if bond.a1_id in atom_ids or bond.a2_id in atom_ids:
                self.update_bond_item(bond.id)
        self._refresh_implicit_h_overlays(atom_ids)
        self._update_aromatic_circles_for_atoms(atom_ids)
        self.recompute_numbering()

    def remove_bond_item(self, bond_id: int) -> None:
        """Elimina bond item.

        Args:
            bond_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = self.bond_items.pop(bond_id, None)
        if item is not None:
            a1_id = item.a1_id
            a2_id = item.a2_id
            self.scene.removeItem(item)
            self._refresh_atom_label(a1_id)
            self._refresh_atom_label(a2_id)
        self._refresh_bond_ring_flags()
        self.recompute_numbering()

    def allocate_ring_id(self) -> int:
        """Método auxiliar para allocate ring id.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        ring_id = self._next_ring_id
        self._next_ring_id += 1
        return ring_id

    def register_ring_center(self, ring_id: int, center: tuple[float, float]) -> None:
        """Método auxiliar para register ring center.

        Args:
            ring_id: Descripción del parámetro.
            center: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._ring_centers[ring_id] = QPointF(center[0], center[1])

    def unregister_ring_center(self, ring_id: int) -> None:
        """Método auxiliar para unregister ring center.

        Args:
            ring_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._ring_centers.pop(ring_id, None)

    def _bond_offset_sign(self, bond) -> int:
        """Método auxiliar para  bond offset sign.

        Args:
            bond: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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

        if bond.is_aromatic:
            aromatic_center = self._aromatic_ring_center_for_bond(bond)
            if aromatic_center is not None:
                vx = aromatic_center.x() - mid.x()
                vy = aromatic_center.y() - mid.y()
                return 1 if (nx * vx + ny * vy) >= 0 else -1
            neighbor_points = []
            for other_bond in self.model.bonds.values():
                if not other_bond.is_aromatic or other_bond.id == bond.id:
                    continue
                if other_bond.a1_id == atom1.id and other_bond.a2_id != atom2.id:
                    other = self.model.get_atom(other_bond.a2_id)
                elif other_bond.a2_id == atom1.id and other_bond.a1_id != atom2.id:
                    other = self.model.get_atom(other_bond.a1_id)
                elif other_bond.a1_id == atom2.id and other_bond.a2_id != atom1.id:
                    other = self.model.get_atom(other_bond.a2_id)
                elif other_bond.a2_id == atom2.id and other_bond.a1_id != atom1.id:
                    other = self.model.get_atom(other_bond.a1_id)
                else:
                    continue
                neighbor_points.append(QPointF(other.x, other.y))
            if neighbor_points:
                avg_x = sum(p.x() for p in neighbor_points) / len(neighbor_points)
                avg_y = sum(p.y() for p in neighbor_points) / len(neighbor_points)
                vx = avg_x - mid.x()
                vy = avg_y - mid.y()
                return 1 if (nx * vx + ny * vy) >= 0 else -1

        def neighbor_score(atom_id: int) -> float:
            """Método auxiliar para neighbor score.

            Args:
                atom_id: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
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
        if is_editing and property_name in ("sub", "sup") and cursor.hasSelection():
            end_pos = cursor.selectionEnd()
            cursor.setPosition(end_pos)
            cursor.setCharFormat(fmt)
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

    def update_text_alignment(self, alignment: Qt.AlignmentFlag) -> None:
        """Actualiza text alignment.

        Args:
            alignment: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        items_to_update = self.scene.selectedItems()
        focus_item = self.scene.focusItem()
        if isinstance(focus_item, TextAnnotationItem) and focus_item not in items_to_update:
            items_to_update.append(focus_item)

        for item in items_to_update:
            if not isinstance(item, TextAnnotationItem):
                continue
            # Ensure there's a text width so alignment has visible effect.
            if alignment == Qt.AlignmentFlag.AlignLeft:
                item.setTextWidth(-1)
            else:
                rect = item.boundingRect()
                min_width = max(60.0, rect.width() + 20.0)
                item.setTextWidth(min_width)
            cursor = item.textCursor()
            is_editing = (
                item.textInteractionFlags()
                & Qt.TextInteractionFlag.TextEditorInteraction
            )
            if not is_editing:
                cursor.select(cursor.SelectionType.Document)
            block_format = cursor.blockFormat()
            block_format.setAlignment(alignment)
            cursor.setBlockFormat(block_format)
            item.setTextCursor(cursor)
            doc = item.document()
            option = doc.defaultTextOption()
            option.setAlignment(alignment)
            doc.setDefaultTextOption(option)

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
        use_circles = self.state.use_aromatic_circles
        rings = self._detect_aromatic_rings(with_atoms=True) if use_circles else []
        ring_pairs_with_circle: set[frozenset[int]] = set()
        for ring in rings:
            ring_pairs_with_circle.update(ring.get("bond_pairs", set()))

        # Update aromatic bond items:
        # - Render as single-edge only when a visible aromatic circle exists
        #   for that ring bond.
        # - Fall back to normal bond rendering otherwise.
        for bond_id, item in self.bond_items.items():
            if item.is_aromatic:
                bond = self.model.bonds.get(bond_id)
                render_as_circle = (
                    use_circles
                    and bond is not None
                    and frozenset({bond.a1_id, bond.a2_id}) in ring_pairs_with_circle
                )
                item.set_render_aromatic_as_circle(render_as_circle)
                self.update_bond_item(bond_id)

        # Remove existing circles
        for circle in self.aromatic_circles:
            self.scene.removeItem(circle)
        self.aromatic_circles.clear()
        
        if not use_circles:
            return
        
        # Draw circles for complete aromatic rings.
        for ring in rings:
            atom_ids = sorted(ring.get("atoms", set()))
            self._draw_aromatic_circle(atom_ids)

    def _trackball_affine_context(
        self,
    ) -> Optional[tuple[float, float, float, float, float, float]]:
        """Devuelve el contexto afín del trackball si la referencia sigue válida."""
        atom_ids = self._rotation_3d_ref_atom_ids
        if not atom_ids or not self._rotation_3d_ref_positions:
            return None
        if any(atom_id not in self.model.atoms for atom_id in atom_ids):
            return None
        if any(atom_id not in self._rotation_3d_ref_positions for atom_id in atom_ids):
            return None

        if not self._is_rotating_3d:
            projected = self._project_trackball_reference(
                atom_ids,
                self._rotation_3d_pitch_deg,
                self._rotation_3d_yaw_deg,
            )
            max_delta = 0.0
            for atom_id in atom_ids:
                px, py = projected[atom_id]
                atom = self.model.get_atom(atom_id)
                max_delta = max(max_delta, math.hypot(px - atom.x, py - atom.y))
            if max_delta > TRACKBALL_REFERENCE_MATCH_TOLERANCE_PX:
                return None

        count = float(len(atom_ids))
        center_x = sum(self._rotation_3d_ref_positions[atom_id][0] for atom_id in atom_ids) / count
        center_y = sum(self._rotation_3d_ref_positions[atom_id][1] for atom_id in atom_ids) / count
        pitch_rad = math.radians(self._rotation_3d_pitch_deg)
        yaw_rad = math.radians(self._rotation_3d_yaw_deg)
        return (
            center_x,
            center_y,
            math.cos(pitch_rad),
            math.sin(pitch_rad),
            math.cos(yaw_rad),
            math.sin(yaw_rad),
        )

    def _aromatic_circle_geometry_from_trackball_reference(
        self,
        atom_ids: Iterable[int],
        affine_ctx: tuple[float, float, float, float, float, float],
    ) -> Optional[tuple[float, float, float, float, float]]:
        """Calcula la elipse aromática aplicando la misma afín pseudo-3D."""
        ref_points: list[tuple[float, float]] = []
        for atom_id in atom_ids:
            point = self._rotation_3d_ref_positions.get(int(atom_id))
            if point is None:
                return None
            ref_points.append(point)
        if len(ref_points) < 3:
            return None

        ring_cx = sum(x for x, _ in ref_points) / len(ref_points)
        ring_cy = sum(y for _, y in ref_points) / len(ref_points)
        avg_dist = (
            sum(math.hypot(x - ring_cx, y - ring_cy) for x, y in ref_points)
            / len(ref_points)
        )
        base_radius = max(2.0, avg_dist * 0.65)

        center_x, center_y, cos_x, sin_x, cos_y, sin_y = affine_ctx
        x0 = ring_cx - center_x
        y0 = ring_cy - center_y
        y_new = y0 * cos_x
        z_new = y0 * sin_x
        proj_cx = x0 * cos_y + z_new * sin_y + center_x
        proj_cy = y_new + center_y

        # Imagen afín del círculo base de radio r en el plano XY.
        v1x = base_radius * cos_y
        v1y = 0.0
        v2x = base_radius * sin_x * sin_y
        v2y = base_radius * cos_x

        c11 = v1x * v1x + v2x * v2x
        c22 = v1y * v1y + v2y * v2y
        c12 = v1x * v1y + v2x * v2y
        trace = c11 + c22
        disc = max(0.0, (c11 - c22) * (c11 - c22) + 4.0 * c12 * c12)
        root = math.sqrt(disc)
        eig_1 = max(0.0, 0.5 * (trace + root))
        eig_2 = max(0.0, 0.5 * (trace - root))

        radius_x = max(2.0, math.sqrt(eig_1))
        radius_y = max(2.0, math.sqrt(eig_2))
        angle_deg = 0.5 * math.degrees(math.atan2(2.0 * c12, c11 - c22))
        return proj_cx, proj_cy, radius_x, radius_y, angle_deg

    def _aromatic_circle_geometry_for_atom_ids(
        self,
        atom_ids: Iterable[int],
        trackball_affine: Optional[tuple[float, float, float, float, float, float]] = None,
    ) -> Optional[tuple[float, float, float, float, float]]:
        """Calcula la geometría del marcador aromático."""
        atom_id_list = [int(atom_id) for atom_id in atom_ids]
        if len(atom_id_list) < 3:
            return None

        ring_atoms = []
        for atom_id in atom_id_list:
            atom = self.model.atoms.get(atom_id)
            if atom is not None:
                ring_atoms.append(atom)
        if len(ring_atoms) < 3:
            return None

        cx = sum(atom.x for atom in ring_atoms) / len(ring_atoms)
        cy = sum(atom.y for atom in ring_atoms) / len(ring_atoms)
        avg_dist = (
            sum(math.hypot(atom.x - cx, atom.y - cy) for atom in ring_atoms)
            / len(ring_atoms)
        )
        circular_radius = max(2.0, avg_dist * 0.65)

        centered = [(atom.x - cx, atom.y - cy) for atom in ring_atoms]
        c11 = sum(dx * dx for dx, _ in centered) / len(centered)
        c22 = sum(dy * dy for _, dy in centered) / len(centered)
        c12 = sum(dx * dy for dx, dy in centered) / len(centered)
        trace = c11 + c22
        disc = max(0.0, (c11 - c22) * (c11 - c22) + 4.0 * c12 * c12)
        root = math.sqrt(disc)
        eig_1 = max(0.0, 0.5 * (trace + root))
        eig_2 = max(0.0, 0.5 * (trace - root))

        max_eig = max(eig_1, eig_2)
        if max_eig <= 1e-9:
            return cx, cy, circular_radius, circular_radius, 0.0

        # Evita jitter angular en anillos prácticamente isotrópicos.
        if abs(eig_1 - eig_2) <= max(1e-6, max_eig * 1e-4):
            return cx, cy, circular_radius, circular_radius, 0.0

        # Para puntos distribuidos sobre una elipse, var ~= r^2 / 2.
        scale = 0.65 * math.sqrt(2.0)
        radius_x = max(2.0, math.sqrt(eig_1) * scale)
        radius_y = max(2.0, math.sqrt(eig_2) * scale)
        angle_deg = 0.5 * math.degrees(math.atan2(2.0 * c12, c11 - c22))
        return cx, cy, radius_x, radius_y, angle_deg

    def _draw_aromatic_circle(
        self,
        atom_ids: Iterable[int],
        trackball_affine: Optional[tuple[float, float, float, float, float, float]] = None,
    ) -> None:
        """Crea y añade a escena el círculo interno de un anillo aromático."""
        atom_ids_tuple = tuple(sorted(int(atom_id) for atom_id in atom_ids))
        if not atom_ids_tuple:
            return
        geometry = self._aromatic_circle_geometry_for_atom_ids(
            atom_ids_tuple,
            trackball_affine=trackball_affine,
        )
        if geometry is None:
            return
        cx, cy, radius_x, radius_y, angle_deg = geometry
        circle = AromaticCircleItem(QRectF(-radius_x, -radius_y, 2.0 * radius_x, 2.0 * radius_y))
        circle.set_geometry(cx, cy, radius_x, radius_y, angle_deg)
        circle.setData(AROMATIC_CIRCLE_ATOMS_ROLE, atom_ids_tuple)
        self.scene.addItem(circle)
        self.aromatic_circles.append(circle)

    def _update_aromatic_circles_for_atoms(self, atom_ids: set[int]) -> None:
        """Actualiza posición/tamaño de círculos aromáticos afectados por átomos movidos."""
        if not self.state.use_aromatic_circles or not self.aromatic_circles:
            return
        if not atom_ids:
            return

        changed_atoms = set(atom_ids)
        force_refresh = False
        for circle in list(self.aromatic_circles):
            try:
                ring_atom_ids = circle.data(AROMATIC_CIRCLE_ATOMS_ROLE)
            except RuntimeError:
                force_refresh = True
                break
            if not ring_atom_ids:
                force_refresh = True
                break
            ring_set = {int(atom_id) for atom_id in ring_atom_ids}
            if changed_atoms.isdisjoint(ring_set):
                continue
            geometry = self._aromatic_circle_geometry_for_atom_ids(
                ring_set,
            )
            if geometry is None:
                force_refresh = True
                break
            cx, cy, radius_x, radius_y, angle_deg = geometry
            circle.set_geometry(cx, cy, radius_x, radius_y, angle_deg)

        if force_refresh:
            self.refresh_aromatic_circles()

    def _detect_aromatic_rings(self, with_atoms: bool = False) -> list:
        """
        Detect complete aromatic rings and return their centers.
        Returns list of dicts with center and radius (and optionally atoms).
        """
        view = MolView(self.model)
        rings = find_rings_simple(view)
        if not rings:
            return []

        bond_by_pair = {
            frozenset({bond.a1_id, bond.a2_id}): bond
            for bond in self.model.bonds.values()
        }
        aromatic_pairs = {
            pair
            for pair, bond in bond_by_pair.items()
            if bond.is_aromatic
        }
        if not aromatic_pairs:
            return []
        
        # Calculate centers and radii for each ring
        result = []
        for ring in rings:
            ring_pairs = ring_bonds(view, ring)
            if not ring_pairs:
                continue
            # Modo estricto: todos los enlaces del ciclo son aromáticos.
            missing_pairs = [pair for pair in ring_pairs if pair not in aromatic_pairs]
            if missing_pairs:
                # Tolerancia visual: permitir hasta dos enlaces no aromáticos
                # cuando fueron estilizados para proyección 2D
                # (cuña/hash/bold). Esto evita perder el círculo aromático por
                # un ajuste de estilo local en anillos aromáticos.
                if len(missing_pairs) > 2:
                    continue
                missing_bonds = [bond_by_pair.get(pair) for pair in missing_pairs]
                if any(bond is None for bond in missing_bonds):
                    continue
                if any(
                    bond.style not in {BondStyle.WEDGE, BondStyle.HASHED, BondStyle.BOLD}
                    for bond in missing_bonds
                ):
                    continue
                aromatic_count = len(ring_pairs) - len(missing_pairs)
                if aromatic_count < max(1, len(ring_pairs) - 2):
                    continue
            if not ring_pairs:
                continue

            atom_ids = list(ring)
            atoms = [self.model.get_atom(aid) for aid in atom_ids if aid in self.model.atoms]
            if len(atoms) < 3:
                continue
            cx = sum(a.x for a in atoms) / len(atoms)
            cy = sum(a.y for a in atoms) / len(atoms)
            avg_dist = sum(((a.x - cx)**2 + (a.y - cy)**2)**0.5 for a in atoms) / len(atoms)
            radius = max(2.0, avg_dist * 0.65)
            entry = {
                "center": QPointF(cx, cy),
                "radius": radius,
            }
            if with_atoms:
                entry["atoms"] = set(atom_ids)
                entry["bond_pairs"] = set(ring_pairs)
            result.append(entry)
        
        return result

    def _aromatic_ring_center_for_bond(self, bond: Bond) -> Optional[QPointF]:
        """Método auxiliar para  aromatic ring center for bond.

        Args:
            bond: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if not bond.is_aromatic:
            return None
        adjacency: dict[int, list[tuple[int, int]]] = {}
        for b in self.model.bonds.values():
            if not b.is_aromatic:
                continue
            adjacency.setdefault(b.a1_id, []).append((b.a2_id, b.id))
            adjacency.setdefault(b.a2_id, []).append((b.a1_id, b.id))
        cycle = self._find_cycle_for_bond(bond, adjacency)
        if not cycle:
            return None
        return self._center_for_atoms(cycle["atom_ids"])

    def _refresh_aromatic_ring_contexts(self) -> None:
        """Método auxiliar para  refresh aromatic ring contexts.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        adjacency: dict[int, list[tuple[int, int]]] = {}
        for bond in self.model.bonds.values():
            if not bond.is_aromatic:
                continue
            adjacency.setdefault(bond.a1_id, []).append((bond.a2_id, bond.id))
            adjacency.setdefault(bond.a2_id, []).append((bond.a1_id, bond.id))
        if not adjacency:
            for bond in self.model.bonds.values():
                if bond.is_aromatic and bond.ring_id is None:
                    item = self.bond_items.get(bond.id)
                    if item is not None:
                        item.set_ring_context(None)
                        atom1 = self.model.get_atom(bond.a1_id)
                        atom2 = self.model.get_atom(bond.a2_id)
                        if atom1 and atom2:
                            item.update_positions(atom1, atom2)
            return
        for bond in self.model.bonds.values():
            if not bond.is_aromatic or bond.ring_id is not None:
                continue
            item = self.bond_items.get(bond.id)
            if item is None:
                continue
            cycle = self._find_cycle_for_bond(bond, adjacency)
            center = self._center_for_atoms(cycle["atom_ids"]) if cycle else None
            item.set_ring_context(center)
            atom1 = self.model.get_atom(bond.a1_id)
            atom2 = self.model.get_atom(bond.a2_id)
            if atom1 and atom2:
                item.update_positions(atom1, atom2)

    def _build_aromatic_adjacency(self) -> dict[int, list[int]]:
        """Método auxiliar para  build aromatic adjacency.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        adjacency: dict[int, list[int]] = {}
        for bond in self.model.bonds.values():
            if not bond.is_aromatic:
                continue
            adjacency.setdefault(bond.a1_id, []).append(bond.a2_id)
            adjacency.setdefault(bond.a2_id, []).append(bond.a1_id)
        return adjacency

    def _center_for_atoms(self, atom_ids: list[int]) -> Optional[QPointF]:
        """Método auxiliar para  center for atoms.

        Args:
            atom_ids: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if not atom_ids:
            return None
        xs = []
        ys = []
        for aid in atom_ids:
            atom = self.model.get_atom(aid)
            if atom is None:
                continue
            xs.append(atom.x)
            ys.append(atom.y)
        if not xs:
            return None
        return QPointF(sum(xs) / len(xs), sum(ys) / len(ys))

    def clean_2d_fallback(self, atom_ids: Optional[set[int]] = None, iterations: int = 50) -> None:
        """Basic 2D cleanup without RDKit using length + angle relaxation."""
        if atom_ids is None:
            atom_ids = set(self.model.atoms.keys())
        if not atom_ids:
            return

        before = {
            aid: (self.model.get_atom(aid).x, self.model.get_atom(aid).y) for aid in atom_ids
        }
        before_center = self._center_for_atoms(list(atom_ids)) or QPointF(0.0, 0.0)

        positions = {aid: QPointF(x, y) for aid, (x, y) in before.items()}

        target_len = float(self.state.bond_length)
        angle_step_deg = 12.0
        iterations = max(1, int(iterations))

        neighbor_map: dict[int, list[int]] = {}
        bonds_in_selection = []
        for bond in self.model.bonds.values():
            if bond.a1_id in atom_ids and bond.a2_id in atom_ids:
                neighbor_map.setdefault(bond.a1_id, []).append(bond.a2_id)
                neighbor_map.setdefault(bond.a2_id, []).append(bond.a1_id)
                bonds_in_selection.append(bond)

        components: list[set[int]] = []
        visited: set[int] = set()
        for aid in atom_ids:
            if aid in visited:
                continue
            stack = [aid]
            comp = set()
            while stack:
                node = stack.pop()
                if node in visited:
                    continue
                visited.add(node)
                comp.add(node)
                for nb in neighbor_map.get(node, []):
                    if nb not in visited:
                        stack.append(nb)
            components.append(comp)

        def rotate_point(p: QPointF, center: QPointF, delta_deg: float) -> QPointF:
            """Método auxiliar para rotate point.

            Args:
                p: Descripción del parámetro.
                center: Descripción del parámetro.
                delta_deg: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
            rad = math.radians(delta_deg)
            dx = p.x() - center.x()
            dy = p.y() - center.y()
            cos_t = math.cos(rad)
            sin_t = math.sin(rad)
            return QPointF(center.x() + dx * cos_t - dy * sin_t, center.y() + dx * sin_t + dy * cos_t)

        def normalize_tree(component: set[int], roots: list[int]) -> None:
            """Método auxiliar para normalize tree.

            Args:
                component: Descripción del parámetro.
                roots: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
            placed = set(roots)
            queue = list(roots)
            while queue:
                parent = queue.pop(0)
                p_parent = positions[parent]
                for neigh in neighbor_map.get(parent, []):
                    if neigh not in component or neigh in placed:
                        continue
                    ox, oy = before[neigh]
                    px, py = before[parent]
                    dx = ox - px
                    dy = oy - py
                    dist = math.hypot(dx, dy)
                    if dist < 1e-3:
                        dx, dy = target_len, 0.0
                        dist = target_len
                    nx = dx / dist
                    ny = dy / dist
                    positions[neigh] = QPointF(p_parent.x() + nx * target_len, p_parent.y() + ny * target_len)
                    placed.add(neigh)
                    queue.append(neigh)

        for component in components:
            if len(component) <= 1:
                continue
            # Find ring core by peeling leaves
            degrees = {aid: len([n for n in neighbor_map.get(aid, []) if n in component]) for aid in component}
            ring_atoms = set(component)
            leaf_queue = [aid for aid, deg in degrees.items() if deg <= 1]
            while leaf_queue:
                node = leaf_queue.pop()
                if node not in ring_atoms:
                    continue
                ring_atoms.remove(node)
                for nb in neighbor_map.get(node, []):
                    if nb in ring_atoms:
                        degrees[nb] -= 1
                        if degrees[nb] <= 1:
                            leaf_queue.append(nb)

            if not ring_atoms:
                # Acyclic: only normalize bond lengths, preserve angles
                root = next(iter(component))
                normalize_tree(component, [root])
                continue

            # Relax ring core
            for _ in range(iterations):
                for bond in bonds_in_selection:
                    if bond.a1_id not in ring_atoms or bond.a2_id not in ring_atoms:
                        continue
                    p1 = positions[bond.a1_id]
                    p2 = positions[bond.a2_id]
                    dx = p2.x() - p1.x()
                    dy = p2.y() - p1.y()
                    dist = math.hypot(dx, dy) or 1e-3
                    delta = (dist - target_len) / dist
                    p1 = QPointF(p1.x() + dx * 0.5 * delta, p1.y() + dy * 0.5 * delta)
                    p2 = QPointF(p2.x() - dx * 0.5 * delta, p2.y() - dy * 0.5 * delta)
                    positions[bond.a1_id] = p1
                    positions[bond.a2_id] = p2

                for aid in ring_atoms:
                    neighbors = [n for n in neighbor_map.get(aid, []) if n in ring_atoms]
                    if len(neighbors) != 2:
                        continue
                    n1, n2 = neighbors
                    center = positions[aid]
                    p1 = positions[n1]
                    p2 = positions[n2]
                    a1 = angle_deg(center, p1)
                    a2 = angle_deg(center, p2)
                    current = angle_distance_deg(a1, a2)

                    has_triple = False
                    has_double = False
                    for bond in bonds_in_selection:
                        if bond.a1_id != aid and bond.a2_id != aid:
                            continue
                        if bond.order >= 3:
                            has_triple = True
                        if bond.order == 2 or bond.is_aromatic:
                            has_double = True
                    if has_triple:
                        target = 180.0
                    elif has_double:
                        target = 120.0
                    else:
                        target = SP3_BOND_ANGLE_DEG

                    delta = target - current
                    if abs(delta) < 1e-2:
                        continue
                    delta = max(-angle_step_deg, min(angle_step_deg, delta))
                    p1_new = rotate_point(p1, center, delta * 0.5)
                    p2_new = rotate_point(p2, center, -delta * 0.5)
                    positions[n1] = p1_new
                    positions[n2] = p2_new

            # Normalize trees attached to ring core
            normalize_tree(component, list(ring_atoms))

        after_center = self._center_for_atoms(list(atom_ids)) or QPointF(0.0, 0.0)
        # Use computed positions for center
        xs = [positions[aid].x() for aid in atom_ids]
        ys = [positions[aid].y() for aid in atom_ids]
        if xs and ys:
            after_center = QPointF(sum(xs) / len(xs), sum(ys) / len(ys))
        shift = QPointF(before_center.x() - after_center.x(), before_center.y() - after_center.y())

        after = {}
        for aid in atom_ids:
            p = positions[aid]
            after[aid] = (p.x() + shift.x(), p.y() + shift.y())

        from gui.commands import MoveAtomsCommand
        cmd = MoveAtomsCommand(self.model, self, before, after)
        self.undo_stack.push(cmd)
        self._update_selection_overlay()

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
        """Método auxiliar para  get atom at.

        Args:
            x: Descripción del parámetro.
            y: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        for atom in self.model.atoms.values():
            dx = atom.x - x
            dy = atom.y - y
            if (dx * dx + dy * dy) < (ATOM_HIT_RADIUS * ATOM_HIT_RADIUS):
                return atom.id
        return None

    def _get_bond_at(self, scene_pos: QPointF) -> Optional[int]:
        """Método auxiliar para  get bond at.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        for item in self.scene.items(scene_pos):
            if isinstance(item, BondItem):
                return item.bond_id
        return None



    def _sync_selection_from_scene(self) -> None:
        """Método auxiliar para  sync selection from scene.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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

        for atom_id, item in list(self.atom_items.items()):
            try:
                item.set_selected(atom_id in selected_atoms or atom_id == self.bond_anchor_id)
            except RuntimeError:
                self.atom_items.pop(atom_id, None)

        self._update_selection_overlay()

        # Emit selection signal
        details = {}
        if len(selected_atoms) == 1 and not selected_bonds and not selected_text_items:
            atom = self.model.get_atom(next(iter(selected_atoms)))
            details = {
                "type": "atom",
                "id": atom.id,
                "element": atom.element,
                "charge": int(getattr(atom, "formal_charge", atom.charge)),
                "x": atom.x,
                "y": atom.y,
            }
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
        """Método auxiliar para  begin selection drag.

        Args:
            start_pos: Descripción del parámetro.
            free_select: Descripción del parámetro.
            additive: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  update selection drag.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._select_drag_mode == "free" and self._select_path is not None:
            self._select_path.lineTo(scene_pos)
            if self._select_preview_path is not None:
                self._select_preview_path.setPath(self._select_path)
            return
        if self._select_drag_mode == "rect" and self._select_start_pos is not None:
            if self._select_preview_rect is not None:
                self._select_preview_rect.setRect(QRectF(self._select_start_pos, scene_pos).normalized())

    def _finalize_selection_drag(self) -> None:
        """Método auxiliar para  finalize selection drag.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
            candidates = [
                item
                for item in self.scene.items(rect)
                if isinstance(
                    item,
                    (AtomItem, BondItem, ArrowItem, BracketItem, TextAnnotationItem),
                )
            ]
            items = [
                item
                for item in candidates
                if rect.contains(item.sceneBoundingRect().center())
            ]

        if not self._select_additive:
            self.scene.clearSelection()
        for item in items:
            item.setSelected(True)
        self._sync_selection_from_scene()
        self._clear_selection_drag()

    def _clear_selection_drag(self) -> None:
        """Método auxiliar para  clear selection drag.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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

    def _structure_bbox(self) -> Optional[QRectF]:
        """Método auxiliar para  structure bbox.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        rect: Optional[QRectF] = None

        def extend(candidate: QRectF) -> None:
            """Método auxiliar para extend.

            Args:
                candidate: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
            nonlocal rect
            if not candidate.isValid() or candidate.isNull():
                return
            rect = candidate if rect is None else rect.united(candidate)

        def extend_atom_bounds(atom_id: int) -> None:
            """Método auxiliar para extend atom bounds.

            Args:
                atom_id: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
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

        for atom_id in self.model.atoms.keys():
            extend_atom_bounds(atom_id)
        for bond_id in self.model.bonds.keys():
            item = self.bond_items.get(bond_id)
            if item is None:
                continue
            extend(item.sceneBoundingRect())
        return rect

    def _analysis_graph_and_bbox(self) -> tuple[Optional[MolGraph], Optional[QRectF]]:
        """Método auxiliar para  analysis graph and bbox.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atom_ids, bonds = self._selected_structure_ids()
        if atom_ids:
            return self._build_selection_graph(atom_ids, bonds), self._selected_items_bbox()
        if self.model.atoms:
            return self.model, self._structure_bbox()
        return None, None

    def _implicit_h_for_graph(self, graph: MolGraph, atom_id: int, element: str) -> int:
        """Método auxiliar para  implicit h for graph.

        Args:
            graph: Descripción del parámetro.
            atom_id: Descripción del parámetro.
            element: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if element not in IMPLICIT_H_ELEMENTS:
            return 0
        try:
            return int(max(0, graph.implicit_h_count(atom_id)))
        except Exception:
            return 0

    def _analysis_atom_counts(self, graph: MolGraph) -> dict[str, int]:
        """Método auxiliar para  analysis atom counts.

        Args:
            graph: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        counts: dict[str, int] = {}
        for atom in graph.atoms.values():
            counts[atom.element] = counts.get(atom.element, 0) + 1
            if atom.explicit_h:
                counts["H"] = counts.get("H", 0) + int(atom.explicit_h)
        for atom in graph.atoms.values():
            if atom.element == "H":
                continue
            implicit_h = self._implicit_h_for_graph(graph, atom.id, atom.element)
            if implicit_h > 0:
                counts["H"] = counts.get("H", 0) + implicit_h
        return {element: count for element, count in counts.items() if count > 0}

    @staticmethod
    def _analysis_formula(counts: dict[str, int]) -> str:
        """Método auxiliar para  analysis formula.

        Args:
            counts: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if not counts:
            return ""
        order: list[str] = []
        if "C" in counts:
            order.append("C")
        if "H" in counts:
            order.append("H")
        for element in sorted(e for e in counts.keys() if e not in {"C", "H"}):
            order.append(element)
        parts = []
        for element in order:
            count = counts.get(element, 0)
            if count <= 0:
                continue
            parts.append(element if count == 1 else f"{element}{count}")
        return "".join(parts)

    @staticmethod
    def _analysis_exact_mass(counts: dict[str, int]) -> Optional[float]:
        """Método auxiliar para  analysis exact mass.

        Args:
            counts: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        total = 0.0
        for element, count in counts.items():
            mass = MONOISOTOPIC_MASSES.get(element)
            if mass is None:
                return None
            total += mass * count
        return total

    @staticmethod
    def _analysis_molecular_weight(counts: dict[str, int]) -> Optional[float]:
        """Método auxiliar para  analysis molecular weight.

        Args:
            counts: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        total = 0.0
        for element, count in counts.items():
            mass = ATOMIC_WEIGHTS.get(element)
            if mass is None:
                return None
            total += mass * count
        return total

    @staticmethod
    def _analysis_convolve(
        distribution: list[tuple[float, float]],
        isotopes: list[tuple[float, float]],
    ) -> list[tuple[float, float]]:
        """Método auxiliar para  analysis convolve.

        Args:
            distribution: Descripción del parámetro.
            isotopes: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if not distribution:
            return []
        new_dist: dict[float, float] = {}
        for mass, prob in distribution:
            for iso_mass, iso_prob in isotopes:
                combined_mass = round(mass + iso_mass, 4)
                new_dist[combined_mass] = new_dist.get(combined_mass, 0.0) + prob * iso_prob
        if not new_dist:
            return []
        max_prob = max(new_dist.values())
        min_prob = max_prob * 1e-8
        pruned = {m: p for m, p in new_dist.items() if p >= min_prob}
        if len(pruned) > ANALYSIS_DIST_KEEP:
            items = sorted(pruned.items(), key=lambda item: item[1], reverse=True)[:ANALYSIS_DIST_KEEP]
            pruned = dict(items)
        return list(pruned.items())

    def _analysis_isotope_peaks(self, counts: dict[str, int]) -> Optional[list[tuple[float, float]]]:
        """Método auxiliar para  analysis isotope peaks.

        Args:
            counts: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        distribution: list[tuple[float, float]] = [(0.0, 1.0)]
        for element in sorted(counts.keys()):
            isotopes = ISOTOPE_ABUNDANCES.get(element)
            if isotopes is None:
                return None
            for _ in range(counts[element]):
                distribution = self._analysis_convolve(distribution, isotopes)
        if not distribution:
            return None
        binned: dict[float, float] = {}
        for mass, prob in distribution:
            mass_key = round(mass, 2)
            binned[mass_key] = binned.get(mass_key, 0.0) + prob
        max_prob = max(binned.values())
        if max_prob <= 0:
            return None
        peaks = [(mass, (prob / max_prob) * 100.0) for mass, prob in binned.items()]
        peaks = [peak for peak in peaks if peak[1] >= ANALYSIS_MIN_PEAK_PERCENT]
        peaks.sort(key=lambda item: item[1], reverse=True)
        if len(peaks) > ANALYSIS_MAX_PEAKS:
            peaks = peaks[:ANALYSIS_MAX_PEAKS]
        return peaks

    def _analysis_elemental_line(self, counts: dict[str, int], molecular_weight: Optional[float]) -> Optional[str]:
        """Método auxiliar para  analysis elemental line.

        Args:
            counts: Descripción del parámetro.
            molecular_weight: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if molecular_weight is None or molecular_weight <= 0:
            return None
        order: list[str] = []
        if "C" in counts:
            order.append("C")
        if "H" in counts:
            order.append("H")
        for element in sorted(e for e in counts.keys() if e not in {"C", "H"}):
            order.append(element)
        parts = []
        for element in order:
            weight = ATOMIC_WEIGHTS.get(element)
            if weight is None:
                continue
            percent = (weight * counts[element] / molecular_weight) * 100.0
            parts.append(f"{element}, {percent:.2f}")
        if not parts:
            return None
        return "Elemental Analysis: " + "; ".join(parts)

    def _analysis_build_text(self, graph: MolGraph, mode: str) -> str:
        """Método auxiliar para  analysis build text.

        Args:
            graph: Descripción del parámetro.
            mode: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        counts = self._analysis_atom_counts(graph)
        if not counts:
            return ""
        formula = self._analysis_formula(counts)
        exact_mass = self._analysis_exact_mass(counts)
        molecular_weight = self._analysis_molecular_weight(counts)
        peaks = self._analysis_isotope_peaks(counts)

        lines: list[str] = []
        if mode in {"name", "all", "iupac"}:
            try:
                iupac = iupac_name(graph, NameOptions())
            except Exception:
                iupac = "N/D"
            if mode in {"name", "all"}:
                lines.append(f"IUPAC: {iupac}")
            else:
                lines.append(iupac)
        if mode in {"name", "all"}:
            smiles = ""
            try:
                smiles = molgraph_to_smiles(graph)
            except Exception:
                smiles = ""
            lines.append(f"SMILES: {smiles or 'N/D'}")
        if mode in {"formula", "all"}:
            lines.append(f"Chemical Formula: {formula}")
        if mode in {"exact", "all"}:
            lines.append(
                f"Exact Mass: {exact_mass:.2f}" if exact_mass is not None else "Exact Mass: N/D"
            )
        if mode in {"weight", "all"}:
            lines.append(
                f"Molecular Weight: {molecular_weight:.2f}"
                if molecular_weight is not None
                else "Molecular Weight: N/D"
            )
        if mode in {"mz", "all"}:
            if peaks:
                formatted = ", ".join(f"{mass:.2f} ({percent:.1f}%)" for mass, percent in peaks)
                lines.append(f"m/z: {formatted}")
            else:
                lines.append("m/z: N/D")
        if mode in {"elemental", "all"}:
            elemental_line = self._analysis_elemental_line(counts, molecular_weight)
            lines.append(elemental_line or "Elemental Analysis: N/D")
        return "\n".join(lines)

    def _insert_analysis_text(
        self,
        text: str,
        bbox: Optional[QRectF],
        scene_pos: QPointF,
    ) -> None:
        """Método auxiliar para  insert analysis text.

        Args:
            text: Descripción del parámetro.
            bbox: Descripción del parámetro.
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = TextAnnotationItem(text, 0.0, 0.0)
        self._apply_text_settings(item)
        item.setTextInteractionFlags(Qt.TextInteractionFlag.NoTextInteraction)
        try:
            item.document().setDocumentMargin(0)
            item.document().setDefaultStyleSheet("body { background: transparent; }")
        except Exception:
            pass
        if bbox is not None and bbox.isValid():
            x = bbox.left()
            y = bbox.bottom() + ANALYSIS_MARGIN_PX
        else:
            x = scene_pos.x() + 10.0
            y = scene_pos.y() + 10.0
        item.setPos(x, y)
        self.scene.clearSelection()
        self.undo_stack.push(AddTextItemCommand(self, item))
        try:
            item.setSelected(True)
        except RuntimeError:
            pass

    def _run_analysis_action(self, mode: str, scene_pos: QPointF) -> None:
        """Método auxiliar para  run analysis action.

        Args:
            mode: Descripción del parámetro.
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        graph, bbox = self._analysis_graph_and_bbox()
        if graph is None:
            return
        text = self._analysis_build_text(graph, mode)
        if not text:
            return
        self._insert_analysis_text(text, bbox, scene_pos)

    def run_analysis(self, mode: str) -> None:
        """Método auxiliar para run analysis.

        Args:
            mode: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        scene_pos = self._last_scene_pos
        if scene_pos is None:
            try:
                center = self.viewport().rect().center()
                scene_pos = self.mapToScene(center)
            except Exception:
                scene_pos = QPointF(0.0, 0.0)
        self._run_analysis_action(mode, scene_pos)

    def _show_context_menu(
        self,
        scene_pos: QPointF,
        global_pos,
        clicked_item: Optional[QGraphicsItem] = None,
    ) -> None:
        """Método auxiliar para  show context menu.

        Args:
            scene_pos: Descripción del parámetro.
            global_pos: Descripción del parámetro.
            clicked_item: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        menu = QMenu(self)
        has_selection = bool(
            self.state.selected_atoms
            or self.state.selected_bonds
            or self._selected_text_items()
            or any(isinstance(item, (ArrowItem, BracketItem)) for item in self.scene.selectedItems())
        )
        has_bond_selection = bool(self.state.selected_bonds)
        coordination_bond_ids: list[int] = []
        if isinstance(clicked_item, BondItem):
            coordination_bond_ids = [clicked_item.bond_id]
        elif has_bond_selection:
            coordination_bond_ids = sorted(
                bond_id for bond_id in self.state.selected_bonds if bond_id in self.model.bonds
            )
        coordination_bond_enable = False
        if coordination_bond_ids:
            coordination_bond_enable = not all(
                self.model.get_bond(bond_id).style == BondStyle.COORDINATION
                for bond_id in coordination_bond_ids
            )
        coordination_atom_ids: list[int] = []
        if isinstance(clicked_item, AtomItem):
            coordination_atom_ids = [clicked_item.atom_id]
        else:
            coordination_atom_ids = sorted(
                atom_id
                for atom_id in self._selected_atom_ids_for_transform()
                if atom_id in self.model.atoms
            )
        coordination_enable = False
        if coordination_atom_ids:
            coordination_enable = not all(
                getattr(self.model.get_atom(atom_id), "is_coordination_center", False)
                for atom_id in coordination_atom_ids
            )
        styled_coordination_ids = [
            atom_id
            for atom_id in coordination_atom_ids
            if bool(getattr(self.model.get_atom(atom_id), "is_coordination_center", False))
        ]
        coordination_all_filled = bool(styled_coordination_ids) and all(
            bool(getattr(self.model.get_atom(atom_id), "sphere_filled", True))
            for atom_id in styled_coordination_ids
        )
        coordination_all_transparent = bool(styled_coordination_ids) and all(
            bool(getattr(self.model.get_atom(atom_id), "sphere_transparent", False))
            for atom_id in styled_coordination_ids
        )
        coordination_has_custom_color = any(
            bool(getattr(self.model.get_atom(atom_id), "sphere_color", None))
            for atom_id in styled_coordination_ids
        )
        charge_atom_ids: list[int] = []
        if isinstance(clicked_item, AtomItem):
            charge_atom_ids = [clicked_item.atom_id]
        elif self.state.selected_atoms:
            charge_atom_ids = sorted(
                atom_id for atom_id in self.state.selected_atoms if atom_id in self.model.atoms
            )
        charge_all_no_implicit = bool(charge_atom_ids) and all(
            bool(getattr(self.model.get_atom(atom_id), "no_implicit", False))
            for atom_id in charge_atom_ids
        )

        act_cut = menu.addAction("Cortar")
        act_copy = menu.addAction("Copiar")
        act_paste = menu.addAction("Pegar")
        menu.addSeparator()
        act_delete = menu.addAction("Eliminar")
        act_thicker = None
        act_thinner = None
        act_reset_thickness = None
        act_color = None
        act_reset_color = None
        if has_bond_selection:
            menu.addSeparator()
            thickness_menu = menu.addMenu("Grosor de enlace")
            act_thicker = thickness_menu.addAction("Incrementar grosor")
            act_thinner = thickness_menu.addAction("Disminuir grosor")
            act_reset_thickness = thickness_menu.addAction("Restablecer grosor")
            act_color = menu.addAction("Color de enlace...")
            act_reset_color = menu.addAction("Restablecer color de enlace")
        act_anchor = None
        act_charge_increment = None
        act_charge_decrement = None
        act_charge_set = None
        act_charge_reset = None
        act_no_implicit = None
        act_toggle_coord_sphere = None
        act_toggle_coord_bond = None
        act_coord_sphere_color = None
        act_coord_sphere_reset_color = None
        act_coord_sphere_toggle_fill = None
        act_coord_sphere_toggle_transparent = None
        act_coord_sphere_radius = None
        act_arrange_coordination = None
        if isinstance(clicked_item, AtomItem):
            atom = self.model.get_atom(clicked_item.atom_id)
            if atom is not None and atom.element not in ELEMENT_SYMBOLS:
                candidates = self._label_anchor_candidates(atom.element)
                if candidates:
                    menu.addSeparator()
                    act_anchor = menu.addAction("Elegir átomo de unión...")
        if charge_atom_ids:
            menu.addSeparator()
            charge_menu = menu.addMenu("Carga formal")
            act_charge_increment = charge_menu.addAction("Incrementar (+1)")
            act_charge_decrement = charge_menu.addAction("Disminuir (-1)")
            act_charge_set = charge_menu.addAction("Establecer...")
            act_charge_reset = charge_menu.addAction("Restablecer a 0")
            charge_menu.addSeparator()
            act_no_implicit = charge_menu.addAction("No implícitos")
            act_no_implicit.setCheckable(True)
            act_no_implicit.setChecked(charge_all_no_implicit)
        if coordination_bond_ids:
            menu.addSeparator()
            bond_text = (
                "Convertir a enlace coordinativo"
                if coordination_bond_enable
                else "Convertir a enlace normal"
            )
            act_toggle_coord_bond = menu.addAction(bond_text)
        if coordination_atom_ids:
            menu.addSeparator()
            coord_text = (
                "Activar esfera de coordinación"
                if coordination_enable
                else "Desactivar esfera de coordinación"
            )
            act_toggle_coord_sphere = menu.addAction(coord_text)
        if styled_coordination_ids:
            sphere_style_menu = menu.addMenu("Estilo de esfera")
            transparent_text = (
                "Mostrar esfera"
                if coordination_all_transparent
                else "Hacer esfera transparente"
            )
            act_coord_sphere_toggle_transparent = sphere_style_menu.addAction(transparent_text)
            act_coord_sphere_color = sphere_style_menu.addAction("Color de esfera...")
            act_coord_sphere_reset_color = sphere_style_menu.addAction("Restablecer color")
            fill_text = (
                "Quitar fondo de esfera"
                if coordination_all_filled
                else "Mostrar fondo de esfera"
            )
            act_coord_sphere_toggle_fill = sphere_style_menu.addAction(fill_text)
            act_coord_sphere_radius = sphere_style_menu.addAction("Tamaño de esfera...")
            act_coord_sphere_reset_color.setEnabled(coordination_has_custom_color)
            act_coord_sphere_color.setEnabled(not coordination_all_transparent)
            act_coord_sphere_reset_color.setEnabled(
                coordination_has_custom_color and not coordination_all_transparent
            )
            act_coord_sphere_toggle_fill.setEnabled(not coordination_all_transparent)
        arrangement_center_ids = [
            atom_id
            for atom_id in coordination_atom_ids
            if bool(getattr(self.model.get_atom(atom_id), "is_coordination_center", False))
        ]
        if arrangement_center_ids:
            can_arrange = any(
                len(self._coordination_ligands_for_center(center_id)) >= 2
                for center_id in arrangement_center_ids
            )
            act_arrange_coordination = menu.addAction("Distribuir ligandos alrededor")
            act_arrange_coordination.setEnabled(can_arrange)
        menu.addSeparator()
        act_select_all = menu.addAction("Seleccionar todo")
        menu.addSeparator()
        analysis_menu = menu.addMenu("Análisis")
        act_analysis_name = analysis_menu.addAction("Nombre (SMILES)")
        act_analysis_iupac = analysis_menu.addAction("Nombre (IUPAC)")
        act_analysis_formula = analysis_menu.addAction("Fórmula química")
        act_analysis_exact = analysis_menu.addAction("Masa exacta")
        act_analysis_weight = analysis_menu.addAction("Peso molecular")
        act_analysis_mz = analysis_menu.addAction("m/z")
        act_analysis_elemental = analysis_menu.addAction("Análisis elemental")
        analysis_menu.addSeparator()
        act_analysis_all = analysis_menu.addAction("Todo")

        act_cut.setEnabled(has_selection)
        act_copy.setEnabled(has_selection)
        act_delete.setEnabled(has_selection)
        analysis_menu.setEnabled(bool(self.model.atoms))

        action = menu.exec(global_pos)
        if action is None:
            return
        if act_thicker is not None and action == act_thicker:
            self._adjust_selected_bond_stroke(self._bond_stroke_step())
            return
        if act_thinner is not None and action == act_thinner:
            self._adjust_selected_bond_stroke(-self._bond_stroke_step())
            return
        if act_reset_thickness is not None and action == act_reset_thickness:
            self._reset_selected_bond_stroke()
            return
        if act_color is not None and action == act_color:
            self._prompt_bond_color()
            return
        if act_reset_color is not None and action == act_reset_color:
            self._set_selected_bond_color(None)
            return
        if act_anchor is not None and action == act_anchor and isinstance(clicked_item, AtomItem):
            self._prompt_anchor_for_atom(clicked_item.atom_id)
            return
        if act_charge_increment is not None and action == act_charge_increment:
            self._adjust_atom_formal_charge(charge_atom_ids, +1)
            return
        if act_charge_decrement is not None and action == act_charge_decrement:
            self._adjust_atom_formal_charge(charge_atom_ids, -1)
            return
        if act_charge_set is not None and action == act_charge_set:
            self._prompt_set_atom_formal_charge(charge_atom_ids)
            return
        if act_charge_reset is not None and action == act_charge_reset:
            self._set_atom_formal_charge(charge_atom_ids, 0)
            return
        if act_no_implicit is not None and action == act_no_implicit:
            self._set_atom_no_implicit(charge_atom_ids, not charge_all_no_implicit)
            return
        if act_toggle_coord_bond is not None and action == act_toggle_coord_bond:
            self._set_bond_coordination_style(coordination_bond_ids, coordination_bond_enable)
            return
        if act_toggle_coord_sphere is not None and action == act_toggle_coord_sphere:
            self._set_coordination_sphere(coordination_atom_ids, coordination_enable)
            return
        if (
            act_coord_sphere_toggle_transparent is not None
            and action == act_coord_sphere_toggle_transparent
        ):
            self._set_coordination_sphere_transparent(
                styled_coordination_ids,
                not coordination_all_transparent,
            )
            return
        if act_coord_sphere_color is not None and action == act_coord_sphere_color:
            self._prompt_coordination_sphere_color(styled_coordination_ids)
            return
        if act_coord_sphere_reset_color is not None and action == act_coord_sphere_reset_color:
            self._set_coordination_sphere_color(styled_coordination_ids, None)
            return
        if act_coord_sphere_toggle_fill is not None and action == act_coord_sphere_toggle_fill:
            self._set_coordination_sphere_fill(styled_coordination_ids, not coordination_all_filled)
            return
        if act_coord_sphere_radius is not None and action == act_coord_sphere_radius:
            self._prompt_coordination_sphere_radius(styled_coordination_ids)
            return
        if act_arrange_coordination is not None and action == act_arrange_coordination:
            self._auto_arrange_coordination_ligands(arrangement_center_ids)
            return
        if action == act_analysis_name:
            self._run_analysis_action("name", scene_pos)
            return
        if action == act_analysis_iupac:
            self._run_analysis_action("iupac", scene_pos)
            return
        if action == act_analysis_formula:
            self._run_analysis_action("formula", scene_pos)
            return
        if action == act_analysis_exact:
            self._run_analysis_action("exact", scene_pos)
            return
        if action == act_analysis_weight:
            self._run_analysis_action("weight", scene_pos)
            return
        if action == act_analysis_mz:
            self._run_analysis_action("mz", scene_pos)
            return
        if action == act_analysis_elemental:
            self._run_analysis_action("elemental", scene_pos)
            return
        if action == act_analysis_all:
            self._run_analysis_action("all", scene_pos)
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

    def _adjust_atom_formal_charge(self, atom_ids: Iterable[int], delta: int) -> None:
        """Ajusta la carga formal de los átomos seleccionados."""
        valid_ids = [atom_id for atom_id in atom_ids if atom_id in self.model.atoms]
        if not valid_ids:
            return
        self.undo_stack.beginMacro("Adjust atom formal charge")
        for atom_id in valid_ids:
            atom = self.model.get_atom(atom_id)
            self.undo_stack.push(
                ChangeChargeCommand(self.model, self, atom_id, int(atom.charge) + int(delta))
            )
        self.undo_stack.endMacro()

    def _set_atom_formal_charge(self, atom_ids: Iterable[int], charge: int) -> None:
        """Establece una carga formal fija para una lista de átomos."""
        valid_ids = [atom_id for atom_id in atom_ids if atom_id in self.model.atoms]
        if not valid_ids:
            return
        normalized_charge = int(charge)
        self.undo_stack.beginMacro("Set atom formal charge")
        for atom_id in valid_ids:
            self.undo_stack.push(ChangeChargeCommand(self.model, self, atom_id, normalized_charge))
        self.undo_stack.endMacro()

    def _prompt_set_atom_formal_charge(self, atom_ids: Iterable[int]) -> None:
        """Solicita la carga formal y la aplica a los átomos seleccionados."""
        valid_ids = [atom_id for atom_id in atom_ids if atom_id in self.model.atoms]
        if not valid_ids:
            return
        first_charge = int(self.model.get_atom(valid_ids[0]).charge)
        value, ok = QInputDialog.getInt(
            self,
            "Carga formal",
            "Carga formal:",
            first_charge,
            -15,
            15,
            1,
        )
        if not ok:
            return
        self._set_atom_formal_charge(valid_ids, value)

    def _set_atom_no_implicit(self, atom_ids: Iterable[int], enabled: bool) -> None:
        """Activa/desactiva la bandera `no_implicit` para átomos seleccionados."""
        valid_ids = [atom_id for atom_id in atom_ids if atom_id in self.model.atoms]
        if not valid_ids:
            return
        enabled_flag = bool(enabled)
        self.undo_stack.beginMacro("Set no implicit hydrogens")
        for atom_id in valid_ids:
            self.undo_stack.push(
                ChangeNoImplicitCommand(self.model, self, atom_id, enabled_flag)
            )
        self.undo_stack.endMacro()

    def _set_coordination_sphere(self, atom_ids: Iterable[int], enabled: bool) -> None:
        """Activa/desactiva esfera de coordinación para una lista de átomos."""
        enabled_flag = bool(enabled)
        valid_ids = [
            atom_id
            for atom_id in atom_ids
            if atom_id in self.model.atoms
            and bool(getattr(self.model.get_atom(atom_id), "is_coordination_center", False))
            != enabled_flag
        ]
        if not valid_ids:
            return
        self.undo_stack.beginMacro(
            "Enable coordination spheres" if enabled_flag else "Disable coordination spheres"
        )
        for atom_id in valid_ids:
            self.undo_stack.push(
                SetCoordinationCenterCommand(
                    self.model,
                    self,
                    atom_id,
                    enabled_flag,
                )
            )
        self.undo_stack.endMacro()

    def _prompt_coordination_sphere_color(self, atom_ids: Iterable[int]) -> None:
        """Solicita color de esfera y lo aplica a los átomos seleccionados."""
        valid_ids = [atom_id for atom_id in atom_ids if atom_id in self.model.atoms]
        if not valid_ids:
            return
        initial = QColor("#D9DDE3")
        for atom_id in valid_ids:
            atom = self.model.get_atom(atom_id)
            color = getattr(atom, "sphere_color", None)
            if color:
                candidate = QColor(color)
                if candidate.isValid():
                    initial = candidate
                    break
        color = QColorDialog.getColor(initial, self, "Seleccionar color de esfera")
        if not color.isValid():
            return
        self._set_coordination_sphere_color(valid_ids, color.name())

    def _set_coordination_sphere_color(
        self, atom_ids: Iterable[int], color: Optional[str]
    ) -> None:
        """Aplica color de esfera a una lista de átomos."""
        valid_ids = [atom_id for atom_id in atom_ids if atom_id in self.model.atoms]
        if not valid_ids:
            return
        self.undo_stack.beginMacro("Set coordination sphere color")
        for atom_id in valid_ids:
            self.undo_stack.push(
                ChangeCoordinationSphereStyleCommand(
                    self.model,
                    self,
                    atom_id,
                    new_color=color,
                )
            )
        self.undo_stack.endMacro()

    def _set_coordination_sphere_transparent(
        self, atom_ids: Iterable[int], transparent: bool
    ) -> None:
        """Activa/desactiva esfera transparente (solo etiqueta)."""
        valid_ids = [atom_id for atom_id in atom_ids if atom_id in self.model.atoms]
        if not valid_ids:
            return
        self.undo_stack.beginMacro("Toggle coordination sphere transparency")
        for atom_id in valid_ids:
            self.undo_stack.push(
                ChangeCoordinationSphereStyleCommand(
                    self.model,
                    self,
                    atom_id,
                    new_transparent=bool(transparent),
                )
            )
        self.undo_stack.endMacro()

    def _set_coordination_sphere_fill(self, atom_ids: Iterable[int], filled: bool) -> None:
        """Activa o desactiva el relleno de la esfera de coordinación."""
        valid_ids = [atom_id for atom_id in atom_ids if atom_id in self.model.atoms]
        if not valid_ids:
            return
        self.undo_stack.beginMacro("Toggle coordination sphere fill")
        for atom_id in valid_ids:
            self.undo_stack.push(
                ChangeCoordinationSphereStyleCommand(
                    self.model,
                    self,
                    atom_id,
                    new_filled=bool(filled),
                )
            )
        self.undo_stack.endMacro()

    def _prompt_coordination_sphere_radius(self, atom_ids: Iterable[int]) -> None:
        """Solicita un radio de esfera y lo aplica a los átomos seleccionados."""
        valid_ids = [atom_id for atom_id in atom_ids if atom_id in self.model.atoms]
        if not valid_ids:
            return
        current = 16.0
        for atom_id in valid_ids:
            atom = self.model.get_atom(atom_id)
            radius = getattr(atom, "sphere_radius", None)
            if radius is not None:
                current = float(radius)
                break
        value, ok = QInputDialog.getDouble(
            self,
            "Tamaño de esfera",
            "Radio (px):",
            current,
            4.0,
            200.0,
            1,
        )
        if not ok:
            return
        self.undo_stack.beginMacro("Set coordination sphere radius")
        for atom_id in valid_ids:
            self.undo_stack.push(
                ChangeCoordinationSphereStyleCommand(
                    self.model,
                    self,
                    atom_id,
                    new_radius=float(value),
                )
            )
        self.undo_stack.endMacro()

    def _set_bond_coordination_style(self, bond_ids: Iterable[int], enabled: bool) -> None:
        """Convierte enlaces seleccionados entre normal y coordinativo."""
        target_style = BondStyle.COORDINATION if enabled else BondStyle.PLAIN
        valid_ids = [bond_id for bond_id in bond_ids if bond_id in self.model.bonds]
        if not valid_ids:
            return
        changed_ids = [
            bond_id
            for bond_id in valid_ids
            if self.model.get_bond(bond_id).style != target_style
        ]
        if not changed_ids:
            return
        self.undo_stack.beginMacro(
            "Set coordination bond style" if enabled else "Set normal bond style"
        )
        for bond_id in changed_ids:
            bond = self.model.get_bond(bond_id)
            donor_atom_id = None
            if enabled:
                donor_atom_id = self._infer_coordination_donor_atom(
                    bond.a1_id,
                    bond.a2_id,
                    preferred=getattr(bond, "donor_atom_id", None),
                )
            self.undo_stack.push(
                ChangeBondCommand(
                    self.model,
                    self,
                    bond_id,
                    new_order=1,
                    new_style=target_style,
                    new_stereo=BondStereo.NONE,
                    new_is_aromatic=False,
                    new_donor_atom_id=donor_atom_id,
                )
            )
        self.undo_stack.endMacro()

    def _coordination_ligands_for_center(self, center_atom_id: int) -> list[int]:
        """Lista ligandos conectados al centro mediante enlaces coordinativos."""
        if center_atom_id not in self.model.atoms:
            return []
        ligands: list[int] = []
        seen: set[int] = set()
        for bond in self.model.bonds.values():
            if bond.style != BondStyle.COORDINATION:
                continue
            other_id: Optional[int] = None
            if bond.a1_id == center_atom_id:
                other_id = bond.a2_id
            elif bond.a2_id == center_atom_id:
                other_id = bond.a1_id
            if other_id is None or other_id == center_atom_id:
                continue
            if other_id not in self.model.atoms or other_id in seen:
                continue
            seen.add(other_id)
            ligands.append(other_id)
        return ligands

    @staticmethod
    def _coordination_arrangement_step_deg(count: int) -> float:
        """Devuelve separación angular para distribuir ligandos."""
        if count <= 1:
            return 360.0
        if count == 2:
            return 180.0
        if count == 3:
            return 120.0
        if count == 4:
            return 90.0
        if count == 6:
            return 60.0
        return 360.0 / float(count)

    def _coordination_arrangement_radius(self, center_id: int, ligand_ids: list[int]) -> float:
        """Calcula radio objetivo para auto-arrange de ligandos."""
        center = self.model.get_atom(center_id)
        distances: list[float] = []
        for ligand_id in ligand_ids:
            ligand = self.model.atoms.get(ligand_id)
            if ligand is None:
                continue
            distances.append(math.hypot(ligand.x - center.x, ligand.y - center.y))
        base_length = max(18.0, float(self.state.bond_length))
        if not distances:
            return base_length
        avg = sum(distances) / float(len(distances))
        target = max(base_length * 0.8, avg)
        return min(target, base_length * 2.2)

    def _auto_arrange_coordination_ligands(self, center_atom_ids: Iterable[int]) -> None:
        """Distribuye ligandos coordinativos alrededor de uno o más centros."""
        centers = [
            atom_id
            for atom_id in center_atom_ids
            if atom_id in self.model.atoms
            and bool(getattr(self.model.get_atom(atom_id), "is_coordination_center", False))
        ]
        if not centers:
            return
        before: Dict[int, Tuple[float, float]] = {}
        after: Dict[int, Tuple[float, float]] = {}
        moved_ligands: set[int] = set()
        for center_id in centers:
            center = self.model.get_atom(center_id)
            ligand_ids = self._coordination_ligands_for_center(center_id)
            if len(ligand_ids) < 2:
                continue
            ordered_ligands = sorted(
                ligand_ids,
                key=lambda atom_id: math.atan2(
                    self.model.get_atom(atom_id).y - center.y,
                    self.model.get_atom(atom_id).x - center.x,
                ),
            )
            step_deg = self._coordination_arrangement_step_deg(len(ordered_ligands))
            first_atom = self.model.get_atom(ordered_ligands[0])
            start_angle_deg = math.degrees(
                math.atan2(first_atom.y - center.y, first_atom.x - center.x)
            )
            radius = self._coordination_arrangement_radius(center_id, ordered_ligands)

            for index, ligand_id in enumerate(ordered_ligands):
                if ligand_id in moved_ligands:
                    continue
                ligand = self.model.get_atom(ligand_id)
                angle_deg = start_angle_deg + step_deg * index
                angle_rad = math.radians(angle_deg)
                new_x = center.x + radius * math.cos(angle_rad)
                new_y = center.y + radius * math.sin(angle_rad)
                if abs(new_x - ligand.x) < 0.05 and abs(new_y - ligand.y) < 0.05:
                    continue
                before[ligand_id] = (ligand.x, ligand.y)
                after[ligand_id] = (new_x, new_y)
                moved_ligands.add(ligand_id)
        if not after:
            return
        self.undo_stack.push(MoveAtomsCommand(self.model, self, before, after))

    def _bond_stroke_step(self) -> float:
        """Método auxiliar para  bond stroke step.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        return max(0.6, self.drawing_style.stroke_px * 0.35)

    def increase_selected_bond_thickness(self) -> None:
        """Método auxiliar para increase selected bond thickness.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._adjust_selected_bond_stroke(self._bond_stroke_step())

    def decrease_selected_bond_thickness(self) -> None:
        """Método auxiliar para decrease selected bond thickness.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._adjust_selected_bond_stroke(-self._bond_stroke_step())

    def reset_selected_bond_thickness(self) -> None:
        """Método auxiliar para reset selected bond thickness.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._reset_selected_bond_stroke()

    def _prompt_bond_color(self) -> None:
        """Método auxiliar para  prompt bond color.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        bond_ids = list(self.state.selected_bonds)
        if not bond_ids:
            return
        initial = QColor(self.drawing_style.bond_color)
        for bond_id in bond_ids:
            bond = self.model.get_bond(bond_id)
            if bond.color:
                initial = QColor(bond.color)
                break
        color = QColorDialog.getColor(initial, self, "Seleccionar color de enlace")
        if not color.isValid():
            return
        self._set_selected_bond_color(color.name())

    def _set_selected_bond_color(self, color: Optional[str]) -> None:
        """Método auxiliar para  set selected bond color.

        Args:
            color: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        bond_ids = list(self.state.selected_bonds)
        if not bond_ids:
            return
        self.undo_stack.beginMacro("Change bond color")
        for bond_id in bond_ids:
            cmd = ChangeBondColorCommand(self.model, self, bond_id, color)
            self.undo_stack.push(cmd)
        self.undo_stack.endMacro()

    def _adjust_selected_bond_stroke(self, delta: float) -> None:
        """Método auxiliar para  adjust selected bond stroke.

        Args:
            delta: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        bond_ids = list(self.state.selected_bonds)
        if not bond_ids:
            return
        default_stroke = self.drawing_style.stroke_px
        self.undo_stack.beginMacro("Change bond thickness")
        for bond_id in bond_ids:
            bond = self.model.get_bond(bond_id)
            current = bond.stroke_px if bond.stroke_px is not None else default_stroke
            new_value = max(0.6, current + delta)
            if abs(new_value - default_stroke) < 0.05:
                new_value = None
            cmd = ChangeBondStrokeCommand(self.model, self, bond_id, new_value)
            self.undo_stack.push(cmd)
        self.undo_stack.endMacro()

    def _reset_selected_bond_stroke(self) -> None:
        """Método auxiliar para  reset selected bond stroke.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        bond_ids = list(self.state.selected_bonds)
        if not bond_ids:
            return
        self.undo_stack.beginMacro("Reset bond thickness")
        for bond_id in bond_ids:
            cmd = ChangeBondStrokeCommand(self.model, self, bond_id, None)
            self.undo_stack.push(cmd)
        self.undo_stack.endMacro()


    def _ensure_selection_overlay(self) -> None:
        """Método auxiliar para  ensure selection overlay.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
            handle = QGraphicsEllipseItem(-radius, -radius, radius * 2.0, radius * 2.0)
            handle.setBrush(QBrush(QColor("#4A90D9")))
            handle.setPen(QPen(QColor("#4A90D9"), 0))
            handle.setZValue(46)
            handle.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIgnoresTransformations, True)
            self.scene.addItem(handle)
            self._selection_handle = handle

        if self._selection_move_handle is None:
            radius = SELECTION_HANDLE_RADIUS_PX
            handle = QGraphicsEllipseItem(-radius, -radius, radius * 2.0, radius * 2.0)
            handle.setBrush(QBrush(QColor("#FFFFFF")))
            handle.setPen(QPen(QColor("#4A90D9"), 1))
            handle.setZValue(46)
            handle.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIgnoresTransformations, True)
            handle.setCursor(Qt.CursorShape.OpenHandCursor)
            self.scene.addItem(handle)
            self._selection_move_handle = handle

        if self._selection_scale_handle is None:
            size = SELECTION_HANDLE_RADIUS_PX * 2.0
            handle = QGraphicsRectItem(-size / 2.0, -size / 2.0, size, size)
            handle.setBrush(QBrush(QColor("#4A90D9")))
            handle.setPen(QPen(QColor("#4A90D9"), 0))
            handle.setZValue(46)
            handle.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIgnoresTransformations, True)
            handle.setCursor(Qt.CursorShape.SizeFDiagCursor)
            self.scene.addItem(handle)
            self._selection_scale_handle = handle

    def _clear_selection_overlay(self) -> None:
        """Método auxiliar para  clear selection overlay.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  selected items bbox.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        rect: Optional[QRectF] = None

        def extend(candidate: QRectF) -> None:
            """Método auxiliar para extend.

            Args:
                candidate: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
            nonlocal rect
            if not candidate.isValid() or candidate.isNull():
                return
            rect = candidate if rect is None else rect.united(candidate)

        def extend_atom_bounds(atom_id: int) -> None:
            """Método auxiliar para extend atom bounds.

            Args:
                atom_id: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
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
            if isinstance(item, (ArrowItem, BracketItem, TextAnnotationItem, WavyAnchorItem)):
                extend(item.sceneBoundingRect())
        return rect

    def _selected_atom_ids_for_transform(self) -> set[int]:
        """Método auxiliar para  selected atom ids for transform.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atom_ids = set(self.state.selected_atoms)
        for bond_id in self.state.selected_bonds:
            if bond_id in self.model.bonds:
                bond = self.model.get_bond(bond_id)
                atom_ids.add(bond.a1_id)
                atom_ids.add(bond.a2_id)
        return atom_ids

    def _update_selection_overlay(self) -> None:
        """Método auxiliar para  update selection overlay.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        bbox = self._selected_items_bbox()
        if bbox is None:
            for attr in (
                "_selection_box",
                "_selection_handle",
                "_selection_move_handle",
                "_selection_scale_handle",
            ):
                item = getattr(self, attr)
                if item is None:
                    continue
                try:
                    item.setVisible(False)
                except RuntimeError:
                    setattr(self, attr, None)
            return
        self._ensure_selection_overlay()
        padded = QRectF(bbox)
        pad = max(2.0, float(self.drawing_style.stroke_px))
        padded.adjust(-pad, -pad, pad, pad)

        def offset_in_scene(base: QPointF, dx_view: float, dy_view: float) -> QPointF:
            """Método auxiliar para offset in scene.

            Args:
                base: Descripción del parámetro.
                dx_view: Descripción del parámetro.
                dy_view: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
            view_pt = self.mapFromScene(base)
            view_x = float(view_pt.x()) + dx_view
            view_y = float(view_pt.y()) + dy_view
            view_pt = QPoint(int(round(view_x)), int(round(view_y)))
            return self.mapToScene(view_pt)

        if self._selection_box is not None:
            try:
                self._selection_box.setRect(padded)
                self._selection_box.setVisible(True)
            except RuntimeError:
                self._selection_box = None
        if self._selection_handle is not None:
            top_center = QPointF(padded.center().x(), padded.top())
            handle_pos = offset_in_scene(top_center, 0.0, -SELECTION_ROTATE_OFFSET_PX)
            try:
                self._selection_handle.setPos(handle_pos)
                self._selection_handle.setVisible(True)
            except RuntimeError:
                self._selection_handle = None

        if self._selection_move_handle is not None:
            top_center = QPointF(padded.center().x(), padded.top())
            handle_pos = offset_in_scene(top_center, 0.0, SELECTION_MOVE_OFFSET_PX)
            try:
                self._selection_move_handle.setPos(handle_pos)
                self._selection_move_handle.setVisible(True)
            except RuntimeError:
                self._selection_move_handle = None

        if self._selection_scale_handle is not None:
            corner = QPointF(padded.right(), padded.bottom())
            handle_pos = offset_in_scene(
                corner, -SELECTION_HANDLE_RADIUS_PX, -SELECTION_HANDLE_RADIUS_PX
            )
            try:
                self._selection_scale_handle.setPos(handle_pos)
                self._selection_scale_handle.setVisible(True)
            except RuntimeError:
                self._selection_scale_handle = None

    def _hit_selection_handle(self, scene_pos: QPointF) -> bool:
        """Método auxiliar para  hit selection handle.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._selection_handle is None:
            return False
        try:
            if not self._selection_handle.isVisible():
                return False
        except RuntimeError:
            self._selection_handle = None
            return False
        return self._hit_handle_item(self._selection_handle, scene_pos)

    def _hit_selection_move_handle(self, scene_pos: QPointF) -> bool:
        """Método auxiliar para  hit selection move handle.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._selection_move_handle is None:
            return False
        try:
            if not self._selection_move_handle.isVisible():
                return False
        except RuntimeError:
            self._selection_move_handle = None
            return False
        return self._hit_handle_item(self._selection_move_handle, scene_pos)

    def _hit_selection_scale_handle(self, scene_pos: QPointF) -> bool:
        """Método auxiliar para  hit selection scale handle.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._selection_scale_handle is None:
            return False
        try:
            if not self._selection_scale_handle.isVisible():
                return False
        except RuntimeError:
            self._selection_scale_handle = None
            return False
        return self._hit_handle_item(self._selection_scale_handle, scene_pos)

    def _hit_handle_item(self, handle: QGraphicsItem, scene_pos: QPointF) -> bool:
        # Use a screen-space hit target so handles stay clickable at low zoom.
        """Método auxiliar para  hit handle item.

        Args:
            handle: Descripción del parámetro.
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        view_pos = self.mapFromScene(scene_pos)
        center = handle.pos() + handle.boundingRect().center()
        center_view = self.mapFromScene(center)
        dx = view_pos.x() - center_view.x()
        dy = view_pos.y() - center_view.y()
        radius = max(SELECTION_HANDLE_RADIUS_PX * 2.0, 12.0)
        return (dx * dx + dy * dy) <= (radius * radius)

    def _trackball_atom_ids(self) -> tuple[int, ...]:
        """Devuelve IDs objetivo para trackball (solo selección activa)."""
        atom_ids = self._selected_atom_ids_for_transform()
        return tuple(sorted(atom_id for atom_id in atom_ids if atom_id in self.model.atoms))

    def _clear_trackball_reference(self) -> None:
        """Limpia la referencia pseudo-3D cuando deja de representar la geometría actual."""
        self._rotation_3d_ref_atom_ids = tuple()
        self._rotation_3d_ref_positions = {}
        self._rotation_3d_pitch_deg = 0.0
        self._rotation_3d_yaw_deg = 0.0

    def _clear_trackball_reference_if_desynced(self, atom_ids: set[int]) -> None:
        """Invalida trackball si un movimiento 2D rompe la referencia afín."""
        if self._is_rotating_3d:
            return
        if not atom_ids:
            return
        ref_atom_ids = self._rotation_3d_ref_atom_ids
        if not ref_atom_ids or not self._rotation_3d_ref_positions:
            return
        ref_set = set(ref_atom_ids)
        moved_tracked = ref_set.intersection(atom_ids)
        if not moved_tracked:
            return
        if any(atom_id not in self.model.atoms for atom_id in ref_atom_ids):
            self._clear_trackball_reference()
            return
        if any(atom_id not in self._rotation_3d_ref_positions for atom_id in ref_atom_ids):
            self._clear_trackball_reference()
            return

        projected = self._project_trackball_reference(
            ref_atom_ids,
            self._rotation_3d_pitch_deg,
            self._rotation_3d_yaw_deg,
        )
        for atom_id in moved_tracked:
            atom = self.model.get_atom(atom_id)
            expected_x, expected_y = projected[atom_id]
            if math.hypot(expected_x - atom.x, expected_y - atom.y) > 1e-4:
                self._clear_trackball_reference()
                return

    def _connected_component_atom_ids(self, seed_atom_id: int) -> set[int]:
        """Obtiene la componente conectada del átomo semilla."""
        if seed_atom_id not in self.model.atoms:
            return set()
        visited: set[int] = set()
        stack = [seed_atom_id]
        while stack:
            current = stack.pop()
            if current in visited:
                continue
            visited.add(current)
            for bond in self.model.bonds.values():
                if bond.a1_id == current and bond.a2_id not in visited:
                    stack.append(bond.a2_id)
                elif bond.a2_id == current and bond.a1_id not in visited:
                    stack.append(bond.a1_id)
        return visited

    def _precise_3d_target_atom_ids(
        self,
        clicked_atom_id: Optional[int] = None,
        clicked_bond_id: Optional[int] = None,
    ) -> tuple[int, ...]:
        """Resuelve el conjunto de átomos objetivo (solo selección activa)."""
        selected_atoms = self._selected_atom_ids_for_transform()
        if selected_atoms:
            return tuple(sorted(atom_id for atom_id in selected_atoms if atom_id in self.model.atoms))

        selected_bonds = {
            bond_id for bond_id in self.state.selected_bonds if bond_id in self.model.bonds
        }
        if selected_bonds:
            atom_ids: set[int] = set()
            for bond_id in selected_bonds:
                bond = self.model.get_bond(bond_id)
                if bond.a1_id in self.model.atoms:
                    atom_ids.add(bond.a1_id)
                if bond.a2_id in self.model.atoms:
                    atom_ids.add(bond.a2_id)
            if atom_ids:
                return tuple(sorted(atom_ids))

        return tuple()

    def _ensure_trackball_reference(
        self,
        atom_ids: tuple[int, ...],
        before_positions: Dict[int, Tuple[float, float]],
    ) -> None:
        """Garantiza referencia estable para proyecciones trackball."""
        reset_reference = atom_ids != self._rotation_3d_ref_atom_ids
        if not reset_reference:
            for atom_id in atom_ids:
                if atom_id not in self._rotation_3d_ref_positions:
                    reset_reference = True
                    break
        if not reset_reference:
            projected = self._project_trackball_reference(
                atom_ids,
                self._rotation_3d_pitch_deg,
                self._rotation_3d_yaw_deg,
            )
            max_delta = 0.0
            for atom_id in atom_ids:
                px, py = projected[atom_id]
                cx, cy = before_positions[atom_id]
                max_delta = max(max_delta, math.hypot(px - cx, py - cy))
            if max_delta > TRACKBALL_REFERENCE_MATCH_TOLERANCE_PX:
                # La geometría cambió por otro flujo (mover, editar, undo, etc.).
                reset_reference = True
        if reset_reference:
            self._rotation_3d_ref_atom_ids = atom_ids
            self._rotation_3d_ref_positions = dict(before_positions)
            self._rotation_3d_pitch_deg = 0.0
            self._rotation_3d_yaw_deg = 0.0

    def _project_trackball_reference(
        self,
        atom_ids: tuple[int, ...],
        pitch_deg: float,
        yaw_deg: float,
    ) -> Dict[int, Tuple[float, float]]:
        """Proyecta la referencia 2D con ángulos pseudo-3D absolutos."""
        if not atom_ids:
            return {}
        points = [self._rotation_3d_ref_positions[atom_id] for atom_id in atom_ids]
        count = len(points)
        cx = sum(x for x, _ in points) / count
        cy = sum(y for _, y in points) / count
        rotated = project_3d_rotation(
            points,
            (cx, cy),
            math.radians(pitch_deg),
            math.radians(yaw_deg),
        )
        return {atom_id: point for atom_id, point in zip(atom_ids, rotated)}

    def _begin_3d_rotation_drag(
        self,
        scene_pos: QPointF,
        click_atom_id: Optional[int] = None,
    ) -> bool:
        """Inicia rotación pseudo-3D usando una referencia estable para evitar colapso."""
        atom_ids = self._trackball_atom_ids()
        if not atom_ids:
            return False

        before_positions = {
            atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
            for atom_id in atom_ids
        }
        self._ensure_trackball_reference(atom_ids, before_positions)

        self._is_rotating_3d = True
        self._rotation_3d_start_view_pos = self.mapFromScene(scene_pos)
        self._rotation_3d_before_positions = before_positions
        self._rotation_3d_click_atom_id = click_atom_id
        self._rotation_3d_has_moved = False
        self._rotation_3d_drag_start_pitch_deg = self._rotation_3d_pitch_deg
        self._rotation_3d_drag_start_yaw_deg = self._rotation_3d_yaw_deg
        self.setCursor(Qt.CursorShape.ClosedHandCursor)
        return True

    def _update_3d_rotation_drag(self, scene_pos: QPointF) -> None:
        """Actualiza rotación pseudo-3D con ángulos absolutos sobre referencia fija."""
        if not self._is_rotating_3d or self._rotation_3d_start_view_pos is None:
            return

        atom_ids = self._rotation_3d_ref_atom_ids
        if not atom_ids:
            return

        view_pos = self.mapFromScene(scene_pos)
        dx = float(view_pos.x() - self._rotation_3d_start_view_pos.x())
        dy = float(view_pos.y() - self._rotation_3d_start_view_pos.y())

        pitch_deg = self._rotation_3d_drag_start_pitch_deg + dy * TRACKBALL_ROTATION_DEG_PER_PIXEL
        yaw_deg = self._rotation_3d_drag_start_yaw_deg + dx * TRACKBALL_ROTATION_DEG_PER_PIXEL
        pitch_deg = max(-TRACKBALL_MAX_TILT_DEG, min(TRACKBALL_MAX_TILT_DEG, pitch_deg))
        yaw_deg = max(-TRACKBALL_MAX_TILT_DEG, min(TRACKBALL_MAX_TILT_DEG, yaw_deg))

        projected = self._project_trackball_reference(atom_ids, pitch_deg, yaw_deg)
        changed = False
        moved_atom_ids: set[int] = set()
        for atom_id in atom_ids:
            if atom_id not in self.model.atoms:
                continue
            new_x, new_y = projected[atom_id]
            atom = self.model.get_atom(atom_id)
            if not changed and (abs(new_x - atom.x) > 1e-9 or abs(new_y - atom.y) > 1e-9):
                changed = True
            self.model.update_atom_position(atom_id, new_x, new_y)
            self.update_atom_item(atom_id, new_x, new_y)
            moved_atom_ids.add(atom_id)

        self._rotation_3d_pitch_deg = pitch_deg
        self._rotation_3d_yaw_deg = yaw_deg

        if changed and moved_atom_ids:
            self.update_bond_items_for_atoms(moved_atom_ids)
            self._rotation_3d_has_moved = True

    def _reset_3d_rotation_perspective(self) -> bool:
        """Restaura la perspectiva base del trackball para selección o toda la molécula."""
        selected_atom_ids = self._selected_atom_ids_for_transform()
        if selected_atom_ids:
            atom_ids = tuple(
                sorted(atom_id for atom_id in selected_atom_ids if atom_id in self.model.atoms)
            )
        elif self._rotation_3d_ref_atom_ids:
            atom_ids = self._rotation_3d_ref_atom_ids
        else:
            atom_ids = self._trackball_atom_ids()
        if not atom_ids:
            return False
        if atom_ids != self._rotation_3d_ref_atom_ids:
            return False
        if any(atom_id not in self.model.atoms for atom_id in atom_ids):
            return False
        if any(atom_id not in self._rotation_3d_ref_positions for atom_id in atom_ids):
            return False

        before = {
            atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
            for atom_id in atom_ids
        }
        after = {
            atom_id: self._rotation_3d_ref_positions[atom_id]
            for atom_id in atom_ids
        }
        changed = any(
            abs(before[atom_id][0] - after[atom_id][0]) > 1e-9
            or abs(before[atom_id][1] - after[atom_id][1]) > 1e-9
            for atom_id in atom_ids
        )
        self._rotation_3d_pitch_deg = 0.0
        self._rotation_3d_yaw_deg = 0.0
        self._rotation_3d_drag_start_pitch_deg = 0.0
        self._rotation_3d_drag_start_yaw_deg = 0.0
        if changed:
            self.undo_stack.push(MoveAtomsCommand(self.model, self, before, after))
        else:
            self._update_selection_overlay()
        return True

    def _finalize_3d_rotation_drag(self) -> None:
        """Finaliza la rotación pseudo-3D y registra un único comando de undo."""
        if not self._is_rotating_3d:
            return

        if self._rotation_3d_has_moved and self._rotation_3d_before_positions:
            before: Dict[int, Tuple[float, float]] = {}
            after: Dict[int, Tuple[float, float]] = {}
            changed = False
            for atom_id, (x0, y0) in self._rotation_3d_before_positions.items():
                if atom_id not in self.model.atoms:
                    continue
                atom = self.model.get_atom(atom_id)
                before[atom_id] = (x0, y0)
                after[atom_id] = (atom.x, atom.y)
                if not changed and (abs(atom.x - x0) > 1e-9 or abs(atom.y - y0) > 1e-9):
                    changed = True
            if changed and before:
                self.undo_stack.push(
                    MoveAtomsCommand(
                        self.model,
                        self,
                        before,
                        after,
                        skip_first_redo=True,
                    )
                )
        elif not self._rotation_3d_has_moved and self._rotation_3d_click_atom_id is not None:
            self._cycle_anchor_override(self._rotation_3d_click_atom_id)

        self._is_rotating_3d = False
        self._rotation_3d_start_view_pos = None
        self._rotation_3d_before_positions = {}
        self._rotation_3d_click_atom_id = None
        self._rotation_3d_has_moved = False
        self._rotation_3d_drag_start_pitch_deg = 0.0
        self._rotation_3d_drag_start_yaw_deg = 0.0
        if self.current_tool in {"tool_select", "tool_select_lasso"}:
            self.setCursor(Qt.CursorShape.ArrowCursor)
        self._update_selection_overlay()

    def _prompt_precise_3d_rotation(
        self,
        clicked_atom_id: Optional[int] = None,
        clicked_bond_id: Optional[int] = None,
    ) -> None:
        """Solicita ángulos precisos y aplica rotación pseudo-3D reproducible."""
        atom_ids = self._precise_3d_target_atom_ids(
            clicked_atom_id=clicked_atom_id,
            clicked_bond_id=clicked_bond_id,
        )
        if not atom_ids:
            return
        before = {
            atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
            for atom_id in atom_ids
            if atom_id in self.model.atoms
        }
        if not before:
            return
        atom_ids = tuple(sorted(before.keys()))
        self._ensure_trackball_reference(atom_ids, before)

        default_pitch = self._rotation_3d_pitch_deg
        default_yaw = self._rotation_3d_yaw_deg

        def apply_preview(pitch_deg: float, yaw_deg: float) -> None:
            """Aplica proyección en vivo para previsualización."""
            pitch_deg = max(-TRACKBALL_MAX_TILT_DEG, min(TRACKBALL_MAX_TILT_DEG, float(pitch_deg)))
            yaw_deg = max(-TRACKBALL_MAX_TILT_DEG, min(TRACKBALL_MAX_TILT_DEG, float(yaw_deg)))
            # Mantener estado trackball activo durante preview para que los
            # círculos aromáticos usen el mismo contexto afín que Alt+Drag.
            self._rotation_3d_pitch_deg = pitch_deg
            self._rotation_3d_yaw_deg = yaw_deg
            projected = self._project_trackball_reference(atom_ids, pitch_deg, yaw_deg)
            moved_atom_ids: set[int] = set()
            for atom_id in atom_ids:
                if atom_id not in self.model.atoms:
                    continue
                new_x, new_y = projected[atom_id]
                self.model.update_atom_position(atom_id, new_x, new_y)
                self.update_atom_item(atom_id, new_x, new_y)
                moved_atom_ids.add(atom_id)
            if moved_atom_ids:
                self.update_bond_items_for_atoms(moved_atom_ids)
            else:
                self._update_selection_overlay()

        dialog = TrackballRotationDialog(
            float(default_pitch),
            float(default_yaw),
            TRACKBALL_MAX_TILT_DEG,
            parent=self,
        )
        dialog.preview_changed.connect(apply_preview)
        result = dialog.exec()
        if result == QDialog.DialogCode.Accepted:
            pitch_deg, yaw_deg = dialog.angles()
            self._rotation_3d_pitch_deg = float(pitch_deg)
            self._rotation_3d_yaw_deg = float(yaw_deg)
            after = {
                atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
                for atom_id in atom_ids
                if atom_id in self.model.atoms
            }
            changed_ids = {
                atom_id
                for atom_id in after.keys()
                if (
                    abs(after[atom_id][0] - before[atom_id][0]) > 1e-9
                    or abs(after[atom_id][1] - before[atom_id][1]) > 1e-9
                )
            }
            if not changed_ids:
                self._update_selection_overlay()
                return
            self.undo_stack.push(
                MoveAtomsCommand(
                    self.model,
                    self,
                    {atom_id: before[atom_id] for atom_id in changed_ids},
                    {atom_id: after[atom_id] for atom_id in changed_ids},
                    skip_first_redo=True,
                )
            )
            return

        restored_ids: set[int] = set()
        for atom_id, (x0, y0) in before.items():
            if atom_id not in self.model.atoms:
                continue
            self.model.update_atom_position(atom_id, x0, y0)
            self.update_atom_item(atom_id, x0, y0)
            restored_ids.add(atom_id)
        self._rotation_3d_pitch_deg = default_pitch
        self._rotation_3d_yaw_deg = default_yaw
        if restored_ids:
            self.update_bond_items_for_atoms(restored_ids)
        else:
            self._update_selection_overlay()

    def _begin_rotation_drag(self, scene_pos: QPointF) -> None:
        """Método auxiliar para  begin rotation drag.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if (
            not self._selected_atom_ids_for_transform()
            and not self._selected_text_items()
            and not self._selected_arrow_items()
        ):
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
        self._rotation_start_arrow_positions = {
            item: (item.start_point(), item.end_point())
            for item in self._selected_arrow_items()
        }

    def _update_rotation_drag(self, scene_pos: QPointF) -> None:
        """Método auxiliar para  update rotation drag.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._rotation_center is None or self._rotation_start_angle is None:
            return
        dx = scene_pos.x() - self._rotation_center.x()
        dy = scene_pos.y() - self._rotation_center.y()
        current_angle = math.atan2(dy, dx)
        delta_angle = current_angle - self._rotation_start_angle
        degrees = math.degrees(delta_angle)

        self.rotate_selection_degrees(degrees, use_start_positions=True)

    def _finalize_rotation_drag(self) -> None:
        """Método auxiliar para  finalize rotation drag.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._rotation_dragging = False
        self._rotation_center = None
        self._rotation_start_angle = None
        self._rotation_start_positions = {}
        self._rotation_start_text_transforms = {}
        self._rotation_start_arrow_positions = {}
        self._update_selection_overlay()

    def _begin_scale_drag(self, scene_pos: QPointF) -> None:
        """Método auxiliar para  begin scale drag.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if (
            not self._selected_atom_ids_for_transform()
            and not self._selected_text_items()
            and not self._selected_arrow_items()
            and not self._selected_bracket_items()
        ):
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
        self._scale_start_arrow_positions = {
            item: (item.start_point(), item.end_point())
            for item in self._selected_arrow_items()
        }
        self._scale_start_bracket_rects = {
            item: item.base_rect()
            for item in self._selected_bracket_items()
        }
        self._scale_has_moved = False

    def _update_scale_drag(self, scene_pos: QPointF) -> None:
        """Método auxiliar para  update scale drag.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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

        if self._scale_start_arrow_positions:
            for item, (start, end) in self._scale_start_arrow_positions.items():
                nsx = self._scale_anchor.x() + (start.x() - self._scale_anchor.x()) * scale
                nsy = self._scale_anchor.y() + (start.y() - self._scale_anchor.y()) * scale
                nex = self._scale_anchor.x() + (end.x() - self._scale_anchor.x()) * scale
                ney = self._scale_anchor.y() + (end.y() - self._scale_anchor.y()) * scale
                item.update_positions(QPointF(nsx, nsy), QPointF(nex, ney))

        if self._scale_start_bracket_rects:
            for item, rect in self._scale_start_bracket_rects.items():
                x0 = self._scale_anchor.x() + (rect.left() - self._scale_anchor.x()) * scale
                y0 = self._scale_anchor.y() + (rect.top() - self._scale_anchor.y()) * scale
                x1 = self._scale_anchor.x() + (rect.right() - self._scale_anchor.x()) * scale
                y1 = self._scale_anchor.y() + (rect.bottom() - self._scale_anchor.y()) * scale
                item.set_rect(QRectF(x0, y0, x1 - x0, y1 - y0))

        self._update_selection_overlay()

    def _finalize_scale_drag(self) -> None:
        """Método auxiliar para  finalize scale drag.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if not self._scale_dragging:
            return

        if self._scale_has_moved:
            move_atoms = bool(self._scale_start_positions)
            move_text = bool(self._scale_start_text_positions)
            move_arrows = bool(self._scale_start_arrow_positions)
            move_brackets = bool(self._scale_start_bracket_rects)

            if sum([move_atoms, move_text, move_arrows, move_brackets]) > 1:
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

            if move_arrows:
                before = self._scale_start_arrow_positions
                after = {item: (item.start_point(), item.end_point()) for item in before.keys()}
                self.undo_stack.push(MoveArrowItemsCommand(self, before, after))

            if move_brackets:
                before = self._scale_start_bracket_rects
                after = {item: item.base_rect() for item in before.keys()}
                self.undo_stack.push(MoveBracketItemsCommand(self, before, after))

            if sum([move_atoms, move_text, move_arrows, move_brackets]) > 1:
                self.undo_stack.endMacro()

        self._scale_dragging = False
        self._scale_anchor = None
        self._scale_start_handle = None
        self._scale_start_length = 0.0
        self._scale_start_positions = {}
        self._scale_start_text_positions = {}
        self._scale_start_arrow_positions = {}
        self._scale_start_bracket_rects = {}
        self._scale_has_moved = False
        self._update_selection_overlay()

    def _selected_text_items(self) -> list[TextAnnotationItem]:
        """Método auxiliar para  selected text items.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        return [item for item in self.scene.selectedItems() if isinstance(item, TextAnnotationItem)]

    def _selected_arrow_items(self) -> list[ArrowItem]:
        """Método auxiliar para  selected arrow items.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        return [item for item in self.scene.selectedItems() if isinstance(item, ArrowItem)]

    def _selected_bracket_items(self) -> list[BracketItem]:
        """Método auxiliar para  selected bracket items.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        return [item for item in self.scene.selectedItems() if isinstance(item, BracketItem)]

    def rotate_selection_degrees(self, degrees: float, use_start_positions: bool = False) -> None:
        """Rotate selected items by degrees around their collective center."""
        selected_atom_ids = self._selected_atom_ids_for_transform()
        selected_text_items = self._selected_text_items()
        selected_arrow_items = self._selected_arrow_items()
        
        if not selected_atom_ids and not selected_text_items and not selected_arrow_items:
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

        def rotate_point(pt: QPointF) -> QPointF:
            dx = pt.x() - cx
            dy = pt.y() - cy
            return QPointF(dx * cos_a - dy * sin_a + cx, dx * sin_a + dy * cos_a + cy)

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

        # Rotate Arrow Items
        if selected_arrow_items:
            if use_start_positions:
                for item in selected_arrow_items:
                    if item in self._rotation_start_arrow_positions:
                        start_pt, end_pt = self._rotation_start_arrow_positions[item]
                    else:
                        start_pt, end_pt = item.start_point(), item.end_point()
                    item.update_positions(rotate_point(start_pt), rotate_point(end_pt))
            else:
                before = {
                    item: (item.start_point(), item.end_point())
                    for item in selected_arrow_items
                }
                after = {
                    item: (rotate_point(start), rotate_point(end))
                    for item, (start, end) in before.items()
                }
                if before:
                    self.undo_stack.push(MoveArrowItemsCommand(self, before, after))

        # Update overlay at the end
        if use_start_positions:
             self._update_selection_overlay()

    def flip_selection_horizontal(self) -> None:
        """Método auxiliar para flip selection horizontal.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para flip selection vertical.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  set bond anchor.

        Args:
            atom_id: Descripción del parámetro.
            reset_angle: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  clear bond anchor.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self.bond_anchor_id is not None and self.bond_anchor_id in self.atom_items:
            if self.bond_anchor_id not in self.state.selected_atoms:
                self.atom_items[self.bond_anchor_id].set_selected(False)
        self.bond_anchor_id = None
        self._bond_last_angle = None
        self._bond_zigzag_sign = 1

    def _infer_coordination_donor_atom(
        self,
        a1_id: int,
        a2_id: int,
        preferred: Optional[int] = None,
    ) -> int:
        """Infere el átomo donador para enlaces coordinativos."""
        if preferred in {a1_id, a2_id}:
            return int(preferred)
        atom1 = self.model.get_atom(a1_id)
        atom2 = self.model.get_atom(a2_id)
        a1_is_center = bool(getattr(atom1, "is_coordination_center", False))
        a2_is_center = bool(getattr(atom2, "is_coordination_center", False))
        if a1_is_center and not a2_is_center:
            return a2_id
        if a2_is_center and not a1_is_center:
            return a1_id
        # Si no hay centro claro, conservar dirección del primer click.
        return a1_id

    @staticmethod
    def _parse_bond_style_payload(bond_data: dict) -> BondStyle:
        """Resuelve estilo de enlace para payloads nuevos y legados."""
        raw_style = bond_data.get("style")
        if raw_style is None:
            raw_style = bond_data.get("type", BondStyle.PLAIN.value)
        try:
            return BondStyle(raw_style)
        except Exception:
            return BondStyle.PLAIN

    @staticmethod
    def _parse_bond_stereo_payload(bond_data: dict) -> BondStereo:
        """Resuelve estereo de enlace para payloads con fallback seguro."""
        raw_stereo = bond_data.get("stereo", BondStereo.NONE.value)
        try:
            return BondStereo(raw_stereo)
        except Exception:
            return BondStereo.NONE

    def _create_or_update_bond(
        self,
        a1_id: int,
        a2_id: int,
        order: int,
        style: BondStyle,
        stereo: BondStereo,
        is_aromatic: bool,
        flex_curve_1: Optional[float] = None,
        flex_curve_2: Optional[float] = None,
    ) -> None:
        """Método auxiliar para  create or update bond.

        Args:
            a1_id: Descripción del parámetro.
            a2_id: Descripción del parámetro.
            order: Descripción del parámetro.
            style: Descripción del parámetro.
            stereo: Descripción del parámetro.
            is_aromatic: Descripción del parámetro.
            flex_curve_1: Curvatura normalizada de control 1 (estilo FLEX).
            flex_curve_2: Curvatura normalizada de control 2 (estilo FLEX).

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
                if existing.is_aromatic and not is_aromatic and style != BondStyle.PLAIN:
                    # Mantener aromaticidad al aplicar estilos visuales no planos
                    # (p. ej. bold para proyecciones tipo Haworth).
                    new_is_aromatic = True
            new_donor_atom_id = None
            if new_style == BondStyle.COORDINATION:
                new_donor_atom_id = self._infer_coordination_donor_atom(
                    existing.a1_id,
                    existing.a2_id,
                    preferred=getattr(existing, "donor_atom_id", None),
                )
            new_flex_curve_1 = None
            new_flex_curve_2 = None
            if new_style == BondStyle.FLEX:
                if flex_curve_1 is not None or flex_curve_2 is not None:
                    new_flex_curve_1 = (
                        float(flex_curve_1) if flex_curve_1 is not None else None
                    )
                    new_flex_curve_2 = (
                        float(flex_curve_2) if flex_curve_2 is not None else None
                    )
                elif getattr(existing, "style", None) == BondStyle.FLEX:
                    new_flex_curve_1 = getattr(existing, "flex_curve_1", None)
                    new_flex_curve_2 = getattr(existing, "flex_curve_2", None)
            cmd = ChangeBondCommand(
                self.model,
                self,
                existing.id,
                new_order=new_order,
                new_style=new_style,
                new_stereo=new_stereo,
                new_is_aromatic=new_is_aromatic,
                new_donor_atom_id=new_donor_atom_id,
                new_flex_curve_1=new_flex_curve_1,
                new_flex_curve_2=new_flex_curve_2,
            )
            self.undo_stack.push(cmd)
            if new_is_aromatic:
                self._kekulize_aromatic_bonds(seed_atoms={a1_id, a2_id})
            return

        donor_atom_id = None
        if style == BondStyle.COORDINATION:
            donor_atom_id = self._infer_coordination_donor_atom(a1_id, a2_id, preferred=a1_id)
        cmd = AddBondCommand(
            self.model,
            self,
            a1_id,
            a2_id,
            order,
            style,
            stereo,
            is_aromatic=is_aromatic,
            donor_atom_id=donor_atom_id,
            flex_curve_1=flex_curve_1,
            flex_curve_2=flex_curve_2,
        )
        self.undo_stack.push(cmd)
        if is_aromatic:
            self._kekulize_aromatic_bonds(seed_atoms={a1_id, a2_id})

    def _add_bond_between(self, a1_id: int, a2_id: int) -> None:
        """Método auxiliar para  add bond between.

        Args:
            a1_id: Descripción del parámetro.
            a2_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        order = self.state.active_bond_order
        style = self.state.active_bond_style
        stereo = self.state.active_bond_stereo
        is_aromatic = self.state.active_bond_aromatic
        if style != BondStyle.PLAIN or is_aromatic:
            order = 1
        self._create_or_update_bond(a1_id, a2_id, order, style, stereo, is_aromatic)

    def _cycle_bond_order(self, bond_id: int) -> None:
        """Método auxiliar para  cycle bond order.

        Args:
            bond_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
                new_donor_atom_id=None,
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
        """Método auxiliar para  apply bond style.

        Args:
            bond_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if bond_id not in self.model.bonds:
            return
        bond = self.model.get_bond(bond_id)
        order = self.state.active_bond_order
        style = self.state.active_bond_style
        stereo = self.state.active_bond_stereo
        requested_aromatic = self.state.active_bond_aromatic
        is_aromatic = requested_aromatic
        if bond.is_aromatic and not requested_aromatic and style != BondStyle.PLAIN:
            # Preservar aromaticidad si solo se cambia el estilo visual.
            is_aromatic = True
        if style != BondStyle.PLAIN:
            order = 1
        if requested_aromatic:
            order = 1
            style = BondStyle.PLAIN
            stereo = BondStereo.NONE
        donor_atom_id = None
        if style == BondStyle.COORDINATION:
            donor_atom_id = self._infer_coordination_donor_atom(
                bond.a1_id,
                bond.a2_id,
                preferred=getattr(bond, "donor_atom_id", None),
            )
        cmd = ChangeBondCommand(
            self.model,
            self,
            bond_id,
            new_order=order,
            new_style=style,
            new_stereo=stereo,
            new_is_aromatic=is_aromatic,
            new_donor_atom_id=donor_atom_id,
        )
        self.undo_stack.push(cmd)
        if is_aromatic:
            self._kekulize_aromatic_bonds(seed_atoms={bond.a1_id, bond.a2_id})

    def _create_first_bond(self, scene_pos: QPointF, modifiers: Qt.KeyboardModifiers) -> None:
        """Método auxiliar para  create first bond.

        Args:
            scene_pos: Descripción del parámetro.
            modifiers: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  extend bond from anchor.

        Args:
            scene_pos: Descripción del parámetro.
            modifiers: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  add bond with new atom.

        Args:
            anchor_id: Descripción del parámetro.
            angle: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  record bond angle between.

        Args:
            a1_id: Descripción del parámetro.
            a2_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        a1 = self.model.get_atom(a1_id)
        a2 = self.model.get_atom(a2_id)
        angle = math.atan2(a2.y - a1.y, a2.x - a1.x)
        self._record_bond_angle(angle)

    def _record_bond_angle(self, angle: float) -> None:
        """Método auxiliar para  record bond angle.

        Args:
            angle: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._bond_last_angle = self._normalize_angle(angle)

    def _default_bond_angle(self, anchor_id: int) -> float:
        """Método auxiliar para  default bond angle.

        Args:
            anchor_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  reset zigzag.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._bond_zigzag_sign = 1

    def _snap_angle_to_environment(self, angle: float, anchor_id: int, step: float) -> float:
        """Método auxiliar para  snap angle to environment.

        Args:
            angle: Descripción del parámetro.
            anchor_id: Descripción del parámetro.
            step: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  snap angle.

        Args:
            angle: Descripción del parámetro.
            step: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if step <= 0:
            return angle
        return self._normalize_angle(round(angle / step) * step)

    def _bond_environment_step(self, anchor_id: int) -> float:
        """Método auxiliar para  bond environment step.

        Args:
            anchor_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        max_order = max(1, self.state.active_bond_order)
        for bond in self.model.bonds.values():
            if bond.a1_id == anchor_id or bond.a2_id == anchor_id:
                order = 2 if bond.is_aromatic else bond.order
                max_order = max(max_order, order)
        if max_order >= 3:
            return math.pi
        if max_order == 2:
            return 2 * math.pi / 3
        if self.state.fixed_angles:
            step_deg = self._angle_snap_step_deg()
            if step_deg > 0:
                return math.radians(snap_angle_deg(SP3_BOND_ANGLE_DEG, step_deg))
        return math.radians(self._sp3_display_angle_deg())

    def _get_anchor_bond_angles(self, anchor_id: int) -> List[float]:
        """Método auxiliar para  get anchor bond angles.

        Args:
            anchor_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  best separated angle.

        Args:
            candidates: Descripción del parámetro.
            existing: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        best = None
        best_sep = -1.0
        for candidate in candidates:
            sep = min(self._angle_distance(candidate, a) for a in existing)
            if sep > best_sep:
                best_sep = sep
                best = candidate
        return self._normalize_angle(best) if best is not None else 0.0

    def _normalize_angle(self, angle: float) -> float:
        """Método auxiliar para  normalize angle.

        Args:
            angle: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        return (angle + math.pi * 2) % (math.pi * 2)

    def _angle_distance(self, a: float, b: float) -> float:
        """Método auxiliar para  angle distance.

        Args:
            a: Descripción del parámetro.
            b: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        diff = (a - b + math.pi) % (2 * math.pi) - math.pi
        return abs(diff)

    def _angle_snap_step_deg(self) -> float:
        """Método auxiliar para  angle snap step deg.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if not self.state.fixed_angles:
            return 0.0
        step = float(self.state.angle_step_deg)
        return step if step > 0 else 0.0

    def _snap_angles_to_grid(self, angles: Iterable[float]) -> list[float]:
        """Método auxiliar para  snap angles to grid.

        Args:
            angles: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  sp3 display angle deg.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        # Mantener 109.5° real para geometría sp3.
        # Si se ajusta a una retícula de 30°, se vuelve 120° (trigonal) y
        # termina bloqueando el 4º enlace en carbonos tetravalentes.
        return SP3_BOND_ANGLE_DEG

    # -------------------------------------------------------------------------
    # Hover + Drag State
    # -------------------------------------------------------------------------
    def _pick_hover_target(self, scene_pos: QPointF) -> tuple[Optional[int], Optional[int]]:
        """Método auxiliar para  pick hover target.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  update hover.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  cancel drag.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._drag_mode = "none"
        self._drag_anchor = None
        self._flex_bond_pending = None
        self._ring_last_vertices = None
        self._chain_last_points = None
        self._drag_free_orientation = False
        self._bond_drag_start = None
        self._rotation_dragging = False
        self._scale_dragging = False
        self._rotation_center = None
        self._rotation_start_angle = None
        self._rotation_start_positions = {}
        self._rotation_start_arrow_positions = {}
        self._is_rotating_3d = False
        self._rotation_3d_start_view_pos = None
        self._rotation_3d_before_positions = {}
        self._rotation_3d_click_atom_id = None
        self._rotation_3d_has_moved = False
        self._rotation_3d_drag_start_pitch_deg = 0.0
        self._rotation_3d_drag_start_yaw_deg = 0.0
        self._clear_selection_drag()
        if self._overlays_ready:
            self._optimize_zone.hide_zone()
            self._preview_bond_item.hide_preview()
            self._preview_ring_item.hide_preview()
            self._preview_chain_item.hide_preview()
            self._preview_chain_label.hide_label()
            self._preview_arrow_item.hide_preview()

    def _on_undo_stack_changed(self, _index: int) -> None:
        """Método auxiliar para  on undo stack changed.

        Args:
            _index: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        self.validate_structure()
        self._update_hover(self._last_scene_pos)
        self._update_selection_overlay()

    def _sync_scene_with_model(self) -> None:
        """Método auxiliar para  sync scene with model.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._suspend_numbering_refresh = True
        try:
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
        finally:
            self._suspend_numbering_refresh = False
        self.recompute_numbering()

    def show_valence_errors(self, errors: Iterable[int]) -> None:
        """Método auxiliar para show valence errors.

        Args:
            errors: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        error_ids = set(errors)
        for atom_id, item in self.atom_items.items():
            item.set_valence_error(atom_id in error_ids)

    def validate_structure(self) -> list[int]:
        """Valida estructura y aplica resaltado de valencias inválidas."""
        errors = list(self.model.validate())
        self.show_valence_errors(errors)
        return errors

    def _get_anchor_bond_angles_deg(self, anchor_id: int) -> list[float]:
        """Método auxiliar para  get anchor bond angles deg.

        Args:
            anchor_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  bond geometry.

        Args:
            anchor_id: Descripción del parámetro.
            bond_order: Descripción del parámetro.
            is_aromatic: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  incoming angle deg.

        Args:
            anchor_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        angles = self._get_anchor_bond_angles_deg(anchor_id)
        if len(angles) == 1:
            return angles[0]
        return None

    def _direction_collision_metrics(
        self, anchor_id: Optional[int], p0: QPointF, p1: QPointF
    ) -> tuple[int, float, float]:
        """Método auxiliar para  direction collision metrics.

        Args:
            anchor_id: Descripción del parámetro.
            p0: Descripción del parámetro.
            p1: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  direction score.

        Args:
            anchor_id: Descripción del parámetro.
            p0: Descripción del parámetro.
            p1: Descripción del parámetro.
            angle_deg_value: Descripción del parámetro.
            mouse_angle_deg: Descripción del parámetro.
            preferred_deg: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  direction is valid.

        Args:
            anchor_id: Descripción del parámetro.
            p0: Descripción del parámetro.
            p1: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  pick bond direction deg.

        Args:
            anchor: Descripción del parámetro.
            anchor_id: Descripción del parámetro.
            mouse_angle_deg: Descripción del parámetro.
            bond_order: Descripción del parámetro.
            is_aromatic: Descripción del parámetro.
            length: Descripción del parámetro.
            apply_collisions: Descripción del parámetro.
            allow_length_boost: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        existing_angles = self._get_anchor_bond_angles_deg(anchor_id) if anchor_id else []
        geometry = (
            self._bond_geometry(anchor_id, bond_order, is_aromatic)
            if anchor_id is not None
            else geometry_for_bond(bond_order, is_aromatic, [])
        )
        # Solo forzamos geometría tetraédrica exacta (109.5°) cuando el átomo
        # ya está congestionado (p. ej., intentando crear el 4º enlace).
        sp3_exact_mode = geometry == "sp3" and anchor_id is not None and len(existing_angles) >= 3
        candidates = candidate_directions_deg(geometry, existing_angles, mouse_angle_deg)
        if self.state.fixed_angles and not sp3_exact_mode:
            candidates = self._snap_angles_to_grid(candidates)
        occupied_tolerance = (
            SP3_OCCUPIED_TOLERANCE_DEG if sp3_exact_mode else ANGLE_OCCUPIED_TOLERANCE_DEG
        )
        candidates = filter_occupied_angles_deg(
            candidates, existing_angles, occupied_tolerance
        )
        if not candidates:
            candidates = candidate_directions_deg(geometry, [], mouse_angle_deg)
            if self.state.fixed_angles and not sp3_exact_mode:
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
                if self.state.fixed_angles and not sp3_exact_mode:
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
        """Método auxiliar para  select preferred angle.

        Args:
            anchor_id: Descripción del parámetro.
            cursor_angle: Descripción del parámetro.
            bond_order: Descripción del parámetro.
            is_aromatic: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  find overlapping bond.

        Args:
            p0: Descripción del parámetro.
            p1: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        for bond in self.model.bonds.values():
            a1 = self.model.get_atom(bond.a1_id)
            a2 = self.model.get_atom(bond.a2_id)
            b0 = QPointF(a1.x, a1.y)
            b1 = QPointF(a2.x, a2.y)
            if segments_nearly_equal(p0, p1, b0, b1, BOND_OVERLAP_TOLERANCE_PX):
                return bond.id
        return None

    def _snap_ring_vertex(self, pos: QPointF, excluded: set[int]) -> Optional[int]:
        """Método auxiliar para  snap ring vertex.

        Args:
            pos: Descripción del parámetro.
            excluded: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        candidates = [
            (atom.id, atom.x, atom.y)
            for atom in self.model.atoms.values()
            if atom.id not in excluded
        ]
        return closest_atom(pos, candidates, ATOM_HIT_RADIUS)

    def _build_ring_vertex_defs(
        self, vertices: List[QPointF], anchor_type: str, anchor_id: Optional[int]
    ) -> List[Tuple[Optional[int], float, float]]:
        """Método auxiliar para  build ring vertex defs.

        Args:
            vertices: Descripción del parámetro.
            anchor_type: Descripción del parámetro.
            anchor_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  bond drag distance.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._bond_drag_start is None:
            return 0.0
        return (scene_pos - self._bond_drag_start).manhattanLength()

    def _update_bond_drag_state(self, scene_pos: QPointF) -> None:
        """Método auxiliar para  update bond drag state.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._drag_free_orientation = (
            self._bond_drag_start is not None
            and self._bond_drag_distance(scene_pos) >= BOND_DRAG_THRESHOLD_PX
        )

    def _should_use_default_bond_angle(
        self, modifiers: Qt.KeyboardModifiers, scene_pos: Optional[QPointF] = None
    ) -> bool:
        """Método auxiliar para  should use default bond angle.

        Args:
            modifiers: Descripción del parámetro.
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  current bond geometry.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        order = self.state.active_bond_order
        is_aromatic = self.state.active_bond_aromatic
        if self.state.active_bond_style != BondStyle.PLAIN or is_aromatic:
            order = 1
        return order, is_aromatic

    def _peek_default_bond_angle(self, anchor_id: Optional[int]) -> float:
        """Método auxiliar para  peek default bond angle.

        Args:
            anchor_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  compute default bond endpoint.

        Args:
            anchor: Descripción del parámetro.
            anchor_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
            allow_length_boost=False,
        )
        p1 = endpoint_from_angle_len(anchor, theta, final_length)
        snap_id = closest_atom(
            p1,
            [(atom.id, atom.x, atom.y) for atom in self.model.atoms.values()],
            ATOM_HIT_RADIUS,
        )
        if snap_id is not None and (anchor_id is None or snap_id != anchor_id):
            # Evita auto-snap al átomo que ya está enlazado con el ancla:
            # en sp3 congestionado (3 enlaces), ese snap fuerza superposición
            # y termina ciclando el enlace existente en vez de crear el 4º.
            if anchor_id is not None and self.model.find_bond_between(anchor_id, snap_id) is not None:
                return p1
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
        """Método auxiliar para  update bond angle state.

        Args:
            p0: Descripción del parámetro.
            p1: Descripción del parámetro.
            used_default: Descripción del parámetro.
            anchor_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  begin place bond.

        Args:
            anchor_id: Descripción del parámetro.
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  update bond preview.

        Args:
            scene_pos: Descripción del parámetro.
            modifiers: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        if self.state.active_bond_style == BondStyle.FLEX:
            curve_1, curve_2 = self._flex_curve_from_pointer(p0, p1, scene_pos, modifiers)
            self._preview_bond_item.update_flex_curve(p0, p1, curve_1, curve_2)
        else:
            self._preview_bond_item.update_line(p0, p1)

    def _flex_curve_from_pointer(
        self,
        start: QPointF,
        end: QPointF,
        pointer: QPointF,
        modifiers: Qt.KeyboardModifiers,
    ) -> tuple[float, float]:
        """Calcula curvaturas normalizadas para enlaces FLEX."""
        dx = end.x() - start.x()
        dy = end.y() - start.y()
        length = math.hypot(dx, dy)
        if length <= 1e-6:
            return 0.22, 0.22
        ux = dx / length
        uy = dy / length
        nx = -uy
        ny = ux
        mid_x = (start.x() + end.x()) * 0.5
        mid_y = (start.y() + end.y()) * 0.5
        signed = ((pointer.x() - mid_x) * nx + (pointer.y() - mid_y) * ny) / length
        magnitude = min(0.95, max(0.12, abs(signed)))
        sign = -1.0 if signed < 0.0 else 1.0
        if abs(signed) < 0.02:
            sign = 1.0
        c1 = sign * magnitude
        # Shift activa curva compleja en "S" (controles opuestos).
        if modifiers & Qt.KeyboardModifier.ShiftModifier:
            c2 = -c1
        else:
            c2 = c1
        return c1, c2

    def _begin_flex_adjust_mode(
        self,
        p0: QPointF,
        p1: QPointF,
        anchor_id: Optional[int],
        target_id: Optional[int],
        use_default: bool,
        modifiers: Qt.KeyboardModifiers,
    ) -> None:
        """Activa ajuste interactivo de curvatura para enlace FLEX."""
        curve_1, curve_2 = self._flex_curve_from_pointer(p0, p1, self._last_scene_pos, modifiers)
        self._flex_bond_pending = {
            "p0": QPointF(p0),
            "p1": QPointF(p1),
            "anchor_id": anchor_id,
            "target_id": target_id,
            "use_default": bool(use_default),
            "curve_1": float(curve_1),
            "curve_2": float(curve_2),
        }
        self._drag_mode = "flex_adjust"
        self._drag_anchor = None
        self._drag_free_orientation = False
        self._bond_drag_start = None
        if self._overlays_ready:
            self._optimize_zone.hide_zone()
            self._preview_bond_item.update_flex_curve(p0, p1, curve_1, curve_2)

    def _update_flex_bond_preview(
        self,
        scene_pos: QPointF,
        modifiers: Qt.KeyboardModifiers,
    ) -> None:
        """Actualiza previsualización durante ajuste de enlace FLEX."""
        pending = self._flex_bond_pending
        if not pending:
            return
        p0 = pending["p0"]
        p1 = pending["p1"]
        curve_1, curve_2 = self._flex_curve_from_pointer(p0, p1, scene_pos, modifiers)
        pending["curve_1"] = float(curve_1)
        pending["curve_2"] = float(curve_2)
        if self._overlays_ready:
            self._preview_bond_item.update_flex_curve(p0, p1, curve_1, curve_2)

    def _commit_flex_bond(self) -> None:
        """Confirma creación del enlace FLEX con la curvatura previsualizada."""
        pending = self._flex_bond_pending
        if not pending:
            self._cancel_drag()
            return

        p0 = QPointF(pending["p0"])
        p1 = QPointF(pending["p1"])
        anchor_id = pending["anchor_id"]
        target_id = pending["target_id"]
        use_default = bool(pending["use_default"])
        curve_1 = float(pending["curve_1"])
        curve_2 = float(pending["curve_2"])

        order = 1
        style = BondStyle.FLEX
        stereo = BondStereo.NONE
        is_aromatic = False

        if anchor_id is None:
            self.undo_stack.beginMacro("Add flexible bond")
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
                    flex_curve_1=curve_1,
                    flex_curve_2=curve_2,
                )
                target_atom = self.model.get_atom(target_id)
                final_p1 = QPointF(target_atom.x, target_atom.y)
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
                    flex_curve_1=curve_1,
                    flex_curve_2=curve_2,
                )
                self.undo_stack.push(cmd)
                final_p1 = QPointF(p1)
            self.undo_stack.endMacro()
            self._update_bond_angle_state(p0, final_p1, use_default, anchor_id)
            self._cancel_drag()
            return

        if target_id is not None:
            self._create_or_update_bond(
                anchor_id,
                target_id,
                order,
                style,
                stereo,
                is_aromatic,
                flex_curve_1=curve_1,
                flex_curve_2=curve_2,
            )
            target_atom = self.model.get_atom(target_id)
            final_p1 = QPointF(target_atom.x, target_atom.y)
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
                flex_curve_1=curve_1,
                flex_curve_2=curve_2,
            )
            self.undo_stack.push(cmd)
            final_p1 = QPointF(p1)
        self._update_bond_angle_state(p0, final_p1, use_default, anchor_id)
        self._cancel_drag()

    def _finalize_bond(self, modifiers: Qt.KeyboardModifiers) -> None:
        """Método auxiliar para  finalize bond.

        Args:
            modifiers: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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

        if style == BondStyle.FLEX:
            self._begin_flex_adjust_mode(
                p0,
                final_p1,
                anchor_id,
                target_id,
                use_default,
                modifiers,
            )
            return

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
        """Método auxiliar para  compute bond endpoint.

        Args:
            anchor: Descripción del parámetro.
            scene_pos: Descripción del parámetro.
            modifiers: Descripción del parámetro.
            bond_order: Descripción del parámetro.
            is_aromatic: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
            allow_length_boost=False,
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
        """Método auxiliar para  begin place ring.

        Args:
            anchor_type: Descripción del parámetro.
            anchor_id: Descripción del parámetro.
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  update ring preview.

        Args:
            scene_pos: Descripción del parámetro.
            modifiers: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        vertices = self._compute_ring_vertices(scene_pos, modifiers)
        self._ring_last_vertices = vertices
        self._preview_ring_item.update_polygon(vertices)

    def _finalize_ring(self, modifiers: Qt.KeyboardModifiers) -> None:
        """Método auxiliar para  finalize ring.

        Args:
            modifiers: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  compute ring vertices.

        Args:
            scene_pos: Descripción del parámetro.
            modifiers: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  regular polygon from edge.

        Args:
            p1: Descripción del parámetro.
            p2: Descripción del parámetro.
            ring_size: Descripción del parámetro.
            direction: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  begin place chain.

        Args:
            anchor_id: Descripción del parámetro.
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  update chain preview.

        Args:
            scene_pos: Descripción del parámetro.
            modifiers: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
            first_angle = raw_angle
        else:
            cursor_angle = raw_angle
            if anchor_id is None:
                # En cadena libre, usar el eje del cursor para permitir un zig-zag
                # "recto" (estilo ChemDraw) en cualquier dirección.
                step_deg = float(self.state.angle_step_deg)
                first_angle = snap_angle_deg(cursor_angle, step_deg) if step_deg > 0 and not use_shift else cursor_angle
            else:
                if use_optimize:
                    cursor_angle = self._select_preferred_angle(anchor_id, cursor_angle, 1, False)
                first_angle, _ = self._pick_bond_direction_deg(
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
        # La herramienta de cadena siempre debe zigzaguear (si ángulos fijos están activos),
        # incluso cuando el ancla está en centros sp2/aromáticos.
        zigzag = not use_alt and self.state.fixed_angles
        sp3_angle = self._sp3_display_angle_deg()
        zigzag_turn = 180.0 - sp3_angle
        turn_pref = ((raw_angle - first_angle + 180.0) % 360.0) - 180.0
        zigzag_sign = 1.0 if turn_pref >= 0.0 else -1.0
        first_step_angle = first_angle
        second_step_angle = first_angle + zigzag_sign * zigzag_turn
        if zigzag and anchor_id is None:
            # Cadena libre: alternar simétricamente alrededor del eje pedido
            # por el cursor para formar una línea zig-zag recta.
            half_turn = zigzag_turn * 0.5
            axis_angle = first_angle
            first_step_angle = axis_angle - zigzag_sign * half_turn
            second_step_angle = axis_angle + zigzag_sign * half_turn
        current = p0
        for i in range(1, n + 1):
            if zigzag and i % 2 == 0:
                step_angle = second_step_angle
            else:
                step_angle = first_step_angle
            next_point = endpoint_from_angle_len(current, step_angle, self.state.bond_length)
            points.append(next_point)
            current = next_point

        self._chain_last_points = points
        self._preview_chain_item.update_polyline(points)
        self._preview_chain_label.update_label(str(n), points[-1] + QPointF(0, -10))

    def _finalize_chain(self, modifiers: Qt.KeyboardModifiers) -> None:
        """Método auxiliar para  finalize chain.

        Args:
            modifiers: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  handle hotkeys.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        key = event.key()
        modifiers = event.modifiers()

        if (
            key == Qt.Key.Key_0
            and (modifiers & Qt.KeyboardModifier.AltModifier)
            and not (
                modifiers
                & (
                    Qt.KeyboardModifier.ControlModifier
                    | Qt.KeyboardModifier.MetaModifier
                    | Qt.KeyboardModifier.ShiftModifier
                )
            )
        ):
            if self._is_rotating_3d:
                self._finalize_3d_rotation_drag()
            self._reset_3d_rotation_perspective()
            return True

        if Qt.Key.Key_0 <= key <= Qt.Key.Key_9:
            digit = key - Qt.Key.Key_0
            if digit == 0:
                return False
            if self.hovered_atom_id is not None:
                if modifiers & Qt.KeyboardModifier.ShiftModifier:
                    if digit >= 3:
                        self._create_ring_hotkey(self.hovered_atom_id, digit)
                        return True
                else:
                    self._create_chain_hotkey(self.hovered_atom_id, digit)
                    return True
            if self.hovered_bond_id is not None:
                if modifiers & Qt.KeyboardModifier.ShiftModifier:
                    if digit >= 3:
                        self._create_ring_from_bond_hotkey(self.hovered_bond_id, digit)
                        return True
                else:
                    self._set_bond_order_hotkey(self.hovered_bond_id, digit)
                    return True
        return False

    def _create_chain_hotkey(self, anchor_id: int, n: int) -> None:
        """Método auxiliar para  create chain hotkey.

        Args:
            anchor_id: Descripción del parámetro.
            n: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  create ring hotkey.

        Args:
            anchor_id: Descripción del parámetro.
            ring_size: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        anchor = self.model.get_atom(anchor_id)
        p0 = QPointF(anchor.x, anchor.y)
        angles = self._get_anchor_bond_angles_deg(anchor_id)
        theta = choose_optimal_direction(angles)
        p1 = endpoint_from_angle_len(p0, theta, self.state.bond_length)
        direction = bond_side(p0, p1, self._last_scene_pos)
        vertices = self._regular_polygon_from_edge(p0, p1, ring_size, direction)
        self._commit_ring_vertices(vertices, anchor_type="atom", anchor_id=anchor_id)

    def _create_ring_from_bond_hotkey(self, bond_id: int, ring_size: int) -> None:
        """Método auxiliar para  create ring from bond hotkey.

        Args:
            bond_id: Descripción del parámetro.
            ring_size: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  commit ring vertices.

        Args:
            vertices: Descripción del parámetro.
            anchor_type: Descripción del parámetro.
            anchor_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  set bond order hotkey.

        Args:
            bond_id: Descripción del parámetro.
            order: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  kekulize aromatic bonds.

        Args:
            seed_atoms: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        for bond in aromatic_bonds:
            if component_atoms is not None and (
                bond.a1_id not in component_atoms or bond.a2_id not in component_atoms
            ):
                continue
            self.update_bond_item(bond.id)
        self._assign_ring_ids_for_aromatic_cycles()
        self._refresh_aromatic_ring_contexts()

    def _assign_ring_ids_for_aromatic_cycles(self) -> None:
        """Método auxiliar para  assign ring ids for aromatic cycles.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        adjacency: dict[int, list[tuple[int, int]]] = {}
        for bond in self.model.bonds.values():
            if not bond.is_aromatic:
                continue
            adjacency.setdefault(bond.a1_id, []).append((bond.a2_id, bond.id))
            adjacency.setdefault(bond.a2_id, []).append((bond.a1_id, bond.id))

        for bond in self.model.bonds.values():
            if not bond.is_aromatic or bond.ring_id is not None:
                continue
            cycle = self._find_cycle_for_bond(bond, adjacency)
            if not cycle:
                continue
            ring_id = None
            for bond_id in cycle["bond_ids"]:
                existing = self.model.get_bond(bond_id)
                if existing.ring_id is not None:
                    ring_id = existing.ring_id
                    break
            if ring_id is None:
                ring_id = self.allocate_ring_id()
            for bond_id in cycle["bond_ids"]:
                b = self.model.get_bond(bond_id)
                if b.ring_id is None:
                    b.ring_id = ring_id
                    self.update_bond_item(bond_id)
            center = self._center_for_atoms(cycle["atom_ids"])
            if center is not None:
                self.register_ring_center(ring_id, (center.x(), center.y()))

    def _find_cycle_for_bond(
        self, bond: Bond, adjacency: dict[int, list[tuple[int, int]]]
    ) -> Optional[dict[str, list[int]]]:
        """Método auxiliar para  find cycle for bond.

        Args:
            bond: Descripción del parámetro.
            adjacency: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        a1 = bond.a1_id
        a2 = bond.a2_id
        if a1 not in adjacency or a2 not in adjacency:
            return None
        from collections import deque

        queue = deque([(a1, [a1], [])])
        while queue:
            node, path_atoms, path_bonds = queue.popleft()
            if len(path_atoms) > 8:
                continue
            for neighbor, bond_id in adjacency.get(node, []):
                if bond_id == bond.id:
                    continue
                if neighbor == a2:
                    if len(path_atoms) >= 4:
                        atom_ids = path_atoms + [a2]
                        bond_ids = path_bonds + [bond_id, bond.id]
                        return {"atom_ids": atom_ids, "bond_ids": bond_ids}
                    continue
                if neighbor in path_atoms:
                    continue
                queue.append((neighbor, path_atoms + [neighbor], path_bonds + [bond_id]))
        return None

    def _build_aromatic_edges(
        self,
        vertex_defs: List[Tuple[Optional[int], float, float]],
        ring_size: int,
    ) -> List[Tuple[int, int, int, BondStyle, BondStereo, bool]]:
        """Método auxiliar para  build aromatic edges.

        Args:
            vertex_defs: Descripción del parámetro.
            ring_size: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
            """Método auxiliar para atom has double.

            Args:
                atom_id: Descripción del parámetro.
                other_id: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
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
            """Método auxiliar para score phase.

            Args:
                phase: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
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
        """Método auxiliar para  begin drag.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if (
            not self._selected_atom_ids_for_transform()
            and not self._selected_text_items()
            and not self._selected_arrow_items()
            and not self._selected_bracket_items()
        ):
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

        # Capture arrows
        self._drag_start_arrow_positions = {}
        for item in self._selected_arrow_items():
            self._drag_start_arrow_positions[item] = (item.start_point(), item.end_point())

        # Capture brackets
        self._drag_start_bracket_rects = {}
        for item in self._selected_bracket_items():
            self._drag_start_bracket_rects[item] = item.base_rect()
            
        self._drag_has_moved = False

    def _is_on_paper(self, x: float, y: float) -> bool:
        """Método auxiliar para  is on paper.

        Args:
            x: Descripción del parámetro.
            y: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        return (
            PAPER_MARGIN <= x <= self.paper_width - PAPER_MARGIN
            and PAPER_MARGIN <= y <= self.paper_height - PAPER_MARGIN
        )

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
