"""
Lienzo principal de Chemuson (QGraphicsView/QGraphicsScene).

Gestiona el dibujo de moléculas, herramientas de edición, anotaciones,
persistencia visual y renderización para exportaciones.
"""
from __future__ import annotations

import base64
import itertools
import json
import math
import os
from copy import deepcopy
from typing import Dict, Iterable, List, Optional, Tuple

from PyQt6.QtWidgets import (
    QApplication,
    QMenu,
    QGraphicsView,
    QGraphicsScene,
    QGraphicsRectItem,
    QGraphicsPathItem,
    QGraphicsDropShadowEffect,
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
    QFontMetrics,
    QTextCharFormat,
    QTextCursor,
    QTextDocument,
    QTextOption,
)
from PyQt6.QtCore import Qt, QPoint, QPointF, QRectF, QRect, QSize, QBuffer, QMimeData, QTimer, pyqtSignal

try:
    from PyQt6.QtSvg import QSvgGenerator
except Exception:  # Optional Qt module at runtime
    QSvgGenerator = None

from chemuson.core.model import (
    Bond,
    BondStyle,
    BondStereo,
    ChemState,
    MolGraph,
    bond_is_structural,
    normalize_opacity,
    normalize_optional_opacity,
)
from chemuson.gui.items import (
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
    EnergyDiagramItem,
    OrbitalAnnotationItem,
    ImageAnnotationItem,
    BaseItem,
    ABBREVIATION_LABELS,
)
from chemuson.gui.composite_diagram_item import CompositeDiagramItem
from chemuson.gui.diagram_presets import build_semantic_diagram_from_preset
from chemuson.gui.diagram_models import SemanticDiagram
from chemuson.gui.energy_diagrams import (
    DEFAULT_ENERGY_DIAGRAM_KIND,
    build_atomic_subshell_diagram,
    build_diatomic_mo_diagram,
    build_ligand_field_diagram,
    default_energy_label,
    default_energy_label_side,
    energy_diagram_default_size,
    energy_diagram_display_name,
    energy_diagram_kind_from_tool_id,
    normalize_energy_occupancies,
)
from chemuson.gui.style import CHEMDOODLE_LIKE, DrawingStyle
from chemuson.gui.numbering import NumberedStructure, compute_atom_numbers, compute_structure_numbers
from chemuson.gui.orbitals import (
    DEFAULT_ORBITAL_KIND,
    orbital_canvas_extent,
    orbital_kind_from_tool_id,
)
from chemuson.gui.templates import build_haworth_template, build_chair_template
from chemuson.gui.geom import (
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
from chemuson.gui.commands import (
    AddAtomCommand,
    AddBondCommand,
    AddRingCommand,
    AddChainCommand,
    AddPlateItemCommand,
    MovePlateItemsCommand,
    ChangeAtomCommand,
    AddArrowCommand,
    AddBracketCommand,
    AddTextItemCommand,
    AddImageItemCommand,
    AddEnergyDiagramItemCommand,
    AddCompositeDiagramItemCommand,
    AddOrbitalItemCommand,
    AddWavyAnchorCommand,
    ChangeBondCommand,
    ChangeBondLengthCommand,
    ChangeBondStrokeCommand,
    ChangeArrowStrokeCommand,
    ChangeBracketStrokeCommand,
    ChangeBondColorCommand,
    ChangeDoubleBondOrientationCommand,
    ChangeChargeCommand,
    ChangeNoImplicitCommand,
    ChangeAtomLabelScaleCommand,
    ChangeCoordinationSphereStyleCommand,
    ChangeCanvasOpacityCommand,
    StyleOrbitalItemsCommand,
    SetCoordinationCenterCommand,
    DeleteSelectionCommand,
    EditSemanticDiagramCommand,
    MoveAtomsCommand,
    MoveTextItemsCommand,
    MoveArrowItemsCommand,
    MoveBracketItemsCommand,
    TransformImageItemsCommand,
    TransformEnergyDiagramItemsCommand,
    TransformOrbitalItemsCommand,
    ScaleTextItemsCommand,
    FormatTextItemsCommand,
    ScaleArrowItemsCommand,
    ScaleBracketItemsCommand,
    ConfigureEnergyDiagramItemsCommand,
    AddPlateItemCommand,
    MovePlateItemsCommand,
)
from chemuson.gui.dialogs import (
    AtomLabelDialog,
    TrackballRotationDialog,
    SelectionScaleDialog,
    TLCInsertDialog,
    GelInsertDialog,
)
from chemuson.chemname.molview import MolView
from chemuson.chemname.rings import find_rings_simple, ring_bonds
from chemuson.gui.plate_items import TLCPlateItem, GelElectrophoresisItem, TLCSpotItem, GelBandItem, PlateItem
from chemuson.chemio.rdkit_io import (
    molgraph_to_molfile,
    molgraph_to_smiles,
    molfile_to_molgraph,
    smiles_to_molgraph,
    kekulize_display_orders,
)
from chemuson.chemname import iupac_name, NameOptions


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
BRANCH_ROTATION_STEP_DEG = 60.0
BRANCH_ROTATION_NOOP_TOLERANCE_DEG = 0.15
FRAGMENT_ROTATION_STEP_DEG = 30.0
CLIPBOARD_RENDER_SCALE = 5.0
CLIPBOARD_LARGE_SELECTION_RENDER_SCALE = 2.5
CLIPBOARD_LARGE_SELECTION_ATOM_THRESHOLD = 18
CLIPBOARD_LARGE_SELECTION_BOND_THRESHOLD = 20
PASTE_IMAGE_OFFSET_PX = 20.0
PASTE_IMAGE_MAX_PAPER_FRACTION = 0.72
TRACKBALL_ROTATION_DEG_PER_PIXEL = 1.0
TRACKBALL_MAX_TILT_DEG = 60.0
TRACKBALL_REFERENCE_MATCH_TOLERANCE_PX = 1.5
SUPPORTED_IMAGE_FILE_MIME_TYPES = {
    ".png": "image/png",
    ".jpg": "image/jpeg",
    ".jpeg": "image/jpeg",
    ".svg": "image/svg+xml",
}
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
SIMPLE_HYDROGEN_GROUP_LABELS: dict[str, tuple[str, int]] = {
    "NH2": ("N", 2),
    "H2N": ("N", 2),
    "NH": ("N", 1),
    "HN": ("N", 1),
    "OH": ("O", 1),
    "HO": ("O", 1),
    "SH": ("S", 1),
    "HS": ("S", 1),
}
CANONICAL_SIMPLE_HYDROGEN_GROUP_LABELS: dict[tuple[str, int], str] = {
    ("N", 2): "NH2",
    ("N", 1): "NH",
    ("O", 1): "OH",
    ("S", 1): "SH",
}
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
IMPLICIT_H_OVERLAY_ANCHOR_ROLE = 5001
IMPLICIT_H_OVERLAY_ANGLE_ROLE = 5002
ITEM_OPACITY_ROLE = 5003
PAPER_ITEM_ROLE = 7001
ITEM_OPACITY_UNSET = object()

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


