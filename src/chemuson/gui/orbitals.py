"""Renderer vectorial y parametrizable para la paleta orbital."""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from enum import Enum
import json
import math
from pathlib import Path
from typing import Any

from PyQt6.QtCore import QPointF, QRectF, QSize, Qt
from PyQt6.QtGui import (
    QBrush,
    QColor,
    QIcon,
    QImage,
    QLinearGradient,
    QPainter,
    QPainterPath,
    QPen,
    QPixmap,
    QRadialGradient,
    QTransform,
)


ORBITAL_ICON_SIZE = 34
DEFAULT_ORBITAL_KIND = "p_shaded"

_STROKE_COLOR = QColor("#111111")
_SOLID_FILL = QColor("#111111")
_PALETTE_BG = QColor("#F2F2F2")
_CELL_BG = QColor("#FFFFFF")
_CELL_LINE = QColor("#D8D8D8")
_DEFAULT_LIGHT_DIR = "northwest"


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


ORBITAL_PRESET_CONFIG_PATH = _repo_root() / "config" / "orbitals_presets.json"

_GRADIENT_STOP_PRESETS: dict[str, tuple[tuple[float, str], ...]] = {
    "linear": (
        (0.0, "#EDEDED"),
        (0.34, "#D6D6D6"),
        (0.72, "#9B9B9B"),
        (1.0, "#6A6A6A"),
    ),
    "radial": (
        (0.0, "#F4F4F4"),
        (0.30, "#DDDDDD"),
        (0.72, "#A1A1A1"),
        (1.0, "#6E6E6E"),
    ),
    "elliptical": (
        (0.0, "#DEDEDE"),
        (0.30, "#C1C1C1"),
        (0.72, "#8F8F8F"),
        (1.0, "#666666"),
    ),
    "negative": (
        (0.0, "#BCBCBC"),
        (0.30, "#8D8D8D"),
        (0.72, "#4B4B4B"),
        (1.0, "#171717"),
    ),
}

_LIGHT_DIRECTION_VECTORS: dict[str, tuple[float, float]] = {
    "north": (0.0, -1.0),
    "south": (0.0, 1.0),
    "east": (1.0, 0.0),
    "west": (-1.0, 0.0),
    "northeast": (0.7071067812, -0.7071067812),
    "northwest": (-0.7071067812, -0.7071067812),
    "southeast": (0.7071067812, 0.7071067812),
    "southwest": (-0.7071067812, 0.7071067812),
    "center": (0.0, 0.0),
}

ORBITAL_PARAMETER_DESCRIPTIONS: dict[str, str] = {
    "radius_x": "Semieje horizontal de la esfera u orbital s.",
    "radius_y": "Semieje vertical de la esfera u orbital s.",
    "visual_padding": "Padding óptico alrededor del glifo al encajar en la palette.",
    "anchor_bias_x": "Desplazamiento horizontal del centro visual del glifo.",
    "anchor_bias_y": "Desplazamiento vertical del centro visual del glifo.",
    "gradient_mode": "Modo base de sombreado: linear, radial o elliptical.",
    "light_dir": "Dirección global de la luz del sombreado.",
    "width": "Ancho total del lóbulo experimental.",
    "height": "Altura total del lóbulo experimental.",
    "tip_roundness": "Qué tan redondeada o aguda se ve la punta del lóbulo.",
    "shoulder_width": "Apertura del hombro superior del lóbulo.",
    "waist_width": "Ancho del cuello o cintura cerca del nodo.",
    "node_gap": "Separación controlada entre lóbulos enfrentados en el nodo.",
    "vertical_offset": "Desplazamiento vertical conjunto de la familia orbital.",
    "geometry_scale_x": "Escala horizontal aplicada a la geometría base antes del estilo.",
    "geometry_scale_y": "Escala vertical aplicada a la geometría base antes del estilo.",
    "shoulder_height_ratio": "Altura relativa del hombro del lóbulo experimental.",
    "belly_width": "Qué tan lleno se ve el vientre o zona media del lóbulo.",
    "belly_height_ratio": "Altura relativa del vientre o zona media del lóbulo.",
    "neck_width": "Ancho del cuello del lóbulo experimental.",
    "neck_height": "Altura relativa del cuello experimental.",
    "cusp_depth": "Profundidad del cierre en el nodo para formas experimentales.",
    "offset_x": "Desplazamiento horizontal local del lóbulo.",
    "offset_y": "Desplazamiento vertical local del lóbulo.",
    "rotation_deg": "Rotación del lóbulo en grados.",
    "orientation": "Orientación base del lóbulo: up, down, left o right.",
    "mirror": "Espeja el lóbulo respecto al eje vertical antes de colocarlo.",
    "outline_scale_x": "Escala horizontal aplicada a la geometría canónica base.",
    "outline_scale_y": "Escala vertical aplicada a la geometría canónica base.",
    "phase": "Fase visual de la parte: positive, negative o neutral.",
    "lobe_width": "Ancho base del lóbulo principal de la familia.",
    "lobe_height": "Altura base del lóbulo principal de la familia.",
    "vertical_width": "Ancho de los lóbulos verticales del trébol.",
    "vertical_height": "Altura de los lóbulos verticales del trébol.",
    "horizontal_width": "Ancho de los lóbulos horizontales del trébol.",
    "horizontal_height": "Altura de los lóbulos horizontales del trébol.",
    "center_roundness": "Redondeo del punto de anclaje central de cada lóbulo.",
    "lobe_gap": "Separación entre lóbulos alrededor del nodo central.",
    "side_rotation_deg": "Rotación de los lóbulos laterales del conjunto sp3.",
    "bottom_rotation_deg": "Rotación del lóbulo inferior del conjunto sp3.",
    "side_scale": "Escala aplicada a los lóbulos laterales del conjunto sp3.",
    "bottom_scale": "Escala aplicada al lóbulo inferior del conjunto sp3.",
    "pair_rotation_deg": "Ángulo entre el lóbulo superior y cada lóbulo lateral del conjunto sp2.",
    "planar_rotation_deg": "Rotación global del conjunto sp2 dentro del plano.",
    "top_scale": "Escala aplicada al lóbulo superior del conjunto sp2.",
    "pair_scale": "Escala aplicada al par lateral del conjunto sp2.",
    "lateral_offset": "Apertura lateral de la nube enlazante del orbital pi.",
    "axial_width": "Ancho de los lóbulos axiales del dz2.",
    "axial_height": "Altura de los lóbulos axiales del dz2.",
    "axial_tip_roundness": "Redondez de la punta axial.",
    "axial_shoulder_width": "Apertura del hombro del lóbulo axial del dz2.",
    "axial_waist_width": "Ancho de la cintura axial cerca del plano nodal.",
    "axial_center_roundness": "Redondeo local del extremo axial hacia el centro.",
    "torus_outer_width": "Ancho exterior del toroide.",
    "torus_outer_height": "Altura exterior del toroide.",
    "torus_inner_width_ratio": "Tamaño relativo del hueco interior del toroide en X.",
    "torus_inner_height_ratio": "Tamaño relativo del hueco interior del toroide en Y.",
    "torus_offset_y": "Desplazamiento vertical del toroide.",
    "major_lobe": "Parámetros del lóbulo grande o principal.",
    "minor_lobe": "Parámetros del lóbulo pequeño o secundario.",
    "upper_lobe": "Parámetros del lóbulo superior.",
    "lower_lobe": "Parámetros del lóbulo inferior.",
    "major_offset_x": "Desplazamiento horizontal del lóbulo principal.",
    "minor_offset_x": "Desplazamiento horizontal del lóbulo secundario.",
    "major_offset_y": "Desplazamiento vertical del lóbulo principal.",
    "minor_offset_y": "Desplazamiento vertical del lóbulo secundario.",
    "bridge_width": "Ancho de la concavidad central de la nube enlazante pi.",
    "bridge_height": "Profundidad vertical de la concavidad central de la nube pi.",
    "cap_width": "Ancho de la corona exterior que redondea la nube pi.",
    "cap_height": "Altura extra de la corona exterior de la nube pi.",
    "torus_gradient_mode": "Modo de sombreado específico del toroide dentro del dz2.",
}


class OrbitalType(str, Enum):
    """Tipos orbitales soportados por la UI."""

    S = "s"
    P = "p"
    D = "d"
    DZ2 = "dz2"
    F = "f"
    FZ3 = "fz3"
    FXZ2 = "fxz2"
    FZX2Y2 = "fzx2y2"
    FXYZ = "fxyz"
    SP3 = "sp3"
    SP2 = "sp2"
    SP_LOBE = "sp_lobe"
    SIGMA_BONDING = "sigma_bonding"
    PI_BONDING = "pi_bonding"
    TORUS = "torus"


class OrbitalStyle(str, Enum):
    """Variantes visuales para un mismo glifo."""

    OUTLINE = "outline"
    SHADED = "shaded"
    SOLID = "solid"


@dataclass(frozen=True)
class GlyphLayer:
    """Capa individual de un glifo con su modo de pintura."""

    path: QPainterPath
    paint: str
    name: str = ""
    gradient: str = "linear"
    stroke: bool = False
    phase: str | None = None
    light_dir: str | None = None


@dataclass(frozen=True)
class GlyphDefinition:
    """Glifo ya resuelto a capas pintables."""

    id: str
    paths_outline: tuple[GlyphLayer, ...]
    paths_shaded: tuple[GlyphLayer, ...]
    paths_solid: tuple[GlyphLayer, ...]
    default_bounds: QRectF
    anchor_center: QPointF
    anchor_bias: QPointF
    visual_padding: float


@dataclass(frozen=True)
class SphereParams:
    """Parámetros geométricos de un orbital s.

    `radius_x` y `radius_y` controlan el radio visual.
    `gradient_mode` y `light_dir` controlan solo el sombreado, no la silueta.
    """

    radius_x: float
    radius_y: float
    visual_padding: float
    anchor_bias_x: float
    anchor_bias_y: float
    gradient_mode: str = "radial"
    light_dir: str = _DEFAULT_LIGHT_DIR


@dataclass(frozen=True)
class CanonicalLobeParams:
    """Parámetros seguros de un lóbulo canónico y convexo."""

    width: float
    height: float
    tip_roundness: float
    shoulder_width: float
    waist_width: float
    belly_width: float = 0.0
    belly_height_ratio: float = 0.0


@dataclass(frozen=True)
class TeardropParams:
    """Parámetros libres para lóbulos experimentales f/fz3."""

    width: float
    height: float
    tip_roundness: float
    shoulder_width: float
    shoulder_height_ratio: float
    belly_width: float
    belly_height_ratio: float
    neck_width: float
    neck_height: float
    cusp_depth: float
    offset_x: float = 0.0
    offset_y: float = 0.0
    rotation_deg: float = 0.0
    orientation: str = "up"
    mirror: bool = False
    outline_scale_x: float = 1.0
    outline_scale_y: float = 1.0
    gradient_mode: str = "linear"
    light_dir: str = _DEFAULT_LIGHT_DIR
    phase: str = "positive"


@dataclass(frozen=True)
class TorusParams:
    """Parámetros de un toroide editable.

    `torus_inner_width_ratio` y `torus_inner_height_ratio` controlan el hueco interior.
    """

    torus_outer_width: float
    torus_outer_height: float
    torus_inner_width_ratio: float
    torus_inner_height_ratio: float
    torus_offset_y: float = 0.0
    visual_padding: float = 6.0
    anchor_bias_x: float = 0.0
    anchor_bias_y: float = 0.0
    gradient_mode: str = "elliptical"
    light_dir: str = _DEFAULT_LIGHT_DIR
    phase: str = "neutral"


@dataclass(frozen=True)
class POrbitalParams:
    """Parámetros seguros para un orbital p canónico."""

    lobe_width: float
    lobe_height: float
    waist_width: float
    shoulder_width: float
    tip_roundness: float
    node_gap: float
    vertical_offset: float
    outline_scale_x: float
    outline_scale_y: float
    visual_padding: float
    anchor_bias_x: float
    anchor_bias_y: float
    gradient_mode: str = "linear"
    light_dir: str = _DEFAULT_LIGHT_DIR


@dataclass(frozen=True)
class CloverParams:
    """Parámetros seguros para un orbital d tipo trébol."""

    vertical_width: float
    vertical_height: float
    horizontal_width: float
    horizontal_height: float
    tip_roundness: float
    shoulder_width: float
    waist_width: float
    center_roundness: float
    lobe_gap: float
    rotation_deg: float
    visual_padding: float
    anchor_bias_x: float
    anchor_bias_y: float
    gradient_mode: str = "linear"
    light_dir: str = _DEFAULT_LIGHT_DIR


@dataclass(frozen=True)
class Dz2Params:
    """Parámetros del dz2 con lóbulos axiales y toroide ecuatorial."""

    axial_width: float
    axial_height: float
    axial_tip_roundness: float
    axial_shoulder_width: float
    axial_waist_width: float
    axial_center_roundness: float
    torus_outer_width: float
    torus_outer_height: float
    torus_inner_width_ratio: float
    torus_inner_height_ratio: float
    torus_offset_y: float
    visual_padding: float
    anchor_bias_x: float
    anchor_bias_y: float
    gradient_mode: str = "linear"
    torus_gradient_mode: str = "elliptical"
    light_dir: str = _DEFAULT_LIGHT_DIR


@dataclass(frozen=True)
class HybridOrbitalParams:
    """Parámetros seguros para híbridos tipo sp_lobe, sp3 y sigma."""

    major_lobe: CanonicalLobeParams
    minor_lobe: CanonicalLobeParams
    major_offset_x: float = 0.0
    minor_offset_x: float = 0.0
    major_offset_y: float = 0.0
    minor_offset_y: float = 0.0
    center_roundness: float = 0.0
    node_gap: float = 0.0
    visual_padding: float = 6.0
    anchor_bias_x: float = 0.0
    anchor_bias_y: float = 0.0
    gradient_mode: str = "linear"
    minor_gradient_mode: str = "linear"
    light_dir: str = _DEFAULT_LIGHT_DIR


@dataclass(frozen=True)
class Sp3Params:
    """Parámetros del conjunto de 4 lóbulos para un sp3 canónico."""

    lobe: CanonicalLobeParams
    center_roundness: float = 0.0
    side_rotation_deg: float = 126.0
    bottom_rotation_deg: float = 180.0
    side_scale: float = 0.96
    bottom_scale: float = 0.88
    visual_padding: float = 7.0
    anchor_bias_x: float = 0.0
    anchor_bias_y: float = 2.0
    gradient_mode: str = "linear"
    light_dir: str = _DEFAULT_LIGHT_DIR


@dataclass(frozen=True)
class Sp2Params:
    """Parámetros del conjunto trigonal plano para un sp2 canónico."""

    lobe: CanonicalLobeParams
    center_roundness: float = 0.0
    pair_rotation_deg: float = 120.0
    planar_rotation_deg: float = 0.0
    top_scale: float = 1.0
    pair_scale: float = 1.0
    visual_padding: float = 7.0
    anchor_bias_x: float = 0.0
    anchor_bias_y: float = 1.5
    gradient_mode: str = "linear"
    light_dir: str = _DEFAULT_LIGHT_DIR


@dataclass(frozen=True)
class PiBondingParams:
    """Parámetros seguros para un orbital pi enlazante."""

    upper_lobe: CanonicalLobeParams
    lower_lobe: CanonicalLobeParams
    node_gap: float = 0.0
    vertical_offset: float = 0.0
    center_roundness: float = 0.0
    lateral_offset: float = 12.0
    bridge_width: float = 0.0
    bridge_height: float = 0.0
    cap_width: float = 0.0
    cap_height: float = 0.0
    visual_padding: float = 6.0
    anchor_bias_x: float = 0.0
    anchor_bias_y: float = 0.0
    gradient_mode: str = "linear"
    light_dir: str = _DEFAULT_LIGHT_DIR


@dataclass(frozen=True)
class FOrbitalParams:
    """Parámetros flexibles para orbitales f/fz3.

    La cantidad de lóbulos se define por la lista `lobes`, no por código rígido.
    """

    lobes: tuple[TeardropParams, ...]
    torus: TorusParams | None = None
    visual_padding: float = 8.0
    anchor_bias_x: float = 0.0
    anchor_bias_y: float = 0.0


@dataclass(frozen=True)
class OrbitalFamilyPreset:
    """Contenedor general de una familia orbital editable."""

    family: str
    builder: str
    params: SphereParams | POrbitalParams | CloverParams | Dz2Params | FOrbitalParams | HybridOrbitalParams | Sp3Params | Sp2Params | PiBondingParams | TorusParams


@dataclass(frozen=True)
class GeometryPart:
    """Parte geométrica reusable antes de estilizar."""

    name: str
    path: QPainterPath
    phase: str = "positive"
    gradient_mode: str = "linear"
    light_dir: str = _DEFAULT_LIGHT_DIR


@dataclass(frozen=True)
class FamilyBuildResult:
    """Resultado base de un builder antes de convertirlo en `GlyphDefinition`."""

    family: str
    parts: tuple[GeometryPart, ...]
    outline_parts: tuple[str, ...]
    shaded_fill_parts: tuple[str, ...]
    shaded_outline_parts: tuple[str, ...]
    solid_fill_parts: tuple[str, ...]
    solid_outline_parts: tuple[str, ...]
    visual_padding: float
    anchor_bias_x: float
    anchor_bias_y: float


@dataclass(frozen=True)
class OrbitalGlyphSpec:
    """Especificación lógica de una variante orbital."""

    kind: str
    orbital_type: OrbitalType
    glyph_id: str
    style: OrbitalStyle
    label: str
    canvas_extent_scale: float = 0.92
    stroke_shaded_lobes: bool = True


@dataclass(frozen=True)
class OrbitalPaletteModel:
    """Grid fijo de la paleta orbital."""

    cells: tuple[tuple[str | None, ...], ...]
    cell_width: int = 70
    cell_height: int = 56
    inner_padding: int = 6

    @property
    def columns(self) -> int:
        return max((len(row) for row in self.cells), default=0)

    @property
    def rows(self) -> int:
        return len(self.cells)

    @property
    def entries(self) -> tuple[str, ...]:
        return tuple(kind for row in self.cells for kind in row if kind)

    def kind_at(self, row: int, column: int) -> str | None:
        if row < 0 or row >= len(self.cells):
            return None
        current_row = self.cells[row]
        if column < 0 or column >= len(current_row):
            return None
        return current_row[column]

    @property
    def image_size(self) -> QSize:
        return QSize(self.columns * self.cell_width, self.rows * self.cell_height)

    @property
    def icon_size(self) -> QSize:
        side = min(self.cell_width, self.cell_height) - 10
        return QSize(max(18, side), max(18, side))


_TYPE_LABELS = {
    OrbitalType.S: "Orbital s",
    OrbitalType.P: "Orbital p",
    OrbitalType.D: "Orbital d",
    OrbitalType.DZ2: "Orbital dz2",
    OrbitalType.F: "Orbital f",
    OrbitalType.FZ3: "Orbital fz3",
    OrbitalType.FXZ2: "Orbital fxz2",
    OrbitalType.FZX2Y2: "Orbital fz(x2-y2)",
    OrbitalType.FXYZ: "Orbital fxyz",
    OrbitalType.SP3: "Orbital sp3",
    OrbitalType.SP2: "Orbital sp2",
    OrbitalType.SP_LOBE: "Lóbulo sp",
    OrbitalType.SIGMA_BONDING: "Sigma enlazante",
    OrbitalType.PI_BONDING: "Pi enlazante",
    OrbitalType.TORUS: "Toroide",
}

_STYLE_LABELS = {
    OrbitalStyle.OUTLINE: "contorno",
    OrbitalStyle.SHADED: "sombreado",
    OrbitalStyle.SOLID: "sólido",
}

_ORBITAL_GLYPH_IDS = {
    OrbitalType.S: "s",
    OrbitalType.P: "p",
    OrbitalType.D: "d",
    OrbitalType.DZ2: "dz2",
    OrbitalType.F: "f",
    OrbitalType.FZ3: "fz3",
    OrbitalType.FXZ2: "fxz2",
    OrbitalType.FZX2Y2: "fzx2y2",
    OrbitalType.FXYZ: "fxyz",
    OrbitalType.SP3: "sp3",
    OrbitalType.SP2: "sp2",
    OrbitalType.SP_LOBE: "sp_lobe",
    OrbitalType.SIGMA_BONDING: "sigma_bonding",
    OrbitalType.PI_BONDING: "pi_bonding",
    OrbitalType.TORUS: "torus",
}

_CANVAS_EXTENT_SCALES = {
    OrbitalType.S: 0.84,
    OrbitalType.P: 0.92,
    OrbitalType.D: 0.92,
    OrbitalType.DZ2: 0.98,
    OrbitalType.F: 0.94,
    OrbitalType.FZ3: 0.98,
    OrbitalType.FXZ2: 0.98,
    OrbitalType.FZX2Y2: 0.98,
    OrbitalType.FXYZ: 0.98,
    OrbitalType.SP3: 0.94,
    OrbitalType.SP2: 0.92,
    OrbitalType.SP_LOBE: 0.88,
    OrbitalType.SIGMA_BONDING: 0.96,
    OrbitalType.PI_BONDING: 0.98,
    OrbitalType.TORUS: 0.84,
}


def _label_for(orbital_type: OrbitalType, style: OrbitalStyle) -> str:
    return f"{_TYPE_LABELS[orbital_type]} ({_STYLE_LABELS[style]})"


def _build_all_specs() -> dict[str, OrbitalGlyphSpec]:
    specs: dict[str, OrbitalGlyphSpec] = {}
    for orbital_type in OrbitalType:
        for style in OrbitalStyle:
            kind = f"{orbital_type.value}_{style.value}"
            specs[kind] = OrbitalGlyphSpec(
                kind=kind,
                orbital_type=orbital_type,
                glyph_id=_ORBITAL_GLYPH_IDS[orbital_type],
                style=style,
                label=_label_for(orbital_type, style),
                canvas_extent_scale=_CANVAS_EXTENT_SCALES[orbital_type],
                stroke_shaded_lobes=True,
            )
    return specs


ORBITAL_SPECS = _build_all_specs()

REFERENCE_PALETTE_CELLS = (
    ("s_outline", "s_shaded", "p_shaded", "p_solid", "d_shaded", "d_solid", None),
    ("torus_outline", "torus_shaded", "torus_solid", "sp3_shaded", "sp3_solid", "dz2_shaded", "dz2_solid"),
    ("sp_lobe_outline", "sp_lobe_shaded", "sp_lobe_solid", "sigma_bonding_shaded", "sigma_bonding_solid", "pi_bonding_shaded", "pi_bonding_solid"),
    ("sp2_outline", "sp2_shaded", "sp2_solid", None, None, None, None),
)

ORBITAL_PALETTE_MODEL = OrbitalPaletteModel(cells=REFERENCE_PALETTE_CELLS)
ORBITAL_MENU_COLUMNS = ORBITAL_PALETTE_MODEL.columns
ORBITAL_MENU_ROWS = ORBITAL_PALETTE_MODEL.rows
ORBITAL_MENU_ORDER = list(ORBITAL_PALETTE_MODEL.entries)


def orbital_tool_id(kind: str) -> str:
    """Devuelve el tool id asociado a una variante orbital."""
    return f"tool_orbital_{kind}"


def orbital_kind_from_tool_id(tool_id: str) -> str | None:
    """Extrae el kind desde un tool id Qt."""
    prefix = "tool_orbital_"
    if not str(tool_id).startswith(prefix):
        return None
    kind = str(tool_id)[len(prefix):]
    return kind if kind in ORBITAL_SPECS else None


def orbital_display_name(kind: str) -> str:
    """Nombre legible del orbital."""
    spec = ORBITAL_SPECS.get(kind) or ORBITAL_SPECS[DEFAULT_ORBITAL_KIND]
    return spec.label


def orbital_canvas_extent(kind: str, bond_length: float) -> float:
    """Longitud sugerida entre anchor0 y anchor1 para el orbital en canvas."""
    spec = ORBITAL_SPECS.get(kind) or ORBITAL_SPECS[DEFAULT_ORBITAL_KIND]
    return max(14.0, float(bond_length or 40.0) * spec.canvas_extent_scale)


def _default_family_presets() -> dict[str, OrbitalFamilyPreset]:
    return {
        "s": OrbitalFamilyPreset(
            family="s",
            builder="build_s_orbital",
            params=SphereParams(
                radius_x=34.0,
                radius_y=34.0,
                visual_padding=6.0,
                anchor_bias_x=0.0,
                anchor_bias_y=0.0,
                gradient_mode="radial",
                light_dir=_DEFAULT_LIGHT_DIR,
            ),
        ),
        "p": OrbitalFamilyPreset(
            family="p",
            builder="build_p_orbital",
            params=POrbitalParams(
                lobe_width=30.4,
                lobe_height=41.0,
                waist_width=0.355,
                shoulder_width=1.0,
                tip_roundness=0.71,
                node_gap=0.0,
                vertical_offset=0.0,
                outline_scale_x=1.0,
                outline_scale_y=1.0,
                visual_padding=4.5,
                anchor_bias_x=0.0,
                anchor_bias_y=0.0,
            ),
        ),
        "d": OrbitalFamilyPreset(
            family="d",
            builder="build_d_clover_orbital",
            params=CloverParams(
                vertical_width=30.4,
                vertical_height=41.0,
                horizontal_width=30.4,
                horizontal_height=41.0,
                tip_roundness=0.71,
                shoulder_width=1.0,
                waist_width=0.35,
                center_roundness=0.12,
                lobe_gap=0.0,
                rotation_deg=45.0,
                visual_padding=5.0,
                anchor_bias_x=0.0,
                anchor_bias_y=0.0,
            ),
        ),
        "dz2": OrbitalFamilyPreset(
            family="dz2",
            builder="build_dz2_orbital",
            params=Dz2Params(
                axial_width=28.0,
                axial_height=41.0,
                axial_tip_roundness=0.71,
                axial_shoulder_width=1.0,
                axial_waist_width=0.24,
                axial_center_roundness=0.0,
                torus_outer_width=62.0,
                torus_outer_height=9.6,
                torus_inner_width_ratio=0.46,
                torus_inner_height_ratio=0.26,
                torus_offset_y=0.0,
                visual_padding=6.0,
                anchor_bias_x=0.0,
                anchor_bias_y=-0.75,
            ),
        ),
        "torus": OrbitalFamilyPreset(
            family="torus",
            builder="build_torus_orbital",
            params=TorusParams(
                torus_outer_width=72.0,
                torus_outer_height=20.0,
                torus_inner_width_ratio=0.42,
                torus_inner_height_ratio=0.38,
                visual_padding=6.0,
                gradient_mode="elliptical",
            ),
        ),
        "sp2": OrbitalFamilyPreset(
            family="sp2",
            builder="build_sp2_orbital",
            params=Sp2Params(
                lobe=CanonicalLobeParams(
                    width=22.0,
                    height=38.0,
                    tip_roundness=0.76,
                    shoulder_width=0.95,
                    waist_width=0.16,
                    belly_width=0.88,
                    belly_height_ratio=0.58,
                ),
                center_roundness=0.0,
                pair_rotation_deg=120.0,
                planar_rotation_deg=0.0,
                top_scale=1.0,
                pair_scale=1.0,
                visual_padding=7.0,
                anchor_bias_x=0.0,
                anchor_bias_y=1.5,
                gradient_mode="linear",
            ),
        ),
        "sp3": OrbitalFamilyPreset(
            family="sp3",
            builder="build_sp3_orbital",
            params=Sp3Params(
                lobe=CanonicalLobeParams(
                    width=22.0,
                    height=38.0,
                    tip_roundness=0.76,
                    shoulder_width=0.95,
                    waist_width=0.16,
                ),
                center_roundness=0.0,
                side_rotation_deg=126.0,
                bottom_rotation_deg=180.0,
                side_scale=0.96,
                bottom_scale=0.88,
                visual_padding=7.0,
                anchor_bias_x=0.0,
                anchor_bias_y=2.0,
                gradient_mode="linear",
            ),
        ),
        "sp_lobe": OrbitalFamilyPreset(
            family="sp_lobe",
            builder="build_sp_lobe_orbital",
            params=HybridOrbitalParams(
                major_lobe=CanonicalLobeParams(
                    width=29.0,
                    height=46.0,
                    tip_roundness=0.56,
                    shoulder_width=0.94,
                    waist_width=0.34,
                    belly_width=0.94,
                    belly_height_ratio=0.56,
                ),
                minor_lobe=CanonicalLobeParams(
                    width=13.0,
                    height=18.0,
                    tip_roundness=0.78,
                    shoulder_width=0.92,
                    waist_width=0.42,
                    belly_width=0.84,
                    belly_height_ratio=0.58,
                ),
                visual_padding=5.0,
                anchor_bias_x=0.0,
                anchor_bias_y=0.0,
            ),
        ),
        "sigma_bonding": OrbitalFamilyPreset(
            family="sigma_bonding",
            builder="build_sigma_bonding_orbital",
            params=HybridOrbitalParams(
                major_lobe=CanonicalLobeParams(
                    width=30.0,
                    height=34.0,
                    tip_roundness=0.78,
                    shoulder_width=0.98,
                    waist_width=0.22,
                    belly_width=0.96,
                    belly_height_ratio=0.48,
                ),
                minor_lobe=CanonicalLobeParams(
                    width=12.0,
                    height=15.0,
                    tip_roundness=0.74,
                    shoulder_width=0.92,
                    waist_width=0.34,
                    belly_width=0.82,
                    belly_height_ratio=0.58,
                ),
                major_offset_x=11.0,
                minor_offset_x=27.0,
                major_offset_y=0.0,
                minor_offset_y=0.0,
                center_roundness=0.18,
                visual_padding=5.0,
                anchor_bias_x=0.0,
                anchor_bias_y=0.0,
                gradient_mode="linear",
                minor_gradient_mode="radial",
            ),
        ),
        "pi_bonding": OrbitalFamilyPreset(
            family="pi_bonding",
            builder="build_pi_bonding_orbital",
            params=PiBondingParams(
                upper_lobe=CanonicalLobeParams(
                    width=24.0,
                    height=24.0,
                    tip_roundness=0.72,
                    shoulder_width=0.96,
                    waist_width=0.26,
                    belly_width=0.90,
                    belly_height_ratio=0.52,
                ),
                lower_lobe=CanonicalLobeParams(
                    width=24.0,
                    height=24.0,
                    tip_roundness=0.72,
                    shoulder_width=0.96,
                    waist_width=0.26,
                    belly_width=0.90,
                    belly_height_ratio=0.52,
                ),
                node_gap=3.0,
                center_roundness=0.16,
                lateral_offset=5.0,
                bridge_width=20.0,
                bridge_height=11.0,
                cap_width=22.0,
                cap_height=6.0,
                visual_padding=6.0,
                anchor_bias_x=0.0,
                anchor_bias_y=0.0,
                gradient_mode="linear",
            ),
        ),
        "f": OrbitalFamilyPreset(
            family="f",
            builder="build_f_orbital",
            params=FOrbitalParams(
                lobes=(
                    TeardropParams(20.0, 26.0, 0.72, 0.96, 0.66, 0.90, 0.46, 0.26, 0.03, 0.0, rotation_deg=0.0, orientation="up", phase="positive"),
                    TeardropParams(20.0, 26.0, 0.72, 0.96, 0.66, 0.90, 0.46, 0.26, 0.03, 0.0, rotation_deg=60.0, orientation="up", phase="negative"),
                    TeardropParams(20.0, 26.0, 0.72, 0.96, 0.66, 0.90, 0.46, 0.26, 0.03, 0.0, rotation_deg=120.0, orientation="up", phase="positive"),
                    TeardropParams(20.0, 26.0, 0.72, 0.96, 0.66, 0.90, 0.46, 0.26, 0.03, 0.0, rotation_deg=180.0, orientation="up", phase="negative"),
                    TeardropParams(20.0, 26.0, 0.72, 0.96, 0.66, 0.90, 0.46, 0.26, 0.03, 0.0, rotation_deg=240.0, orientation="up", phase="positive"),
                    TeardropParams(20.0, 26.0, 0.72, 0.96, 0.66, 0.90, 0.46, 0.26, 0.03, 0.0, rotation_deg=300.0, orientation="up", phase="negative"),
                ),
                visual_padding=8.0,
            ),
        ),
        "fz3": OrbitalFamilyPreset(
            family="fz3",
            builder="build_f_orbital",
            params=FOrbitalParams(
                lobes=(
                    TeardropParams(18.0, 24.0, 0.72, 0.96, 0.66, 0.90, 0.44, 0.28, 0.03, 0.0, offset_y=-30.0, orientation="up", phase="positive"),
                    TeardropParams(18.0, 24.0, 0.72, 0.96, 0.66, 0.90, 0.44, 0.28, 0.03, 0.0, offset_y=30.0, orientation="down", phase="positive"),
                ),
                torus=TorusParams(
                    torus_outer_width=40.0,
                    torus_outer_height=10.0,
                    torus_inner_width_ratio=0.50,
                    torus_inner_height_ratio=0.42,
                    visual_padding=8.0,
                    phase="neutral",
                ),
                visual_padding=8.0,
            ),
        ),
        "fxz2": OrbitalFamilyPreset(
            family="fxz2",
            builder="build_f_orbital",
            params=FOrbitalParams(
                lobes=(
                    TeardropParams(16.0, 22.0, 0.72, 0.96, 0.66, 0.90, 0.44, 0.28, 0.03, 0.0, offset_y=-26.0, orientation="up", phase="positive"),
                    TeardropParams(16.0, 22.0, 0.72, 0.96, 0.66, 0.90, 0.44, 0.28, 0.03, 0.0, offset_y=26.0, orientation="down", phase="negative"),
                    TeardropParams(14.0, 26.0, 0.72, 0.96, 0.66, 0.88, 0.42, 0.24, 0.03, 0.0, rotation_deg=50.0, orientation="up", phase="positive"),
                    TeardropParams(14.0, 26.0, 0.72, 0.96, 0.66, 0.88, 0.42, 0.24, 0.03, 0.0, rotation_deg=130.0, orientation="up", phase="negative"),
                    TeardropParams(14.0, 26.0, 0.72, 0.96, 0.66, 0.88, 0.42, 0.24, 0.03, 0.0, rotation_deg=230.0, orientation="up", phase="positive"),
                    TeardropParams(14.0, 26.0, 0.72, 0.96, 0.66, 0.88, 0.42, 0.24, 0.03, 0.0, rotation_deg=310.0, orientation="up", phase="negative"),
                ),
                visual_padding=8.0,
            ),
        ),
        "fzx2y2": OrbitalFamilyPreset(
            family="fzx2y2",
            builder="build_f_orbital",
            params=FOrbitalParams(
                lobes=(
                    TeardropParams(18.0, 26.0, 0.72, 0.96, 0.66, 0.90, 0.46, 0.26, 0.03, 0.0, rotation_deg=0.0, orientation="up", phase="positive"),
                    TeardropParams(15.0, 21.0, 0.72, 0.94, 0.66, 0.86, 0.42, 0.24, 0.03, 0.0, rotation_deg=45.0, orientation="up", phase="negative"),
                    TeardropParams(18.0, 26.0, 0.72, 0.96, 0.66, 0.90, 0.46, 0.26, 0.03, 0.0, rotation_deg=90.0, orientation="up", phase="positive"),
                    TeardropParams(15.0, 21.0, 0.72, 0.94, 0.66, 0.86, 0.42, 0.24, 0.03, 0.0, rotation_deg=135.0, orientation="up", phase="negative"),
                    TeardropParams(18.0, 26.0, 0.72, 0.96, 0.66, 0.90, 0.46, 0.26, 0.03, 0.0, rotation_deg=180.0, orientation="up", phase="positive"),
                    TeardropParams(15.0, 21.0, 0.72, 0.94, 0.66, 0.86, 0.42, 0.24, 0.03, 0.0, rotation_deg=225.0, orientation="up", phase="negative"),
                    TeardropParams(18.0, 26.0, 0.72, 0.96, 0.66, 0.90, 0.46, 0.26, 0.03, 0.0, rotation_deg=270.0, orientation="up", phase="positive"),
                    TeardropParams(15.0, 21.0, 0.72, 0.94, 0.66, 0.86, 0.42, 0.24, 0.03, 0.0, rotation_deg=315.0, orientation="up", phase="negative"),
                ),
                visual_padding=8.0,
            ),
        ),
        "fxyz": OrbitalFamilyPreset(
            family="fxyz",
            builder="build_f_orbital",
            params=FOrbitalParams(
                lobes=(
                    TeardropParams(15.0, 20.0, 0.72, 0.94, 0.66, 0.88, 0.42, 0.24, 0.03, 0.0, rotation_deg=0.0, orientation="up", phase="positive"),
                    TeardropParams(15.0, 20.0, 0.72, 0.94, 0.66, 0.88, 0.42, 0.24, 0.03, 0.0, rotation_deg=45.0, orientation="up", phase="negative"),
                    TeardropParams(15.0, 20.0, 0.72, 0.94, 0.66, 0.88, 0.42, 0.24, 0.03, 0.0, rotation_deg=90.0, orientation="up", phase="positive"),
                    TeardropParams(15.0, 20.0, 0.72, 0.94, 0.66, 0.88, 0.42, 0.24, 0.03, 0.0, rotation_deg=135.0, orientation="up", phase="negative"),
                    TeardropParams(15.0, 20.0, 0.72, 0.94, 0.66, 0.88, 0.42, 0.24, 0.03, 0.0, rotation_deg=180.0, orientation="up", phase="positive"),
                    TeardropParams(15.0, 20.0, 0.72, 0.94, 0.66, 0.88, 0.42, 0.24, 0.03, 0.0, rotation_deg=225.0, orientation="up", phase="negative"),
                    TeardropParams(15.0, 20.0, 0.72, 0.94, 0.66, 0.88, 0.42, 0.24, 0.03, 0.0, rotation_deg=270.0, orientation="up", phase="positive"),
                    TeardropParams(15.0, 20.0, 0.72, 0.94, 0.66, 0.88, 0.42, 0.24, 0.03, 0.0, rotation_deg=315.0, orientation="up", phase="negative"),
                ),
                visual_padding=8.0,
            ),
        ),
    }


_DEFAULT_ORBITAL_FAMILY_PRESETS = _default_family_presets()
ORBITAL_FAMILY_ORDER = tuple(_DEFAULT_ORBITAL_FAMILY_PRESETS.keys())


def _dataclass_payload(obj: Any) -> Any:
    if isinstance(obj, tuple):
        return [_dataclass_payload(item) for item in obj]
    if isinstance(obj, list):
        return [_dataclass_payload(item) for item in obj]
    if isinstance(obj, dict):
        return {key: _dataclass_payload(value) for key, value in obj.items()}
    if hasattr(obj, "__dataclass_fields__"):
        return {key: _dataclass_payload(value) for key, value in asdict(obj).items()}
    return obj


def _docs_for_value(value: Any) -> Any:
    if hasattr(value, "__dataclass_fields__"):
        docs: dict[str, Any] = {}
        for key in value.__dataclass_fields__:
            child = getattr(value, key)
            nested = _docs_for_value(child)
            docs[key] = nested if nested else ORBITAL_PARAMETER_DESCRIPTIONS.get(key, key)
        return docs
    if isinstance(value, (tuple, list)) and value and hasattr(value[0], "__dataclass_fields__"):
        return {
            key: ORBITAL_PARAMETER_DESCRIPTIONS.get(key, key)
            for key in value[0].__dataclass_fields__
        }
    return {}


def _family_docs(family_key: str) -> dict[str, Any]:
    preset = _DEFAULT_ORBITAL_FAMILY_PRESETS[family_key]
    docs = _docs_for_value(preset.params)
    return docs if isinstance(docs, dict) else {}


def default_orbital_presets_payload(*, include_docs: bool = True) -> dict[str, Any]:
    """Payload serializable con los presets por defecto."""
    payload: dict[str, Any] = {
        "_meta": {
            "format_version": 1,
            "path_hint": str(ORBITAL_PRESET_CONFIG_PATH),
            "notes": "Edita estos parámetros y usa tools/orbital_tuner.py o tools/render_orbital_family_preview.py para previsualizar.",
        }
    }
    for family_key, preset in _DEFAULT_ORBITAL_FAMILY_PRESETS.items():
        family_payload = {
            "builder": preset.builder,
            "params": _dataclass_payload(preset.params),
        }
        if include_docs:
            family_payload["_docs"] = _family_docs(family_key)
        payload[family_key] = family_payload
    return payload


def _strip_meta(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: _strip_meta(item) for key, item in value.items() if not str(key).startswith("_")}
    if isinstance(value, list):
        return [_strip_meta(item) for item in value]
    return value


def _deep_merge(base: Any, override: Any) -> Any:
    if isinstance(base, dict) and isinstance(override, dict):
        merged = dict(base)
        for key, value in override.items():
            if key in merged:
                merged[key] = _deep_merge(merged[key], value)
            else:
                merged[key] = value
        return merged
    return override


def _load_json_payload(path: Path, *, include_docs: bool) -> dict[str, Any]:
    base = default_orbital_presets_payload(include_docs=include_docs)
    if not path.exists():
        return base
    raw = json.loads(path.read_text(encoding="utf-8"))
    return _deep_merge(base, _strip_meta(raw))


def load_orbital_presets_payload(
    path: Path = ORBITAL_PRESET_CONFIG_PATH,
    *,
    include_docs: bool = False,
) -> dict[str, Any]:
    """Carga el payload editable completo desde disco con fallback a defaults."""
    return _load_json_payload(path, include_docs=include_docs)


def _parse_teardrop(payload: dict[str, Any]) -> TeardropParams:
    data = _strip_meta(payload)
    return TeardropParams(**data)


def _parse_canonical_lobe(payload: dict[str, Any]) -> CanonicalLobeParams:
    data = _strip_meta(payload)
    allowed = set(CanonicalLobeParams.__dataclass_fields__)
    filtered = {key: value for key, value in data.items() if key in allowed}
    return CanonicalLobeParams(**filtered)


def _parse_torus(payload: dict[str, Any]) -> TorusParams:
    data = _strip_meta(payload)
    allowed = set(TorusParams.__dataclass_fields__)
    filtered = {key: value for key, value in data.items() if key in allowed}
    return TorusParams(**filtered)


def _canonical_from_legacy_teardrop(payload: dict[str, Any]) -> dict[str, Any]:
    data = _strip_meta(payload)
    return {
        "width": float(data.get("width", 24.0)),
        "height": float(data.get("height", 32.0)),
        "tip_roundness": float(data.get("tip_roundness", 0.68)),
        "shoulder_width": float(data.get("shoulder_width", 1.0)),
        "waist_width": float(data.get("waist_width", data.get("neck_width", 0.35))),
        "belly_width": float(data.get("belly_width", 0.0)),
        "belly_height_ratio": float(data.get("belly_height_ratio", 0.0)),
    }


def _normalize_p_params(payload: dict[str, Any]) -> dict[str, Any]:
    data = dict(_strip_meta(payload))
    if "waist_width" not in data:
        data["waist_width"] = data.get("neck_width", 0.35)
    if "node_gap" not in data:
        upper = float(data.get("upper_offset", 0.0))
        lower = float(data.get("lower_offset", 0.0))
        data["node_gap"] = abs(lower - upper)
    if "vertical_offset" not in data:
        upper = float(data.get("upper_offset", 0.0))
        lower = float(data.get("lower_offset", 0.0))
        data["vertical_offset"] = (upper + lower) * 0.5
    allowed = set(POrbitalParams.__dataclass_fields__)
    return {key: value for key, value in data.items() if key in allowed}


def _normalize_clover_params(payload: dict[str, Any]) -> dict[str, Any]:
    data = dict(_strip_meta(payload))
    data.setdefault("vertical_width", data.get("vertical_lobe_width", 29.6))
    data.setdefault("vertical_height", data.get("vertical_lobe_height", 38.0))
    data.setdefault("horizontal_width", data.get("horizontal_lobe_width", data.get("vertical_width", 29.6)))
    data.setdefault("horizontal_height", data.get("horizontal_lobe_height", data.get("vertical_height", 38.0)))
    data.setdefault("waist_width", data.get("neck_width", 0.35))
    data.setdefault("center_roundness", data.get("node_roundness", 0.12))
    data.setdefault("rotation_deg", data.get("clover_rotation_deg", 45.0))
    allowed = set(CloverParams.__dataclass_fields__)
    return {key: value for key, value in data.items() if key in allowed}


def _normalize_dz2_params(payload: dict[str, Any]) -> dict[str, Any]:
    data = dict(_strip_meta(payload))
    data.setdefault("axial_width", data.get("axial_lobe_width", 28.0))
    data.setdefault("axial_height", data.get("axial_lobe_height", 41.0))
    data.setdefault("axial_shoulder_width", data.get("axial_shoulder", 1.0))
    data.setdefault("axial_waist_width", data.get("axial_neck_width", 0.24))
    data.setdefault("axial_center_roundness", data.get("axial_node_roundness", 0.0))
    allowed = set(Dz2Params.__dataclass_fields__)
    return {key: value for key, value in data.items() if key in allowed}


def _parse_hybrid_payload(payload: dict[str, Any], *, family_key: str) -> HybridOrbitalParams:
    data = dict(_strip_meta(payload))
    if "major_lobe" in data and "minor_lobe" in data:
        return HybridOrbitalParams(
            major_lobe=_parse_canonical_lobe(data["major_lobe"]),
            minor_lobe=_parse_canonical_lobe(data["minor_lobe"]),
            major_offset_x=float(data.get("major_offset_x", 0.0)),
            minor_offset_x=float(data.get("minor_offset_x", 0.0)),
            major_offset_y=float(data.get("major_offset_y", 0.0)),
            minor_offset_y=float(data.get("minor_offset_y", 0.0)),
            center_roundness=float(data.get("center_roundness", 0.0)),
            node_gap=float(data.get("node_gap", 0.0)),
            visual_padding=float(data.get("visual_padding", 6.0)),
            anchor_bias_x=float(data.get("anchor_bias_x", 0.0)),
            anchor_bias_y=float(data.get("anchor_bias_y", 0.0)),
            gradient_mode=str(data.get("gradient_mode", "linear")),
            minor_gradient_mode=str(data.get("minor_gradient_mode", "linear")),
            light_dir=str(data.get("light_dir", _DEFAULT_LIGHT_DIR)),
        )

    primary = _strip_meta(data.get("primary_lobe", {}))
    secondary = _strip_meta(data.get("secondary_lobe", {}))
    if family_key == "sp3":
        major_payload = _canonical_from_legacy_teardrop(secondary)
        minor_payload = _canonical_from_legacy_teardrop(primary)
        major_offset_y = float(secondary.get("offset_y", 0.0))
        minor_offset_y = float(primary.get("offset_y", 0.0))
    else:
        major_payload = _canonical_from_legacy_teardrop(primary)
        minor_payload = _canonical_from_legacy_teardrop(secondary)
        major_offset_y = float(primary.get("offset_y", 0.0))
        minor_offset_y = float(secondary.get("offset_y", 0.0))

    return HybridOrbitalParams(
        major_lobe=_parse_canonical_lobe(major_payload),
        minor_lobe=_parse_canonical_lobe(minor_payload),
        major_offset_x=float(data.get("major_offset_x", 0.0)),
        minor_offset_x=float(data.get("minor_offset_x", 0.0)),
        major_offset_y=float(data.get("major_offset_y", major_offset_y)),
        minor_offset_y=float(data.get("minor_offset_y", minor_offset_y)),
        center_roundness=float(data.get("center_roundness", 0.0)),
        node_gap=float(data.get("node_gap", 0.0)),
        visual_padding=float(data.get("visual_padding", 6.0)),
        anchor_bias_x=float(data.get("anchor_bias_x", 0.0)),
        anchor_bias_y=float(data.get("anchor_bias_y", 0.0)),
        gradient_mode=str(data.get("gradient_mode", primary.get("gradient_mode", "linear"))),
        minor_gradient_mode=str(data.get("minor_gradient_mode", secondary.get("gradient_mode", "linear"))),
        light_dir=str(data.get("light_dir", _DEFAULT_LIGHT_DIR)),
    )


def _parse_sp3_payload(payload: dict[str, Any]) -> Sp3Params:
    data = dict(_strip_meta(payload))
    if "lobe" in data:
        return Sp3Params(
            lobe=_parse_canonical_lobe(data["lobe"]),
            center_roundness=float(data.get("center_roundness", 0.0)),
            side_rotation_deg=float(data.get("side_rotation_deg", 126.0)),
            bottom_rotation_deg=float(data.get("bottom_rotation_deg", 180.0)),
            side_scale=float(data.get("side_scale", 0.96)),
            bottom_scale=float(data.get("bottom_scale", 0.88)),
            visual_padding=float(data.get("visual_padding", 7.0)),
            anchor_bias_x=float(data.get("anchor_bias_x", 0.0)),
            anchor_bias_y=float(data.get("anchor_bias_y", 2.0)),
            gradient_mode=str(data.get("gradient_mode", "linear")),
            light_dir=str(data.get("light_dir", _DEFAULT_LIGHT_DIR)),
        )

    major_payload = _strip_meta(data.get("major_lobe", {}))
    legacy_lobe = major_payload if major_payload else {
        "width": 22.0,
        "height": 38.0,
        "tip_roundness": 0.76,
        "shoulder_width": 0.95,
        "waist_width": 0.16,
    }
    return Sp3Params(
        lobe=_parse_canonical_lobe(legacy_lobe),
        center_roundness=float(data.get("center_roundness", 0.0)),
        side_rotation_deg=float(data.get("side_rotation_deg", 126.0)),
        bottom_rotation_deg=float(data.get("bottom_rotation_deg", 180.0)),
        side_scale=float(data.get("side_scale", 0.96)),
        bottom_scale=float(data.get("bottom_scale", 0.88)),
        visual_padding=float(data.get("visual_padding", 7.0)),
        anchor_bias_x=float(data.get("anchor_bias_x", 0.0)),
        anchor_bias_y=float(data.get("anchor_bias_y", 2.0)),
        gradient_mode=str(data.get("gradient_mode", "linear")),
        light_dir=str(data.get("light_dir", _DEFAULT_LIGHT_DIR)),
    )


def _parse_sp2_payload(payload: dict[str, Any]) -> Sp2Params:
    data = dict(_strip_meta(payload))
    if "lobe" in data:
        return Sp2Params(
            lobe=_parse_canonical_lobe(data["lobe"]),
            center_roundness=float(data.get("center_roundness", 0.0)),
            pair_rotation_deg=float(data.get("pair_rotation_deg", 120.0)),
            planar_rotation_deg=float(data.get("planar_rotation_deg", 0.0)),
            top_scale=float(data.get("top_scale", 1.0)),
            pair_scale=float(data.get("pair_scale", 1.0)),
            visual_padding=float(data.get("visual_padding", 7.0)),
            anchor_bias_x=float(data.get("anchor_bias_x", 0.0)),
            anchor_bias_y=float(data.get("anchor_bias_y", 1.5)),
            gradient_mode=str(data.get("gradient_mode", "linear")),
            light_dir=str(data.get("light_dir", _DEFAULT_LIGHT_DIR)),
        )

    major_payload = _strip_meta(data.get("major_lobe", {}))
    legacy_lobe = major_payload if major_payload else {
        "width": 22.0,
        "height": 38.0,
        "tip_roundness": 0.76,
        "shoulder_width": 0.95,
        "waist_width": 0.16,
    }
    return Sp2Params(
        lobe=_parse_canonical_lobe(legacy_lobe),
        center_roundness=float(data.get("center_roundness", 0.0)),
        pair_rotation_deg=float(data.get("pair_rotation_deg", 120.0)),
        planar_rotation_deg=float(data.get("planar_rotation_deg", 0.0)),
        top_scale=float(data.get("top_scale", 1.0)),
        pair_scale=float(data.get("pair_scale", 1.0)),
        visual_padding=float(data.get("visual_padding", 7.0)),
        anchor_bias_x=float(data.get("anchor_bias_x", 0.0)),
        anchor_bias_y=float(data.get("anchor_bias_y", 1.5)),
        gradient_mode=str(data.get("gradient_mode", "linear")),
        light_dir=str(data.get("light_dir", _DEFAULT_LIGHT_DIR)),
    )


def _parse_pi_payload(payload: dict[str, Any]) -> PiBondingParams:
    data = dict(_strip_meta(payload))
    legacy_ring = _parse_torus(data["ring"]) if data.get("ring") else (_parse_torus(data["torus"]) if data.get("torus") else None)
    lateral_offset = float(data.get("lateral_offset", legacy_ring.torus_outer_width * 0.25 if legacy_ring is not None else 10.0))
    if "upper_lobe" in data and "lower_lobe" in data:
        return PiBondingParams(
            upper_lobe=_parse_canonical_lobe(data["upper_lobe"]),
            lower_lobe=_parse_canonical_lobe(data["lower_lobe"]),
            node_gap=float(data.get("node_gap", 0.0)),
            vertical_offset=float(data.get("vertical_offset", 0.0)),
            center_roundness=float(data.get("center_roundness", 0.16)),
            lateral_offset=lateral_offset,
            bridge_width=float(data.get("bridge_width", 20.0)),
            bridge_height=float(data.get("bridge_height", 11.0)),
            cap_width=float(data.get("cap_width", 22.0)),
            cap_height=float(data.get("cap_height", 6.0)),
            visual_padding=float(data.get("visual_padding", 6.0)),
            anchor_bias_x=float(data.get("anchor_bias_x", 0.0)),
            anchor_bias_y=float(data.get("anchor_bias_y", 0.0)),
            gradient_mode=str(data.get("gradient_mode", "linear")),
            light_dir=str(data.get("light_dir", _DEFAULT_LIGHT_DIR)),
        )

    primary = _strip_meta(data.get("primary_lobe", {}))
    secondary = _strip_meta(data.get("secondary_lobe", {}))
    upper = _canonical_from_legacy_teardrop(primary)
    lower = _canonical_from_legacy_teardrop(secondary)
    upper_offset = float(primary.get("offset_y", 0.0))
    lower_offset = float(secondary.get("offset_y", 0.0))
    return PiBondingParams(
        upper_lobe=_parse_canonical_lobe(upper),
        lower_lobe=_parse_canonical_lobe(lower),
        node_gap=float(data.get("node_gap", abs(lower_offset - upper_offset))),
        vertical_offset=float(data.get("vertical_offset", (upper_offset + lower_offset) * 0.5)),
        center_roundness=float(data.get("center_roundness", 0.16)),
        lateral_offset=lateral_offset,
        bridge_width=float(data.get("bridge_width", 20.0)),
        bridge_height=float(data.get("bridge_height", 11.0)),
        cap_width=float(data.get("cap_width", 22.0)),
        cap_height=float(data.get("cap_height", 6.0)),
        visual_padding=float(data.get("visual_padding", 6.0)),
        anchor_bias_x=float(data.get("anchor_bias_x", 0.0)),
        anchor_bias_y=float(data.get("anchor_bias_y", 0.0)),
        gradient_mode=str(data.get("gradient_mode", primary.get("gradient_mode", "linear"))),
        light_dir=str(data.get("light_dir", _DEFAULT_LIGHT_DIR)),
    )


def _parse_preset_payload(payload: dict[str, Any]) -> dict[str, OrbitalFamilyPreset]:
    families: dict[str, OrbitalFamilyPreset] = {}
    for family_key, default_preset in _DEFAULT_ORBITAL_FAMILY_PRESETS.items():
        family_payload = payload.get(family_key, {})
        params_payload = _strip_meta(family_payload.get("params", {}))
        builder = str(family_payload.get("builder") or default_preset.builder)
        if family_key == "s":
            params = SphereParams(**params_payload)
        elif family_key == "p":
            params = POrbitalParams(**_normalize_p_params(params_payload))
        elif family_key == "d":
            params = CloverParams(**_normalize_clover_params(params_payload))
        elif family_key == "dz2":
            params = Dz2Params(**_normalize_dz2_params(params_payload))
        elif family_key == "torus":
            params = _parse_torus(params_payload)
        elif family_key == "sp2":
            params = _parse_sp2_payload(params_payload)
        elif family_key == "sp3":
            params = _parse_sp3_payload(params_payload)
        elif family_key in {"sp_lobe", "sigma_bonding"}:
            params = _parse_hybrid_payload(params_payload, family_key=family_key)
        elif family_key == "pi_bonding":
            params = _parse_pi_payload(params_payload)
        elif family_key in {"f", "fz3", "fxz2", "fzx2y2", "fxyz"}:
            params = FOrbitalParams(
                lobes=tuple(_parse_teardrop(item) for item in params_payload.get("lobes", [])),
                torus=_parse_torus(params_payload["torus"]) if params_payload.get("torus") else None,
                visual_padding=float(params_payload.get("visual_padding", 8.0)),
                anchor_bias_x=float(params_payload.get("anchor_bias_x", 0.0)),
                anchor_bias_y=float(params_payload.get("anchor_bias_y", 0.0)),
            )
        else:
            raise KeyError(family_key)
        families[family_key] = OrbitalFamilyPreset(family=family_key, builder=builder, params=params)
    return families


def load_orbital_family_presets(path: Path = ORBITAL_PRESET_CONFIG_PATH) -> dict[str, OrbitalFamilyPreset]:
    """Carga presets externos con fallback a defaults internos."""
    return _parse_preset_payload(_load_json_payload(path, include_docs=False))


def save_orbital_presets_payload(payload: dict[str, Any], path: Path = ORBITAL_PRESET_CONFIG_PATH) -> Path:
    """Guarda un payload de presets externos en disco."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    return path


def save_default_orbital_presets(path: Path = ORBITAL_PRESET_CONFIG_PATH) -> Path:
    """Escribe en disco el archivo editable de presets por defecto."""
    return save_orbital_presets_payload(default_orbital_presets_payload(include_docs=True), path)


def save_current_orbital_presets(path: Path = ORBITAL_PRESET_CONFIG_PATH) -> Path:
    """Guarda el estado activo actual como archivo externo editable."""
    return save_orbital_presets_payload(current_orbital_presets_payload(include_docs=True), path)


def current_orbital_presets_payload(*, include_docs: bool = False) -> dict[str, Any]:
    """Devuelve el estado activo actual serializado a payload."""
    payload: dict[str, Any] = {"_meta": {"format_version": 1}}
    for family_key, preset in _ACTIVE_ORBITAL_FAMILY_PRESETS.items():
        family_payload = {
            "builder": preset.builder,
            "params": _dataclass_payload(preset.params),
        }
        if include_docs:
            family_payload["_docs"] = _family_docs(family_key)
        payload[family_key] = family_payload
    return payload


def _ellipse_path(x: float, y: float, w: float, h: float) -> QPainterPath:
    path = QPainterPath()
    path.addEllipse(QRectF(x, y, w, h))
    return path


def _ring_path(x: float, y: float, w: float, h: float, inner_x: float, inner_y: float) -> QPainterPath:
    path = QPainterPath()
    path.setFillRule(Qt.FillRule.OddEvenFill)
    path.addEllipse(QRectF(x, y, w, h))
    inset_x = w * (1.0 - inner_x) * 0.5
    inset_y = h * (1.0 - inner_y) * 0.5
    path.addEllipse(QRectF(x + inset_x, y + inset_y, w * inner_x, h * inner_y))
    return path


def _transform_path(
    path: QPainterPath,
    *,
    translate: tuple[float, float] = (0.0, 0.0),
    rotate: float = 0.0,
    scale_x: float = 1.0,
    scale_y: float = 1.0,
) -> QPainterPath:
    transform = QTransform()
    transform.translate(translate[0], translate[1])
    if rotate:
        transform.rotate(rotate)
    transform.scale(scale_x, scale_y)
    return transform.map(path)


def mirror_path(path: QPainterPath, *, axis: str) -> QPainterPath:
    """Espeja un path sobre el eje `x` o `y`."""
    if axis == "x":
        return _transform_path(path, scale_y=-1.0)
    if axis == "y":
        return _transform_path(path, scale_x=-1.0)
    raise ValueError(axis)


def compose_parts(*groups: GeometryPart | tuple[GeometryPart, ...]) -> tuple[GeometryPart, ...]:
    """Aplana grupos de partes manteniendo un orden explícito."""
    parts: list[GeometryPart] = []
    for group in groups:
        if isinstance(group, GeometryPart):
            parts.append(group)
        else:
            parts.extend(group)
    return tuple(parts)


def _bounds_for_parts(parts: tuple[GeometryPart, ...]) -> QRectF:
    bounds: QRectF | None = None
    for part in parts:
        rect = part.path.boundingRect()
        bounds = rect if bounds is None else bounds.united(rect)
    return bounds if bounds is not None else QRectF(-1.0, -1.0, 2.0, 2.0)


def scale_to_visual_box(parts: tuple[GeometryPart, ...], box: QRectF) -> tuple[GeometryPart, ...]:
    """Escala un conjunto de partes a una caja visual explícita."""
    source = _bounds_for_parts(parts)
    if source.width() <= 1e-6 or source.height() <= 1e-6:
        return parts
    scale_x = box.width() / source.width()
    scale_y = box.height() / source.height()
    transform = QTransform()
    transform.translate(box.center().x(), box.center().y())
    transform.scale(scale_x, scale_y)
    transform.translate(-source.center().x(), -source.center().y())
    return tuple(replace(part, path=transform.map(part.path)) for part in parts)


def _catmull_rom_path(points: list[QPointF]) -> QPainterPath:
    path = QPainterPath(points[0])
    if len(points) < 2:
        return path
    for index in range(len(points) - 1):
        p0 = points[index - 1] if index > 0 else points[index]
        p1 = points[index]
        p2 = points[index + 1]
        p3 = points[index + 2] if index + 2 < len(points) else points[index + 1]
        c1 = QPointF(
            p1.x() + (p2.x() - p0.x()) / 6.0,
            p1.y() + (p2.y() - p0.y()) / 6.0,
        )
        c2 = QPointF(
            p2.x() - (p3.x() - p1.x()) / 6.0,
            p2.y() - (p3.y() - p1.y()) / 6.0,
        )
        path.cubicTo(c1, c2, p2)
    return path


def _closed_catmull_rom_path(points: list[QPointF]) -> QPainterPath:
    path = QPainterPath(points[0])
    if len(points) < 3:
        return path
    count = len(points)
    for index in range(count):
        p0 = points[(index - 1) % count]
        p1 = points[index % count]
        p2 = points[(index + 1) % count]
        p3 = points[(index + 2) % count]
        c1 = QPointF(
            p1.x() + (p2.x() - p0.x()) / 6.0,
            p1.y() + (p2.y() - p0.y()) / 6.0,
        )
        c2 = QPointF(
            p2.x() - (p3.x() - p1.x()) / 6.0,
            p2.y() - (p3.y() - p1.y()) / 6.0,
        )
        path.cubicTo(c1, c2, p2)
    path.closeSubpath()
    return path


def build_teardrop(params: TeardropParams) -> GeometryPart:
    """Construye un lóbulo lágrima desde parámetros explícitos."""
    half_width = params.width * 0.5
    tip = QPointF(0.0, -params.height)
    shoulder_y = -params.height * params.shoulder_height_ratio
    belly_y = -params.height * params.belly_height_ratio
    neck_y = -params.height * params.neck_height
    node_y = params.height * params.cusp_depth
    tip_round = QPointF(half_width * params.tip_roundness, (tip.y() + shoulder_y) * 0.5)
    shoulder = QPointF(half_width * params.shoulder_width, shoulder_y)
    belly = QPointF(half_width * params.belly_width, belly_y)
    neck = QPointF(half_width * params.neck_width, neck_y)
    node = QPointF(0.0, node_y)

    right_points = [tip, tip_round, shoulder, belly, neck, node]
    right_curve = _catmull_rom_path(right_points)
    left_points = [node] + [QPointF(-point.x(), point.y()) for point in reversed(right_points[1:-1])] + [tip]
    left_curve = _catmull_rom_path(left_points)
    path = QPainterPath()
    path.addPath(right_curve)
    path.connectPath(left_curve)
    path.closeSubpath()

    if params.outline_scale_x != 1.0 or params.outline_scale_y != 1.0:
        path = _transform_path(path, scale_x=params.outline_scale_x, scale_y=params.outline_scale_y)
    if params.mirror:
        path = mirror_path(path, axis="y")
    if params.orientation == "down":
        path = mirror_path(path, axis="x")
    elif params.orientation == "left":
        path = _transform_path(path, rotate=-90.0)
    elif params.orientation == "right":
        path = _transform_path(path, rotate=90.0)
    elif params.orientation != "up":
        raise ValueError(f"Orientación no soportada: {params.orientation}")
    if params.rotation_deg:
        path = _transform_path(path, rotate=params.rotation_deg)
    if params.offset_x or params.offset_y:
        path = _transform_path(path, translate=(params.offset_x, params.offset_y))

    return GeometryPart(
        name="lobe",
        path=path,
        phase=params.phase,
        gradient_mode=params.gradient_mode,
        light_dir=params.light_dir,
    )


def build_teardrop_as_canonical_lobe(
    params: TeardropParams,
    *,
    template: "_LobeTemplate",
    smooth_profile: bool = True,
    root_roundness: float = 0.0,
) -> GeometryPart:
    """Reinterpreta un teardrop legacy usando la silueta canónica moderna."""
    canonical = CanonicalLobeParams(
        width=params.width,
        height=params.height,
        tip_roundness=params.tip_roundness,
        shoulder_width=params.shoulder_width,
        waist_width=params.neck_width,
        belly_width=params.belly_width,
        belly_height_ratio=params.belly_height_ratio,
    )
    path = _canonical_lobe_path(
        canonical,
        template,
        orientation=params.orientation,
        offset_x=params.offset_x,
        offset_y=params.offset_y,
        scale_x=params.outline_scale_x,
        scale_y=params.outline_scale_y,
        root_roundness=root_roundness,
        rotation_deg=params.rotation_deg,
        smooth_profile=smooth_profile,
    )
    if params.mirror:
        path = mirror_path(path, axis="y")
    return GeometryPart(
        name="lobe",
        path=path,
        phase=params.phase,
        gradient_mode=params.gradient_mode,
        light_dir=params.light_dir,
    )


@dataclass(frozen=True)
class _LobeTemplate:
    shoulder_y: float
    belly_y: float
    belly_floor: float
    belly_gap: float
    neck_y: float


_P_TEMPLATE = _LobeTemplate(shoulder_y=0.71, belly_y=0.48, belly_floor=0.86, belly_gap=0.18, neck_y=0.06)
_D_VERTICAL_TEMPLATE = _P_TEMPLATE
_D_HORIZONTAL_TEMPLATE = _P_TEMPLATE
_DZ2_TEMPLATE = _LobeTemplate(shoulder_y=0.67, belly_y=0.44, belly_floor=0.84, belly_gap=0.17, neck_y=0.16)
_HYBRID_MAJOR_TEMPLATE = _LobeTemplate(shoulder_y=0.70, belly_y=0.46, belly_floor=0.84, belly_gap=0.16, neck_y=0.08)
_HYBRID_MINOR_TEMPLATE = _LobeTemplate(shoulder_y=0.76, belly_y=0.58, belly_floor=0.82, belly_gap=0.16, neck_y=0.14)
_SIGMA_MAJOR_TEMPLATE = _LobeTemplate(shoulder_y=0.55, belly_y=0.26, belly_floor=0.86, belly_gap=0.18, neck_y=0.06)
_SIGMA_MINOR_TEMPLATE = _LobeTemplate(shoulder_y=0.76, belly_y=0.54, belly_floor=0.86, belly_gap=0.14, neck_y=0.12)
_PI_TEMPLATE = _LobeTemplate(shoulder_y=0.67, belly_y=0.44, belly_floor=0.84, belly_gap=0.16, neck_y=0.10)


def _clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, float(value)))


def _canonical_lobe_path(
    params: CanonicalLobeParams,
    template: _LobeTemplate,
    *,
    orientation: str = "up",
    offset_x: float = 0.0,
    offset_y: float = 0.0,
    scale_x: float = 1.0,
    scale_y: float = 1.0,
    root_roundness: float = 0.0,
    rotation_deg: float = 0.0,
    smooth_profile: bool = False,
) -> QPainterPath:
    width = _clamp(params.width, 6.0, 120.0)
    height = _clamp(params.height, 8.0, 140.0)
    tip_roundness = _clamp(params.tip_roundness, 0.18, 0.90)
    waist_width = _clamp(params.waist_width, 0.12, 0.74)
    shoulder_width = _clamp(params.shoulder_width, waist_width + 0.12, 1.0)
    custom_belly = params.belly_width > 1e-6 or params.belly_height_ratio > 1e-6
    belly_y = _clamp(
        params.belly_height_ratio if params.belly_height_ratio > 1e-6 else template.belly_y,
        template.neck_y + 0.14,
        template.shoulder_y - 0.04,
    )
    auto_belly_width = max(template.belly_floor, shoulder_width * 0.88, waist_width + template.belly_gap)
    belly_width = _clamp(
        params.belly_width if params.belly_width > 1e-6 else auto_belly_width,
        max(shoulder_width * 0.80, waist_width + template.belly_gap, 0.55),
        0.98,
    )

    half_width = width * 0.5
    tip = QPointF(0.0, -height)
    belly = QPointF(half_width * belly_width, -height * belly_y)
    root_roundness = _clamp(root_roundness, 0.0, 0.38)
    root_half_width = half_width * root_roundness
    root_depth = height * root_roundness * 0.14
    node = QPointF(0.0, 0.0)
    node_right = QPointF(root_half_width, 0.0)
    node_left = QPointF(-root_half_width, 0.0)

    tip_shoulder_width = shoulder_width
    if custom_belly:
        tip_shoulder_width = _clamp(max(shoulder_width, belly_width * 0.64), shoulder_width, belly_width)
    tip_ctrl_1 = QPointF(half_width * tip_roundness, -height * 0.97)
    tip_ctrl_2 = QPointF(half_width * tip_shoulder_width, -height * template.shoulder_y)
    root_ctrl_1 = QPointF(
        half_width * min(1.0, belly_width + 0.05),
        -height * max(template.neck_y + 0.08, belly_y * (0.54 if custom_belly else 0.48)),
    )
    root_ctrl_2 = QPointF(half_width * waist_width, -height * template.neck_y)

    if smooth_profile:
        tip_ctrl_2 = QPointF(belly.x(), tip_ctrl_2.y())
        root_ctrl_1 = QPointF(belly.x(), 2.0 * belly.y() - tip_ctrl_2.y())
    path = QPainterPath(tip)
    path.cubicTo(tip_ctrl_1, tip_ctrl_2, belly)
    if root_roundness > 1e-6:
        path.cubicTo(
            root_ctrl_1,
            QPointF(max(root_ctrl_2.x(), root_half_width * 1.05), root_ctrl_2.y()),
            node_right,
        )
        path.cubicTo(
            QPointF(root_half_width * 0.35, root_depth),
            QPointF(-root_half_width * 0.35, root_depth),
            node_left,
        )
        path.cubicTo(
            QPointF(min(-root_ctrl_2.x(), -root_half_width * 1.05), root_ctrl_2.y()),
            QPointF(-root_ctrl_1.x(), root_ctrl_1.y()),
            QPointF(-belly.x(), belly.y()),
        )
    else:
        path.cubicTo(root_ctrl_1, root_ctrl_2, node)
        path.cubicTo(
            QPointF(-root_ctrl_2.x(), root_ctrl_2.y()),
            QPointF(-root_ctrl_1.x(), root_ctrl_1.y()),
            QPointF(-belly.x(), belly.y()),
        )
    path.cubicTo(
        QPointF(-tip_ctrl_2.x(), tip_ctrl_2.y()),
        QPointF(-tip_ctrl_1.x(), tip_ctrl_1.y()),
        tip,
    )
    path.closeSubpath()

    if scale_x != 1.0 or scale_y != 1.0:
        path = _transform_path(path, scale_x=scale_x, scale_y=scale_y)
    if orientation == "down":
        path = mirror_path(path, axis="x")
    elif orientation == "left":
        path = _transform_path(path, rotate=-90.0)
    elif orientation == "right":
        path = _transform_path(path, rotate=90.0)
    elif orientation != "up":
        raise ValueError(f"Orientación no soportada: {orientation}")
    if rotation_deg:
        path = _transform_path(path, rotate=rotation_deg)
    if offset_x or offset_y:
        path = _transform_path(path, translate=(offset_x, offset_y))
    return path


def _canonical_part(
    name: str,
    params: CanonicalLobeParams,
    template: _LobeTemplate,
    *,
    orientation: str,
    offset_x: float = 0.0,
    offset_y: float = 0.0,
    scale_x: float = 1.0,
    scale_y: float = 1.0,
    root_roundness: float = 0.0,
    rotation_deg: float = 0.0,
    smooth_profile: bool = False,
    phase: str,
    gradient_mode: str,
    light_dir: str,
) -> GeometryPart:
    return GeometryPart(
        name=name,
        path=_canonical_lobe_path(
            params,
            template,
            orientation=orientation,
            offset_x=offset_x,
            offset_y=offset_y,
            scale_x=scale_x,
            scale_y=scale_y,
            root_roundness=root_roundness,
            rotation_deg=rotation_deg,
            smooth_profile=smooth_profile,
        ),
        phase=phase,
        gradient_mode=gradient_mode,
        light_dir=light_dir,
    )


def build_s_orbital(params: SphereParams) -> FamilyBuildResult:
    """Builder canónico para el orbital s."""
    path = _ellipse_path(-params.radius_x, -params.radius_y, params.radius_x * 2.0, params.radius_y * 2.0)
    part = GeometryPart("sphere", path, phase="positive", gradient_mode=params.gradient_mode, light_dir=params.light_dir)
    return FamilyBuildResult(
        family="s",
        parts=(part,),
        outline_parts=("sphere",),
        shaded_fill_parts=("sphere",),
        shaded_outline_parts=(),
        solid_fill_parts=("sphere",),
        solid_outline_parts=(),
        visual_padding=params.visual_padding,
        anchor_bias_x=params.anchor_bias_x,
        anchor_bias_y=params.anchor_bias_y,
    )


def build_torus_orbital(params: TorusParams) -> FamilyBuildResult:
    """Builder canónico para un toroide limpio."""
    outer_width = _clamp(params.torus_outer_width, 18.0, 140.0)
    outer_height = _clamp(params.torus_outer_height, 6.0, 64.0)
    inner_width_ratio = _clamp(params.torus_inner_width_ratio, 0.18, 0.82)
    inner_height_ratio = _clamp(params.torus_inner_height_ratio, 0.12, 0.78)
    ring = _ring_path(
        -outer_width * 0.5,
        params.torus_offset_y - outer_height * 0.5,
        outer_width,
        outer_height,
        inner_width_ratio,
        inner_height_ratio,
    )
    part = GeometryPart("ring", ring, phase=params.phase, gradient_mode=params.gradient_mode, light_dir=params.light_dir)
    return FamilyBuildResult(
        family="torus",
        parts=(part,),
        outline_parts=("ring",),
        shaded_fill_parts=("ring",),
        shaded_outline_parts=(),
        solid_fill_parts=("ring",),
        solid_outline_parts=(),
        visual_padding=params.visual_padding,
        anchor_bias_x=params.anchor_bias_x,
        anchor_bias_y=params.anchor_bias_y,
    )


def build_p_orbital(params: POrbitalParams) -> FamilyBuildResult:
    """Builder canónico para un orbital p reconocible."""
    base = CanonicalLobeParams(
        width=params.lobe_width,
        height=params.lobe_height,
        tip_roundness=params.tip_roundness,
        shoulder_width=params.shoulder_width,
        waist_width=params.waist_width,
    )
    gap = _clamp(params.node_gap, 0.0, 12.0) * 0.5
    top = _canonical_part(
        "top",
        base,
        _P_TEMPLATE,
        orientation="up",
        offset_y=params.vertical_offset - gap,
        scale_x=_clamp(params.outline_scale_x, 0.75, 1.35),
        scale_y=_clamp(params.outline_scale_y, 0.75, 1.35),
        smooth_profile=True,
        phase="positive",
        gradient_mode=params.gradient_mode,
        light_dir=params.light_dir,
    )
    bottom = _canonical_part(
        "bottom",
        base,
        _P_TEMPLATE,
        orientation="down",
        offset_y=params.vertical_offset + gap,
        scale_x=_clamp(params.outline_scale_x, 0.75, 1.35),
        scale_y=_clamp(params.outline_scale_y, 0.75, 1.35),
        smooth_profile=True,
        phase="negative",
        gradient_mode=params.gradient_mode,
        light_dir=params.light_dir,
    )
    return FamilyBuildResult(
        family="p",
        parts=(top, bottom),
        outline_parts=("top", "bottom"),
        shaded_fill_parts=("top",),
        shaded_outline_parts=("bottom",),
        solid_fill_parts=("top",),
        solid_outline_parts=("bottom",),
        visual_padding=params.visual_padding,
        anchor_bias_x=params.anchor_bias_x,
        anchor_bias_y=params.anchor_bias_y,
    )


def build_d_clover_orbital(params: CloverParams) -> FamilyBuildResult:
    """Builder canónico para un trébol d de cuatro lóbulos independientes."""
    vertical = CanonicalLobeParams(
        width=params.vertical_width,
        height=params.vertical_height,
        tip_roundness=params.tip_roundness,
        shoulder_width=params.shoulder_width,
        waist_width=params.waist_width,
    )
    horizontal = CanonicalLobeParams(
        width=params.horizontal_width,
        height=params.horizontal_height,
        tip_roundness=params.tip_roundness,
        shoulder_width=params.shoulder_width,
        waist_width=params.waist_width,
    )
    gap = _clamp(params.lobe_gap, 0.0, 14.0) * 0.5
    parts = (
        _canonical_part("top", vertical, _D_VERTICAL_TEMPLATE, orientation="up", offset_y=-gap, root_roundness=params.center_roundness, smooth_profile=True, phase="positive", gradient_mode=params.gradient_mode, light_dir=params.light_dir),
        _canonical_part("bottom", vertical, _D_VERTICAL_TEMPLATE, orientation="down", offset_y=gap, root_roundness=params.center_roundness, smooth_profile=True, phase="positive", gradient_mode=params.gradient_mode, light_dir=params.light_dir),
        _canonical_part("left", horizontal, _D_HORIZONTAL_TEMPLATE, orientation="left", offset_x=-gap, root_roundness=params.center_roundness, smooth_profile=True, phase="negative", gradient_mode=params.gradient_mode, light_dir=params.light_dir),
        _canonical_part("right", horizontal, _D_HORIZONTAL_TEMPLATE, orientation="right", offset_x=gap, root_roundness=params.center_roundness, smooth_profile=True, phase="negative", gradient_mode=params.gradient_mode, light_dir=params.light_dir),
    )
    rotation_deg = float(params.rotation_deg)
    if rotation_deg:
        parts = tuple(replace(part, path=_transform_path(part.path, rotate=rotation_deg)) for part in parts)
    return FamilyBuildResult(
        family="d",
        parts=parts,
        outline_parts=("top", "bottom", "left", "right"),
        shaded_fill_parts=("top", "bottom"),
        shaded_outline_parts=("left", "right"),
        solid_fill_parts=("top", "bottom"),
        solid_outline_parts=("left", "right"),
        visual_padding=params.visual_padding,
        anchor_bias_x=params.anchor_bias_x,
        anchor_bias_y=params.anchor_bias_y,
    )


def build_dz2_orbital(params: Dz2Params) -> FamilyBuildResult:
    """Builder canónico para un dz2 con lóbulos axiales y toroide ecuatorial."""
    axial = CanonicalLobeParams(
        width=params.axial_width,
        height=params.axial_height,
        tip_roundness=params.axial_tip_roundness,
        shoulder_width=params.axial_shoulder_width,
        waist_width=params.axial_waist_width,
    )
    top = _canonical_part(
        "top",
        axial,
        _P_TEMPLATE,
        orientation="up",
        root_roundness=params.axial_center_roundness,
        phase="positive",
        gradient_mode=params.gradient_mode,
        light_dir=params.light_dir,
    )
    bottom = _canonical_part(
        "bottom",
        axial,
        _P_TEMPLATE,
        orientation="down",
        root_roundness=params.axial_center_roundness,
        phase="positive",
        gradient_mode=params.gradient_mode,
        light_dir=params.light_dir,
    )
    ring = build_torus_orbital(
        TorusParams(
            torus_outer_width=params.torus_outer_width,
            torus_outer_height=params.torus_outer_height,
            torus_inner_width_ratio=params.torus_inner_width_ratio,
            torus_inner_height_ratio=params.torus_inner_height_ratio,
            torus_offset_y=params.torus_offset_y,
            visual_padding=params.visual_padding,
            anchor_bias_x=params.anchor_bias_x,
            anchor_bias_y=params.anchor_bias_y,
            gradient_mode=params.torus_gradient_mode,
            light_dir=params.light_dir,
            phase="neutral",
        )
    ).parts[0]
    ring = replace(ring, name="ring")
    return FamilyBuildResult(
        family="dz2",
        parts=(top, ring, bottom),
        outline_parts=("top", "ring", "bottom"),
        shaded_fill_parts=("top", "bottom", "ring"),
        shaded_outline_parts=(),
        solid_fill_parts=("top", "bottom", "ring"),
        solid_outline_parts=(),
        visual_padding=params.visual_padding,
        anchor_bias_x=params.anchor_bias_x,
        anchor_bias_y=params.anchor_bias_y,
    )


def _build_hybrid_result(
    family: str,
    params: HybridOrbitalParams,
    *,
    major_orientation: str,
    minor_orientation: str,
    major_template: _LobeTemplate,
    minor_template: _LobeTemplate,
    outline_parts: tuple[str, ...],
    shaded_fill_parts: tuple[str, ...],
    shaded_outline_parts: tuple[str, ...],
    solid_fill_parts: tuple[str, ...],
    solid_outline_parts: tuple[str, ...],
) -> FamilyBuildResult:
    gap = _clamp(params.node_gap, 0.0, 12.0) * 0.5
    major = _canonical_part(
        "major",
        params.major_lobe,
        major_template,
        orientation=major_orientation,
        offset_y=params.major_offset_y + gap,
        phase="positive",
        gradient_mode=params.gradient_mode,
        light_dir=params.light_dir,
    )
    minor = _canonical_part(
        "minor",
        params.minor_lobe,
        minor_template,
        orientation=minor_orientation,
        offset_y=params.minor_offset_y - gap,
        phase="negative",
        gradient_mode=params.minor_gradient_mode,
        light_dir=params.light_dir,
    )
    return FamilyBuildResult(
        family=family,
        parts=(minor, major),
        outline_parts=outline_parts,
        shaded_fill_parts=shaded_fill_parts,
        shaded_outline_parts=shaded_outline_parts,
        solid_fill_parts=solid_fill_parts,
        solid_outline_parts=solid_outline_parts,
        visual_padding=params.visual_padding,
        anchor_bias_x=params.anchor_bias_x,
        anchor_bias_y=params.anchor_bias_y,
    )


def build_sp_lobe_orbital(params: HybridOrbitalParams) -> FamilyBuildResult:
    """Builder canónico para un lóbulo sp."""
    gap = _clamp(params.node_gap, 0.0, 12.0) * 0.5
    major = _canonical_part(
        "major",
        params.major_lobe,
        _HYBRID_MAJOR_TEMPLATE,
        orientation="right",
        offset_x=gap,
        offset_y=params.major_offset_y,
        smooth_profile=True,
        phase="positive",
        gradient_mode=params.gradient_mode,
        light_dir=params.light_dir,
    )
    minor = _canonical_part(
        "minor",
        params.minor_lobe,
        _HYBRID_MINOR_TEMPLATE,
        orientation="left",
        offset_x=-gap,
        offset_y=params.minor_offset_y,
        smooth_profile=True,
        phase="negative",
        gradient_mode=params.minor_gradient_mode,
        light_dir=params.light_dir,
    )
    return FamilyBuildResult(
        family="sp_lobe",
        parts=(minor, major),
        outline_parts=("minor", "major"),
        shaded_fill_parts=("minor", "major"),
        shaded_outline_parts=(),
        solid_fill_parts=("minor", "major"),
        solid_outline_parts=(),
        visual_padding=params.visual_padding,
        anchor_bias_x=params.anchor_bias_x,
        anchor_bias_y=params.anchor_bias_y,
    )


def build_sp3_orbital(params: Sp3Params) -> FamilyBuildResult:
    """Builder canónico para un orbital sp3."""
    parts = (
        _canonical_part(
            "bottom",
            params.lobe,
            _P_TEMPLATE,
            orientation="up",
            scale_x=params.bottom_scale,
            scale_y=params.bottom_scale,
            root_roundness=params.center_roundness,
            rotation_deg=params.bottom_rotation_deg,
            smooth_profile=True,
            phase="positive",
            gradient_mode=params.gradient_mode,
            light_dir=params.light_dir,
        ),
        _canonical_part(
            "left",
            params.lobe,
            _P_TEMPLATE,
            orientation="up",
            scale_x=params.side_scale,
            scale_y=params.side_scale,
            root_roundness=params.center_roundness,
            rotation_deg=-params.side_rotation_deg,
            smooth_profile=True,
            phase="positive",
            gradient_mode=params.gradient_mode,
            light_dir=params.light_dir,
        ),
        _canonical_part(
            "right",
            params.lobe,
            _P_TEMPLATE,
            orientation="up",
            scale_x=params.side_scale,
            scale_y=params.side_scale,
            root_roundness=params.center_roundness,
            rotation_deg=params.side_rotation_deg,
            smooth_profile=True,
            phase="positive",
            gradient_mode=params.gradient_mode,
            light_dir=params.light_dir,
        ),
        _canonical_part(
            "top",
            params.lobe,
            _P_TEMPLATE,
            orientation="up",
            root_roundness=params.center_roundness,
            smooth_profile=True,
            phase="positive",
            gradient_mode=params.gradient_mode,
            light_dir=params.light_dir,
        ),
    )
    part_names = tuple(part.name for part in parts)
    return FamilyBuildResult(
        family="sp3",
        parts=parts,
        outline_parts=part_names,
        shaded_fill_parts=part_names,
        shaded_outline_parts=(),
        solid_fill_parts=part_names,
        solid_outline_parts=(),
        visual_padding=params.visual_padding,
        anchor_bias_x=params.anchor_bias_x,
        anchor_bias_y=params.anchor_bias_y,
    )


def build_sp2_orbital(params: Sp2Params) -> FamilyBuildResult:
    """Builder canónico para un orbital sp2 trigonal plano."""
    parts = (
        _canonical_part(
            "top",
            params.lobe,
            _P_TEMPLATE,
            orientation="up",
            scale_x=params.top_scale,
            scale_y=params.top_scale,
            root_roundness=params.center_roundness,
            rotation_deg=params.planar_rotation_deg,
            smooth_profile=True,
            phase="positive",
            gradient_mode=params.gradient_mode,
            light_dir=params.light_dir,
        ),
        _canonical_part(
            "left",
            params.lobe,
            _P_TEMPLATE,
            orientation="up",
            scale_x=params.pair_scale,
            scale_y=params.pair_scale,
            root_roundness=params.center_roundness,
            rotation_deg=params.planar_rotation_deg - params.pair_rotation_deg,
            smooth_profile=True,
            phase="positive",
            gradient_mode=params.gradient_mode,
            light_dir=params.light_dir,
        ),
        _canonical_part(
            "right",
            params.lobe,
            _P_TEMPLATE,
            orientation="up",
            scale_x=params.pair_scale,
            scale_y=params.pair_scale,
            root_roundness=params.center_roundness,
            rotation_deg=params.planar_rotation_deg + params.pair_rotation_deg,
            smooth_profile=True,
            phase="positive",
            gradient_mode=params.gradient_mode,
            light_dir=params.light_dir,
        ),
    )
    part_names = tuple(part.name for part in parts)
    return FamilyBuildResult(
        family="sp2",
        parts=parts,
        outline_parts=part_names,
        shaded_fill_parts=part_names,
        shaded_outline_parts=(),
        solid_fill_parts=part_names,
        solid_outline_parts=(),
        visual_padding=params.visual_padding,
        anchor_bias_x=params.anchor_bias_x,
        anchor_bias_y=params.anchor_bias_y,
    )


def build_sigma_bonding_orbital(params: HybridOrbitalParams) -> FamilyBuildResult:
    """Builder canónico para sigma enlazante."""
    major_offset_x = _clamp(params.major_offset_x, 0.0, 40.0)
    minor_offset_x = _clamp(max(params.minor_offset_x, major_offset_x + params.major_lobe.height * 0.46), 0.0, 90.0)
    bond = GeometryPart(
        "bond",
        _ellipse_path(
            -max(params.major_lobe.height + major_offset_x * 1.8, params.major_lobe.height * 1.3) * 0.5,
            params.major_offset_y - params.major_lobe.width * 0.52,
            max(params.major_lobe.height + major_offset_x * 1.8, params.major_lobe.height * 1.3),
            params.major_lobe.width * 1.04,
        ),
        phase="positive",
        gradient_mode=params.gradient_mode,
        light_dir=params.light_dir,
    )
    minor_left = _canonical_part(
        "minor_left",
        params.minor_lobe,
        _SIGMA_MINOR_TEMPLATE,
        orientation="left",
        offset_x=-minor_offset_x,
        offset_y=params.minor_offset_y,
        smooth_profile=True,
        phase="negative",
        gradient_mode=params.minor_gradient_mode,
        light_dir=params.light_dir,
    )
    minor_right = _canonical_part(
        "minor_right",
        params.minor_lobe,
        _SIGMA_MINOR_TEMPLATE,
        orientation="right",
        offset_x=minor_offset_x,
        offset_y=params.minor_offset_y,
        smooth_profile=True,
        phase="negative",
        gradient_mode=params.minor_gradient_mode,
        light_dir=params.light_dir,
    )
    return FamilyBuildResult(
        family="sigma_bonding",
        parts=(minor_left, bond, minor_right),
        outline_parts=("minor_left", "bond", "minor_right"),
        shaded_fill_parts=("minor_left", "bond", "minor_right"),
        shaded_outline_parts=(),
        solid_fill_parts=("minor_left", "bond", "minor_right"),
        solid_outline_parts=(),
        visual_padding=params.visual_padding,
        anchor_bias_x=params.anchor_bias_x,
        anchor_bias_y=params.anchor_bias_y,
    )


def _build_pi_cloud_path(
    lobe: CanonicalLobeParams,
    *,
    orientation: str,
    base_y: float,
    lateral: float,
    center_roundness: float,
    bridge_width: float,
    bridge_height: float,
    cap_width: float,
    cap_height: float,
) -> QPainterPath:
    """Construye una nube pi continua tipo riñón para la mitad superior o inferior."""
    shoulder_width = _clamp(lobe.shoulder_width, lobe.waist_width + 0.12, 1.0)
    belly_width = _clamp(
        lobe.belly_width if lobe.belly_width > 1e-6 else max(shoulder_width * 0.92, 0.84),
        0.55,
        0.98,
    )
    belly_height = _clamp(
        lobe.belly_height_ratio if lobe.belly_height_ratio > 1e-6 else _PI_TEMPLATE.belly_y,
        _PI_TEMPLATE.neck_y + 0.14,
        _PI_TEMPLATE.shoulder_y - 0.04,
    )
    waist_width = _clamp(lobe.waist_width, 0.12, 0.74)
    tip_roundness = _clamp(lobe.tip_roundness, 0.18, 0.90)
    center_roundness = _clamp(center_roundness, 0.0, 0.38)

    outer_width = _clamp(
        max(
            lobe.width * (0.84 + 0.54 * max(shoulder_width, belly_width)),
            lobe.width + lateral * 2.8,
            (cap_width * 1.2) if cap_width > 1e-6 else 0.0,
        ),
        10.0,
        150.0,
    )
    outer_height = _clamp(
        max(
            lobe.height * (0.78 + 0.18 * tip_roundness),
            lobe.height * (0.68 + 0.24 * belly_height),
        ),
        8.0,
        150.0,
    )
    default_notch_width = outer_width * (0.16 + 0.30 * waist_width)
    default_notch_height = outer_height * (0.14 + 0.18 * belly_height + 0.10 * (1.0 - center_roundness))
    notch_width = _clamp(
        bridge_width if bridge_width > 1e-6 else default_notch_width,
        outer_width * 0.12,
        outer_width * 0.72,
    )
    notch_height = _clamp(
        bridge_height if bridge_height > 1e-6 else default_notch_height,
        outer_height * 0.10,
        outer_height * 0.62,
    )
    crown_width = _clamp(
        cap_width if cap_width > 1e-6 else outer_width * (0.72 + 0.12 * tip_roundness),
        outer_width * 0.45,
        outer_width * 1.05,
    )
    crown_height = _clamp(
        cap_height if cap_height > 1e-6 else outer_height * (0.10 + 0.08 * tip_roundness),
        outer_height * 0.05,
        outer_height * 0.36,
    )
    half_outer = outer_width * 0.5
    half_crown = crown_width * 0.5
    half_notch = notch_width * 0.5
    apex_y = -outer_height - crown_height * (0.08 + 0.08 * tip_roundness)
    shoulder_y = -outer_height * (0.82 + 0.04 * (1.0 - tip_roundness))
    belly_y = -outer_height * (0.46 + 0.08 * belly_height)
    lower_y = -outer_height * (0.10 + 0.04 * center_roundness)
    inner_y = -notch_height * (0.16 + 0.08 * center_roundness)
    notch_y = -notch_height * (0.82 + 0.08 * (1.0 - center_roundness))

    top_arc_x = half_crown * (0.46 + 0.06 * tip_roundness)
    shoulder_x = half_outer * (0.70 + 0.08 * shoulder_width)
    belly_x = half_outer * (0.92 + 0.06 * belly_width)
    lower_x = half_outer * (0.82 + 0.04 * belly_width)
    inner_x = max(half_notch * (0.88 + 0.08 * center_roundness), half_outer * (0.18 + 0.08 * waist_width))
    notch_side_x = half_notch * (0.56 + 0.10 * center_roundness)
    notch_side_y = (inner_y + notch_y) * 0.5

    points = [
        QPointF(top_arc_x, apex_y),
        QPointF(shoulder_x, shoulder_y),
        QPointF(belly_x, belly_y),
        QPointF(lower_x, lower_y),
        QPointF(inner_x, inner_y),
        QPointF(notch_side_x, notch_side_y),
        QPointF(0.0, notch_y),
        QPointF(-notch_side_x, notch_side_y),
        QPointF(-inner_x, inner_y),
        QPointF(-lower_x, lower_y),
        QPointF(-belly_x, belly_y),
        QPointF(-shoulder_x, shoulder_y),
        QPointF(-top_arc_x, apex_y),
    ]
    path = _closed_catmull_rom_path(points).simplified()

    if orientation == "down":
        path = mirror_path(path, axis="x")
    elif orientation != "up":
        raise ValueError(f"Orientación no soportada: {orientation}")
    if base_y:
        path = _transform_path(path, translate=(0.0, base_y))
    return path


def build_pi_bonding_orbital(params: PiBondingParams) -> FamilyBuildResult:
    """Builder canónico para pi enlazante."""
    gap = _clamp(params.node_gap, 0.0, 12.0) * 0.5
    lateral = _clamp(params.lateral_offset, 0.0, 40.0)
    upper = GeometryPart(
        "upper",
        _build_pi_cloud_path(
            params.upper_lobe,
            orientation="up",
            base_y=params.vertical_offset - gap,
            lateral=lateral,
            center_roundness=params.center_roundness,
            bridge_width=params.bridge_width,
            bridge_height=params.bridge_height,
            cap_width=params.cap_width,
            cap_height=params.cap_height,
        ),
        phase="positive",
        gradient_mode=params.gradient_mode,
        light_dir=params.light_dir,
    )
    lower = GeometryPart(
        "lower",
        _build_pi_cloud_path(
            params.lower_lobe,
            orientation="down",
            base_y=params.vertical_offset + gap,
            lateral=lateral,
            center_roundness=params.center_roundness,
            bridge_width=params.bridge_width,
            bridge_height=params.bridge_height,
            cap_width=params.cap_width,
            cap_height=params.cap_height,
        ),
        phase="negative",
        gradient_mode=params.gradient_mode,
        light_dir=params.light_dir,
    )
    parts = (upper, lower)
    return FamilyBuildResult(
        family="pi_bonding",
        parts=parts,
        outline_parts=("upper", "lower"),
        shaded_fill_parts=("upper", "lower"),
        shaded_outline_parts=(),
        solid_fill_parts=("upper", "lower"),
        solid_outline_parts=(),
        visual_padding=params.visual_padding,
        anchor_bias_x=params.anchor_bias_x,
        anchor_bias_y=params.anchor_bias_y,
    )


def build_f_orbital(params: FOrbitalParams, *, family: str) -> FamilyBuildResult:
    """Builder explícito para familias f/fz3 a partir de una lista de lóbulos."""
    if family in {"f", "fz3", "fxz2", "fzx2y2", "fxyz"}:
        root_roundness = 0.05 if family == "f" else (0.03 if family in {"fxz2", "fzx2y2", "fxyz"} else 0.0)
        parts = [
            replace(
                build_teardrop_as_canonical_lobe(
                    lobe,
                    template=_P_TEMPLATE,
                    smooth_profile=True,
                    root_roundness=root_roundness,
                ),
                name=f"lobe_{index}",
            )
            for index, lobe in enumerate(params.lobes)
        ]
    else:
        parts = [replace(build_teardrop(lobe), name=f"lobe_{index}") for index, lobe in enumerate(params.lobes)]
    if params.torus is not None:
        if family == "fz3":
            ring_delta = params.torus.torus_outer_height * 0.75
            upper_ring = replace(
                build_torus(replace(params.torus, torus_offset_y=params.torus.torus_offset_y - ring_delta)).parts[0],
                name="ring_upper",
            )
            lower_ring = replace(
                build_torus(replace(params.torus, torus_offset_y=params.torus.torus_offset_y + ring_delta)).parts[0],
                name="ring_lower",
            )
            parts.extend((upper_ring, lower_ring))
        else:
            parts.append(replace(build_torus(params.torus).parts[0], name="ring"))
    part_names = tuple(part.name for part in parts)
    return FamilyBuildResult(
        family=family,
        parts=tuple(parts),
        outline_parts=part_names,
        shaded_fill_parts=part_names,
        shaded_outline_parts=(),
        solid_fill_parts=part_names,
        solid_outline_parts=(),
        visual_padding=params.visual_padding,
        anchor_bias_x=params.anchor_bias_x,
        anchor_bias_y=params.anchor_bias_y,
    )


def build_sphere(params: SphereParams) -> FamilyBuildResult:
    return build_s_orbital(params)


def build_torus(params: TorusParams) -> FamilyBuildResult:
    return build_torus_orbital(params)


def build_d_clover(params: CloverParams) -> FamilyBuildResult:
    return build_d_clover_orbital(params)


def build_dz2(params: Dz2Params) -> FamilyBuildResult:
    return build_dz2_orbital(params)


def build_sp_lobe(params: HybridOrbitalParams) -> FamilyBuildResult:
    return build_sp_lobe_orbital(params)


def build_sp2(params: Sp2Params) -> FamilyBuildResult:
    return build_sp2_orbital(params)


def build_sp3(params: Sp3Params) -> FamilyBuildResult:
    return build_sp3_orbital(params)


def build_sigma_bonding(params: HybridOrbitalParams) -> FamilyBuildResult:
    return build_sigma_bonding_orbital(params)


def build_pi_bonding(params: PiBondingParams) -> FamilyBuildResult:
    return build_pi_bonding_orbital(params)


_BUILDER_MAP = {
    "build_s_orbital": lambda preset: build_s_orbital(preset.params),  # type: ignore[arg-type]
    "build_sphere": lambda preset: build_s_orbital(preset.params),  # type: ignore[arg-type]
    "build_p_orbital": lambda preset: build_p_orbital(preset.params),  # type: ignore[arg-type]
    "build_d_clover_orbital": lambda preset: build_d_clover_orbital(preset.params),  # type: ignore[arg-type]
    "build_d_clover": lambda preset: build_d_clover_orbital(preset.params),  # type: ignore[arg-type]
    "build_dz2_orbital": lambda preset: build_dz2_orbital(preset.params),  # type: ignore[arg-type]
    "build_dz2": lambda preset: build_dz2_orbital(preset.params),  # type: ignore[arg-type]
    "build_torus_orbital": lambda preset: build_torus_orbital(preset.params),  # type: ignore[arg-type]
    "build_torus": lambda preset: build_torus_orbital(preset.params),  # type: ignore[arg-type]
    "build_sp2_orbital": lambda preset: build_sp2_orbital(preset.params),  # type: ignore[arg-type]
    "build_sp2": lambda preset: build_sp2_orbital(preset.params),  # type: ignore[arg-type]
    "build_sp3_orbital": lambda preset: build_sp3_orbital(preset.params),  # type: ignore[arg-type]
    "build_sp3": lambda preset: build_sp3_orbital(preset.params),  # type: ignore[arg-type]
    "build_sp_lobe_orbital": lambda preset: build_sp_lobe_orbital(preset.params),  # type: ignore[arg-type]
    "build_sp_lobe": lambda preset: build_sp_lobe_orbital(preset.params),  # type: ignore[arg-type]
    "build_sigma_bonding_orbital": lambda preset: build_sigma_bonding_orbital(preset.params),  # type: ignore[arg-type]
    "build_sigma_bonding": lambda preset: build_sigma_bonding_orbital(preset.params),  # type: ignore[arg-type]
    "build_pi_bonding_orbital": lambda preset: build_pi_bonding_orbital(preset.params),  # type: ignore[arg-type]
    "build_pi_bonding": lambda preset: build_pi_bonding_orbital(preset.params),  # type: ignore[arg-type]
    "build_f_orbital": lambda preset: build_f_orbital(preset.params, family=preset.family),  # type: ignore[arg-type]
}


def _parts_by_name(parts: tuple[GeometryPart, ...]) -> dict[str, GeometryPart]:
    return {part.name: part for part in parts}


def apply_outline_style(parts: tuple[GeometryPart, ...], names: tuple[str, ...] | None = None) -> tuple[GlyphLayer, ...]:
    """Aplica el estilo outline sin alterar la silueta base."""
    part_map = _parts_by_name(parts)
    target_names = names or tuple(part_map.keys())
    return tuple(GlyphLayer(part_map[name].path, "outline", name=name) for name in target_names)


def apply_shaded_style(
    parts: tuple[GeometryPart, ...],
    fill_names: tuple[str, ...],
    *,
    outline_names: tuple[str, ...] = (),
) -> tuple[GlyphLayer, ...]:
    """Aplica sombreado usando exactamente la misma geometría base."""
    part_map = _parts_by_name(parts)
    layers = [
        GlyphLayer(
            part_map[name].path,
            "shaded",
            name=name,
            gradient=part_map[name].gradient_mode,
            phase=part_map[name].phase,
            light_dir=part_map[name].light_dir,
        )
        for name in fill_names
    ]
    layers.extend(GlyphLayer(part_map[name].path, "outline") for name in outline_names)
    return tuple(layers)


def apply_solid_style(
    parts: tuple[GeometryPart, ...],
    fill_names: tuple[str, ...],
    *,
    outline_names: tuple[str, ...] = (),
) -> tuple[GlyphLayer, ...]:
    """Aplica relleno sólido usando exactamente la misma geometría base."""
    part_map = _parts_by_name(parts)
    layers = [
        GlyphLayer(
            part_map[name].path,
            "solid",
            name=name,
            phase=part_map[name].phase,
            light_dir=part_map[name].light_dir,
        )
        for name in fill_names
    ]
    layers.extend(GlyphLayer(part_map[name].path, "outline") for name in outline_names)
    return tuple(layers)


def _bounds_for_layers(*groups: tuple[GlyphLayer, ...]) -> QRectF:
    bounds: QRectF | None = None
    for group in groups:
        for layer in group:
            rect = layer.path.boundingRect()
            bounds = rect if bounds is None else bounds.united(rect)
    return bounds if bounds is not None else QRectF(-1.0, -1.0, 2.0, 2.0)


def _glyph_from_build_result(result: FamilyBuildResult) -> GlyphDefinition:
    outline = apply_outline_style(result.parts, result.outline_parts)
    shaded = apply_shaded_style(result.parts, result.shaded_fill_parts, outline_names=result.shaded_outline_parts)
    solid = apply_solid_style(result.parts, result.solid_fill_parts, outline_names=result.solid_outline_parts)
    return GlyphDefinition(
        id=result.family,
        paths_outline=outline,
        paths_shaded=shaded,
        paths_solid=solid,
        default_bounds=_bounds_for_layers(outline, shaded, solid),
        anchor_center=QPointF(0.0, 0.0),
        anchor_bias=QPointF(result.anchor_bias_x, result.anchor_bias_y),
        visual_padding=result.visual_padding,
    )


def _build_glyph_library_from_presets(presets: dict[str, OrbitalFamilyPreset]) -> dict[str, GlyphDefinition]:
    library: dict[str, GlyphDefinition] = {}
    for family_key, preset in presets.items():
        builder = _BUILDER_MAP[preset.builder]
        library[family_key] = _glyph_from_build_result(builder(preset))
    return library


def build_orbital_renderer(payload: dict[str, Any] | None = None) -> "OrbitalRenderer":
    """Crea un renderer aislado a partir de un payload o de los presets activos."""
    presets = _parse_preset_payload(_deep_merge(default_orbital_presets_payload(include_docs=False), _strip_meta(payload or {})))
    library = _build_glyph_library_from_presets(presets)
    return OrbitalRenderer(glyph_library=library, family_presets=presets)


_ACTIVE_ORBITAL_FAMILY_PRESETS = load_orbital_family_presets()
ORBITAL_FAMILY_PRESETS = _ACTIVE_ORBITAL_FAMILY_PRESETS
ORBITAL_GEOMETRY_PRESETS = ORBITAL_FAMILY_PRESETS
ORBITAL_GLYPH_LIBRARY = _build_glyph_library_from_presets(_ACTIVE_ORBITAL_FAMILY_PRESETS)


class OrbitalRenderer:
    """Renderer compartido por toolbar, documento y herramientas de tuning."""

    def __init__(
        self,
        *,
        glyph_library: dict[str, GlyphDefinition] | None = None,
        family_presets: dict[str, OrbitalFamilyPreset] | None = None,
    ) -> None:
        self._glyph_library = glyph_library or dict(ORBITAL_GLYPH_LIBRARY)
        self._family_presets = family_presets or dict(_ACTIVE_ORBITAL_FAMILY_PRESETS)

    def set_glyph_library(
        self,
        glyph_library: dict[str, GlyphDefinition],
        family_presets: dict[str, OrbitalFamilyPreset],
    ) -> None:
        self._glyph_library = glyph_library
        self._family_presets = family_presets

    def spec_for_kind(self, kind: str) -> OrbitalGlyphSpec:
        return ORBITAL_SPECS.get(kind) or ORBITAL_SPECS[DEFAULT_ORBITAL_KIND]

    def glyph_for_spec(self, spec: OrbitalGlyphSpec) -> GlyphDefinition:
        return self._glyph_library[spec.glyph_id]

    def family_preset(self, family: str) -> OrbitalFamilyPreset:
        return self._family_presets[family]

    @staticmethod
    def _layers_for_style(glyph: GlyphDefinition, style: OrbitalStyle) -> tuple[GlyphLayer, ...]:
        if style == OrbitalStyle.OUTLINE:
            return glyph.paths_outline
        if style == OrbitalStyle.SHADED:
            return glyph.paths_shaded
        return glyph.paths_solid

    @staticmethod
    def _reference_extent(glyph: GlyphDefinition) -> float:
        center = glyph.anchor_center
        bounds = glyph.default_bounds
        return max(
            abs(center.x() - bounds.left()),
            abs(center.x() - bounds.right()),
            abs(center.y() - bounds.top()),
            abs(center.y() - bounds.bottom()),
            1.0,
        )

    def _anchor_transform(self, glyph: GlyphDefinition, anchor0: QPointF, anchor1: QPointF) -> QTransform:
        dx = anchor1.x() - anchor0.x()
        dy = anchor1.y() - anchor0.y()
        extent = math.hypot(dx, dy)
        if extent <= 1e-6:
            extent = 1.0
            dx = 0.0
            dy = -1.0
        axis_x = dx / extent
        axis_y = dy / extent
        perp_x = -axis_y
        perp_y = axis_x
        scale = extent / self._reference_extent(glyph)
        translated_center = glyph.anchor_center
        center_to_origin = QTransform()
        center_to_origin.translate(-translated_center.x(), -translated_center.y())
        basis = QTransform(
            perp_x * scale,
            perp_y * scale,
            -axis_x * scale,
            -axis_y * scale,
            anchor0.x(),
            anchor0.y(),
        )
        return center_to_origin * basis

    def transformed_layers(
        self,
        spec: OrbitalGlyphSpec,
        anchor0: QPointF,
        anchor1: QPointF,
        *,
        part_styles: dict[str, dict[str, Any]] | None = None,
    ) -> tuple[GlyphLayer, ...]:
        glyph = self.glyph_for_spec(spec)
        transform = self._anchor_transform(glyph, anchor0, anchor1)
        return tuple(
            GlyphLayer(
                transform.map(self._layer_path_with_offset(layer, part_styles)),
                layer.paint,
                name=layer.name,
                gradient=layer.gradient,
                stroke=layer.stroke,
                phase=layer.phase,
                light_dir=layer.light_dir,
            )
            for layer in self._layers_for_style(glyph, spec.style)
        )

    def combined_path(
        self,
        spec: OrbitalGlyphSpec,
        anchor0: QPointF,
        anchor1: QPointF,
        *,
        part_styles: dict[str, dict[str, Any]] | None = None,
    ) -> QPainterPath:
        path = QPainterPath()
        for layer in self.transformed_layers(spec, anchor0, anchor1, part_styles=part_styles):
            path.addPath(layer.path)
        return path

    def bounding_rect(
        self,
        spec: OrbitalGlyphSpec,
        anchor0: QPointF,
        anchor1: QPointF,
        *,
        part_styles: dict[str, dict[str, Any]] | None = None,
    ) -> QRectF:
        return self.combined_path(spec, anchor0, anchor1, part_styles=part_styles).boundingRect()

    def _stroke_width(self, glyph: GlyphDefinition, anchor0: QPointF, anchor1: QPointF) -> float:
        extent = math.hypot(anchor1.x() - anchor0.x(), anchor1.y() - anchor0.y())
        scale = extent / self._reference_extent(glyph)
        return max(1.1, min(2.2, scale * 1.9))

    @staticmethod
    def _gradient_stops(layer: GlyphLayer) -> tuple[tuple[float, str], ...]:
        if layer.phase == "negative":
            return _GRADIENT_STOP_PRESETS["negative"]
        return _GRADIENT_STOP_PRESETS.get(layer.gradient, _GRADIENT_STOP_PRESETS["linear"])

    @staticmethod
    def _light_vector(light_dir: str | None) -> tuple[float, float]:
        return _LIGHT_DIRECTION_VECTORS.get(light_dir or _DEFAULT_LIGHT_DIR, _LIGHT_DIRECTION_VECTORS[_DEFAULT_LIGHT_DIR])

    @classmethod
    def _gradient_for_path(
        cls,
        layer: GlyphLayer,
        *,
        color_override: QColor | None = None,
        opacity: float = 1.0,
    ) -> QBrush:
        rect = layer.path.boundingRect()
        center = rect.center()
        vx, vy = cls._light_vector(layer.light_dir)
        opacity = max(0.0, min(1.0, float(opacity)))
        if layer.gradient in {"radial", "elliptical"}:
            radius = max(rect.width(), rect.height(), 1.0) * 0.92
            focus = QPointF(center.x() - vx * radius * 0.18, center.y() - vy * radius * 0.18)
            gradient = QRadialGradient(center, radius, focus)
            for stop, color in cls._gradient_stops(layer):
                gradient.setColorAt(stop, cls._gradient_stop_color(color, stop, layer.phase, color_override, opacity))
            brush = QBrush(gradient)
            if layer.gradient == "elliptical":
                transform = QTransform()
                transform.translate(center.x(), center.y())
                if rect.width() >= rect.height():
                    transform.scale(1.0, max(rect.width() / max(rect.height(), 1.0), 1.0))
                else:
                    transform.scale(max(rect.height() / max(rect.width(), 1.0), 1.0), 1.0)
                transform.translate(-center.x(), -center.y())
                brush.setTransform(transform)
            return brush
        start = QPointF(center.x() - vx * rect.width() * 0.58, center.y() - vy * rect.height() * 0.58)
        end = QPointF(center.x() + vx * rect.width() * 0.58, center.y() + vy * rect.height() * 0.58)
        gradient = QLinearGradient(start, end)
        for stop, color in cls._gradient_stops(layer):
            gradient.setColorAt(stop, cls._gradient_stop_color(color, stop, layer.phase, color_override, opacity))
        return QBrush(gradient)

    @staticmethod
    def _outline_pen(width: float, *, color: QColor | None = None, opacity: float = 1.0) -> QPen:
        stroke_color = QColor(color or _STROKE_COLOR)
        stroke_color.setAlphaF(max(0.0, min(1.0, float(opacity))))
        return QPen(
            stroke_color,
            width,
            Qt.PenStyle.SolidLine,
            Qt.PenCapStyle.RoundCap,
            Qt.PenJoinStyle.RoundJoin,
        )

    @staticmethod
    def _gradient_stop_color(
        fallback: str,
        stop: float,
        phase: str | None,
        color_override: QColor | None,
        opacity: float,
    ) -> QColor:
        if color_override is None or not color_override.isValid():
            color = QColor(fallback)
            color.setAlphaF(opacity)
            return color
        factor = (0.82 + stop * 0.42) if phase == "negative" else (1.38 - stop * 0.62)
        factor = max(0.60, min(1.60, factor))
        color = QColor(color_override)
        if factor >= 1.0:
            color = color.lighter(int(round(factor * 100.0)))
        else:
            color = color.darker(int(round((1.0 / max(factor, 0.05)) * 100.0)))
        color.setAlphaF(opacity)
        return color

    @staticmethod
    def _part_style(part_styles: dict[str, dict[str, Any]] | None, name: str) -> dict[str, Any]:
        if not part_styles or not name:
            return {}
        raw = part_styles.get(name)
        return raw if isinstance(raw, dict) else {}

    @classmethod
    def _layer_path_with_offset(
        cls,
        layer: GlyphLayer,
        part_styles: dict[str, dict[str, Any]] | None,
    ) -> QPainterPath:
        part_style = cls._part_style(part_styles, layer.name)
        try:
            offset_x = float(part_style.get("offset_x", 0.0) or 0.0)
        except Exception:
            offset_x = 0.0
        try:
            offset_y = float(part_style.get("offset_y", 0.0) or 0.0)
        except Exception:
            offset_y = 0.0
        if abs(offset_x) <= 1e-6 and abs(offset_y) <= 1e-6:
            return layer.path
        return _transform_path(layer.path, translate=(offset_x, offset_y))

    def part_names(self, spec: OrbitalGlyphSpec) -> tuple[str, ...]:
        glyph = self.glyph_for_spec(spec)
        ordered: list[str] = []
        seen: set[str] = set()
        for group in (glyph.paths_outline, glyph.paths_shaded, glyph.paths_solid):
            for layer in group:
                if not layer.name or layer.name in seen:
                    continue
                seen.add(layer.name)
                ordered.append(layer.name)
        return tuple(ordered)

    def _paint_layer(
        self,
        painter: QPainter,
        layer: GlyphLayer,
        stroke_width: float,
        *,
        part_styles: dict[str, dict[str, Any]] | None = None,
    ) -> None:
        part_style = self._part_style(part_styles, layer.name)
        color_override = QColor(str(part_style["color"])) if part_style.get("color") else None
        if color_override is not None and not color_override.isValid():
            color_override = None
        try:
            opacity = max(0.0, min(1.0, float(part_style.get("opacity", 1.0))))
        except Exception:
            opacity = 1.0
        if layer.paint == "outline":
            painter.setPen(self._outline_pen(stroke_width, color=color_override, opacity=opacity))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawPath(layer.path)
            return
        if layer.paint == "shaded":
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(self._gradient_for_path(layer, color_override=color_override, opacity=opacity))
            painter.drawPath(layer.path)
            return
        painter.setPen(Qt.PenStyle.NoPen)
        solid_color = QColor(color_override or _SOLID_FILL)
        solid_color.setAlphaF(opacity)
        painter.setBrush(QBrush(solid_color))
        painter.drawPath(layer.path)

    def paint_glyph(
        self,
        painter: QPainter,
        spec: OrbitalGlyphSpec,
        anchor0: QPointF,
        anchor1: QPointF,
        *,
        stroke_shaded_lobes: bool | None = None,
        part_styles: dict[str, dict[str, Any]] | None = None,
    ) -> None:
        del stroke_shaded_lobes
        glyph = self.glyph_for_spec(spec)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
        stroke_width = self._stroke_width(glyph, anchor0, anchor1)
        for layer in self.transformed_layers(spec, anchor0, anchor1, part_styles=part_styles):
            self._paint_layer(painter, layer, stroke_width, part_styles=part_styles)

    def canonical_anchors(self, spec: OrbitalGlyphSpec, rect: QRectF) -> tuple[QPointF, QPointF]:
        glyph = self.glyph_for_spec(spec)
        inner = rect.adjusted(4.0, 4.0, -4.0, -4.0)
        padded = glyph.default_bounds.adjusted(
            -glyph.visual_padding,
            -glyph.visual_padding,
            glyph.visual_padding,
            glyph.visual_padding,
        )
        scale = min(inner.width() / padded.width(), inner.height() / padded.height())
        scale = max(0.1, float(scale))
        anchor0 = QPointF(
            inner.center().x() - glyph.anchor_center.x() * scale,
            inner.center().y() - glyph.anchor_center.y() * scale,
        )
        anchor1 = QPointF(anchor0.x(), anchor0.y() - self._reference_extent(glyph) * scale)
        bbox = self.bounding_rect(spec, anchor0, anchor1)
        target_center = inner.center() + QPointF(glyph.anchor_bias.x() * scale, glyph.anchor_bias.y() * scale)
        delta = target_center - bbox.center()
        return anchor0 + delta, anchor1 + delta

    def draw_icon(self, kind: str, size: int = ORBITAL_ICON_SIZE) -> QIcon:
        spec = self.spec_for_kind(kind)
        pixmap = QPixmap(max(1, int(size)), max(1, int(size)))
        pixmap.fill(Qt.GlobalColor.transparent)
        painter = QPainter(pixmap)
        anchor0, anchor1 = self.canonical_anchors(spec, QRectF(0.0, 0.0, float(size), float(size)))
        self.paint_glyph(painter, spec, anchor0, anchor1)
        painter.end()
        return QIcon(pixmap)

    def render_palette_image(
        self,
        model: OrbitalPaletteModel = ORBITAL_PALETTE_MODEL,
        *,
        background: QColor = _PALETTE_BG,
    ) -> QImage:
        image = QImage(model.image_size, QImage.Format.Format_ARGB32_Premultiplied)
        image.fill(background)
        painter = QPainter(image)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        for row in range(model.rows):
            for column in range(model.columns):
                cell = QRectF(
                    float(column * model.cell_width),
                    float(row * model.cell_height),
                    float(model.cell_width),
                    float(model.cell_height),
                )
                painter.fillRect(cell, _CELL_BG)
                painter.setPen(QPen(_CELL_LINE, 1.0))
                painter.setBrush(Qt.BrushStyle.NoBrush)
                painter.drawRect(cell.adjusted(0.0, 0.0, -1.0, -1.0))
                kind = model.kind_at(row, column)
                if not kind:
                    continue
                spec = self.spec_for_kind(kind)
                icon_rect = cell.adjusted(
                    float(model.inner_padding),
                    float(model.inner_padding),
                    float(-model.inner_padding),
                    float(-model.inner_padding),
                )
                anchor0, anchor1 = self.canonical_anchors(spec, icon_rect)
                self.paint_glyph(painter, spec, anchor0, anchor1)
        painter.end()
        return image

    def render_style_triptych(self, family: str, panel_size: int = 180, gap: int = 18) -> QImage:
        image = QImage(panel_size * 3 + gap * 4, panel_size + gap * 2, QImage.Format.Format_ARGB32_Premultiplied)
        image.fill(QColor("#FFFFFF"))
        painter = QPainter(image)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        styles = (OrbitalStyle.OUTLINE, OrbitalStyle.SHADED, OrbitalStyle.SOLID)
        for index, style in enumerate(styles):
            kind = f"{family}_{style.value}"
            spec = self.spec_for_kind(kind)
            rect = QRectF(float(gap + index * (panel_size + gap)), float(gap), float(panel_size), float(panel_size))
            painter.fillRect(rect, QColor("#FFFFFF"))
            painter.setPen(QPen(QColor("#D8D8D8"), 1.0))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawRect(rect.adjusted(0.0, 0.0, -1.0, -1.0))
            anchor0, anchor1 = self.canonical_anchors(spec, rect.adjusted(12.0, 12.0, -12.0, -12.0))
            self.paint_glyph(painter, spec, anchor0, anchor1)
        painter.end()
        return image


_RENDERER = OrbitalRenderer(glyph_library=ORBITAL_GLYPH_LIBRARY, family_presets=_ACTIVE_ORBITAL_FAMILY_PRESETS)


def reload_orbital_presets(path: Path = ORBITAL_PRESET_CONFIG_PATH) -> dict[str, OrbitalFamilyPreset]:
    """Recarga presets desde disco y actualiza el renderer compartido."""
    global _ACTIVE_ORBITAL_FAMILY_PRESETS, ORBITAL_FAMILY_PRESETS, ORBITAL_GEOMETRY_PRESETS, ORBITAL_GLYPH_LIBRARY
    _ACTIVE_ORBITAL_FAMILY_PRESETS = load_orbital_family_presets(path)
    ORBITAL_FAMILY_PRESETS = _ACTIVE_ORBITAL_FAMILY_PRESETS
    ORBITAL_GEOMETRY_PRESETS = ORBITAL_FAMILY_PRESETS
    ORBITAL_GLYPH_LIBRARY = _build_glyph_library_from_presets(_ACTIVE_ORBITAL_FAMILY_PRESETS)
    _RENDERER.set_glyph_library(ORBITAL_GLYPH_LIBRARY, _ACTIVE_ORBITAL_FAMILY_PRESETS)
    return _ACTIVE_ORBITAL_FAMILY_PRESETS


def draw_orbital_icon(kind: str, size: int = ORBITAL_ICON_SIZE) -> QIcon:
    """Icono Qt para toolbar/paleta."""
    return _RENDERER.draw_icon(kind, size=size)


def render_orbital_palette_image(model: OrbitalPaletteModel = ORBITAL_PALETTE_MODEL) -> QImage:
    """Renderiza el grid orbital completo usando los presets activos."""
    return _RENDERER.render_palette_image(model)


def render_orbital_family_triptych(family: str, panel_size: int = 180, gap: int = 18) -> QImage:
    """Preview grande outline/shaded/solid para una familia orbital."""
    return _RENDERER.render_style_triptych(family, panel_size=panel_size, gap=gap)


def orbital_renderer() -> OrbitalRenderer:
    """Devuelve la instancia compartida del renderer orbital."""
    return _RENDERER
