"""Renderer vectorial para una paleta orbital tipo ChemDraw."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math

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

_GRADIENT_STOP_PRESETS: dict[str, tuple[tuple[float, str], ...]] = {
    "default": (
        (0.0, "#FFFFFF"),
        (0.38, "#D9D9D9"),
        (0.72, "#676767"),
        (1.0, "#101010"),
    ),
    "dual_positive": (
        (0.0, "#FFFFFF"),
        (0.34, "#E8E8E8"),
        (0.72, "#747474"),
        (1.0, "#111111"),
    ),
    "dual_negative": (
        (0.0, "#202020"),
        (0.34, "#5A5A5A"),
        (0.74, "#BFBFBF"),
        (1.0, "#F4F4F4"),
    ),
    "inverted": (
        (0.0, "#101010"),
        (0.34, "#676767"),
        (0.72, "#D9D9D9"),
        (1.0, "#FFFFFF"),
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


class OrbitalType(str, Enum):
    """Tipos orbitales soportados por la UI."""

    S = "s"
    P = "p"
    D = "d"
    DZ2 = "dz2"
    F = "f"
    FZ3 = "fz3"
    SP3 = "sp3"
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
    gradient: str = "linear"
    stroke: bool = True
    phase: str | None = None
    light_dir: str | None = None
    light_variant: str | None = None


@dataclass(frozen=True)
class GlyphDefinition:
    """Glifo manual trazado en coordenadas normalizadas."""

    id: str
    paths_outline: tuple[GlyphLayer, ...]
    paths_shaded: tuple[GlyphLayer, ...]
    paths_solid: tuple[GlyphLayer, ...]
    default_bounds: QRectF
    anchor_center: QPointF
    anchor_bias: QPointF
    visual_padding: float


@dataclass(frozen=True)
class LobeProfile:
    """Coeficientes normalizados para un lóbulo cúbico simétrico."""

    tip_ctrl_x: float
    tip_ctrl_y: float
    shoulder_ctrl_x: float
    shoulder_ctrl_y: float
    waist_end_x: float
    waist_end_y: float
    neck_ctrl_x: float
    neck_ctrl_y: float
    neck_inner_y: float
    end_y: float
    neck_inner_x: float | None = None


@dataclass(frozen=True)
class OrbitalGeometryPreset:
    """Preset geométrico reusable para construir un glifo orbital."""

    main_lobe_width: float
    main_lobe_height: float
    secondary_lobe_width: float
    secondary_lobe_height: float
    neck_ratio: float
    lobe_separation: float
    torus_thickness: float
    upper_lobe_offset: float
    lower_lobe_offset: float
    visual_padding: float
    anchor_bias: tuple[float, float] = (0.0, 0.0)
    default_bounds_override: tuple[float, float, float, float] | None = None
    main_profile: str = "large"
    secondary_profile: str | None = None
    outline_main_profile: str | None = None
    outline_secondary_profile: str | None = None
    outline_main_scale: tuple[float, float] = (1.0, 1.0)
    outline_secondary_scale: tuple[float, float] = (1.0, 1.0)
    torus_outer_width: float = 0.0
    torus_outer_height: float = 0.0
    torus_inner_scale: tuple[float, float] = (0.0, 0.0)
    outline_torus_scale: tuple[float, float] = (1.0, 1.0)
    outline_torus_inner_scale: tuple[float, float] | None = None


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
    OrbitalType.SP3: "Orbital sp3",
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
    OrbitalType.S: "s_round",
    OrbitalType.P: "p_vertical",
    OrbitalType.D: "d_clover",
    OrbitalType.DZ2: "dz2_axial",
    OrbitalType.F: "f_clover",
    OrbitalType.FZ3: "fz3_axial",
    OrbitalType.SP3: "sp3_hybrid",
    OrbitalType.SP_LOBE: "sp_lobe",
    OrbitalType.SIGMA_BONDING: "sigma_bonding",
    OrbitalType.PI_BONDING: "pi_bonding",
    OrbitalType.TORUS: "torus_flat",
}

_CANVAS_EXTENT_SCALES = {
    OrbitalType.S: 0.84,
    OrbitalType.P: 0.92,
    OrbitalType.D: 0.92,
    OrbitalType.DZ2: 0.98,
    OrbitalType.F: 0.94,
    OrbitalType.FZ3: 0.98,
    OrbitalType.SP3: 0.94,
    OrbitalType.SP_LOBE: 0.88,
    OrbitalType.SIGMA_BONDING: 0.96,
    OrbitalType.PI_BONDING: 0.98,
    OrbitalType.TORUS: 0.84,
}

ORBITAL_GEOMETRY_PRESETS: dict[str, OrbitalGeometryPreset] = {
    "s_round": OrbitalGeometryPreset(
        main_lobe_width=68.0,
        main_lobe_height=68.0,
        secondary_lobe_width=68.0,
        secondary_lobe_height=68.0,
        neck_ratio=1.0,
        lobe_separation=0.0,
        torus_thickness=0.0,
        upper_lobe_offset=0.0,
        lower_lobe_offset=0.0,
        visual_padding=6.0,
        default_bounds_override=(-34.0, -34.0, 68.0, 68.0),
    ),
    "torus_flat": OrbitalGeometryPreset(
        main_lobe_width=72.0,
        main_lobe_height=20.0,
        secondary_lobe_width=28.8,
        secondary_lobe_height=6.8,
        neck_ratio=0.0,
        lobe_separation=0.0,
        torus_thickness=6.6,
        upper_lobe_offset=0.0,
        lower_lobe_offset=0.0,
        visual_padding=6.0,
        default_bounds_override=(-36.0, -10.0, 72.0, 20.0),
        torus_outer_width=72.0,
        torus_outer_height=20.0,
        torus_inner_scale=(0.40, 0.34),
    ),
    "p_vertical": OrbitalGeometryPreset(
        main_lobe_width=30.4,
        main_lobe_height=41.0,
        secondary_lobe_width=30.4,
        secondary_lobe_height=41.0,
        neck_ratio=0.3552631579,
        lobe_separation=0.0,
        torus_thickness=0.0,
        upper_lobe_offset=0.0,
        lower_lobe_offset=0.0,
        visual_padding=4.5,
        default_bounds_override=(-16.0, -41.0, 32.0, 82.0),
        main_profile="p_refined",
        secondary_profile="p_refined",
        outline_main_profile="p_outline",
        outline_secondary_profile="p_outline",
        outline_main_scale=(1.3157894737, 0.9268292683),
        outline_secondary_scale=(1.3157894737, 0.9268292683),
    ),
    "d_clover": OrbitalGeometryPreset(
        main_lobe_width=29.6,
        main_lobe_height=38.0,
        secondary_lobe_width=20.4,
        secondary_lobe_height=36.0,
        neck_ratio=0.3552631579,
        lobe_separation=0.0,
        torus_thickness=0.0,
        upper_lobe_offset=0.0,
        lower_lobe_offset=0.0,
        visual_padding=5.0,
        default_bounds_override=(-36.0, -38.0, 72.0, 76.0),
        main_profile="d_refined_vertical",
        secondary_profile="d_refined_horizontal",
        outline_main_profile="d_outline_vertical",
        outline_secondary_profile="d_outline_horizontal",
        outline_main_scale=(1.2162162162, 0.9210526316),
        outline_secondary_scale=(1.2745098039, 0.8888888889),
    ),
    "dz2_axial": OrbitalGeometryPreset(
        main_lobe_width=24.4,
        main_lobe_height=39.0,
        secondary_lobe_width=24.4,
        secondary_lobe_height=39.0,
        neck_ratio=0.3442622951,
        lobe_separation=0.0,
        torus_thickness=3.552,
        upper_lobe_offset=0.0,
        lower_lobe_offset=0.0,
        visual_padding=6.0,
        anchor_bias=(0.0, -0.75),
        default_bounds_override=(-31.0, -39.0, 62.0, 78.0),
        main_profile="dz2_refined",
        secondary_profile="dz2_refined",
        outline_main_profile="dz2_outline",
        outline_secondary_profile="dz2_outline",
        outline_main_scale=(1.0655737705, 0.9487179487),
        outline_secondary_scale=(1.0655737705, 0.9487179487),
        torus_outer_width=62.0,
        torus_outer_height=9.6,
        torus_inner_scale=(0.46, 0.26),
        outline_torus_scale=(1.0645161290, 1.4583333333),
        outline_torus_inner_scale=(0.44, 0.36),
    ),
    "sp3_hybrid": OrbitalGeometryPreset(
        main_lobe_width=48.0,
        main_lobe_height=48.0,
        secondary_lobe_width=24.0,
        secondary_lobe_height=28.0,
        neck_ratio=0.34,
        lobe_separation=23.5,
        torus_thickness=0.0,
        upper_lobe_offset=-22.0,
        lower_lobe_offset=1.5,
        visual_padding=6.0,
        default_bounds_override=(-21.0, -50.0, 42.0, 101.0),
        main_profile="large",
        secondary_profile="small",
    ),
    "sp_lobe": OrbitalGeometryPreset(
        main_lobe_width=29.6,
        main_lobe_height=44.0,
        secondary_lobe_width=11.4,
        secondary_lobe_height=39.0,
        neck_ratio=0.2982456140,
        lobe_separation=0.0,
        torus_thickness=0.0,
        upper_lobe_offset=0.0,
        lower_lobe_offset=0.0,
        visual_padding=4.0,
        anchor_bias=(0.0, -4.5),
        default_bounds_override=(-17.0, -38.0, 34.0, 76.0),
        main_profile="sp_refined_large",
        secondary_profile="sp_refined_small",
        outline_main_profile="sp_outline",
        outline_main_scale=(1.1486486486, 0.8636363636),
    ),
    "sigma_bonding": OrbitalGeometryPreset(
        main_lobe_width=32.0,
        main_lobe_height=38.0,
        secondary_lobe_width=13.0,
        secondary_lobe_height=12.0,
        neck_ratio=0.4375,
        lobe_separation=22.0,
        torus_thickness=0.0,
        upper_lobe_offset=0.0,
        lower_lobe_offset=22.0,
        visual_padding=5.0,
        anchor_bias=(0.0, 2.75),
        default_bounds_override=(-16.0, -38.0, 32.0, 72.0),
        main_profile="sigma_large",
        secondary_profile="sigma_small",
    ),
    "pi_bonding": OrbitalGeometryPreset(
        main_lobe_width=26.0,
        main_lobe_height=37.0,
        secondary_lobe_width=26.0,
        secondary_lobe_height=37.0,
        neck_ratio=0.3076923077,
        lobe_separation=4.0,
        torus_thickness=3.3,
        upper_lobe_offset=-2.0,
        lower_lobe_offset=2.0,
        visual_padding=6.0,
        default_bounds_override=(-24.0, -39.0, 48.0, 78.0),
        main_profile="dz2_outline",
        secondary_profile="dz2_outline",
        torus_outer_width=48.0,
        torus_outer_height=11.0,
        torus_inner_scale=(0.44, 0.40),
    ),
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

# Correspondencia visual de la referencia del usuario:
# fila 1, col 1 -> s_outline
# fila 1, col 2 -> s_shaded
# fila 1, col 3 -> p_shaded
# fila 1, col 4 -> p_solid
# fila 1, col 5 -> d_shaded
# fila 1, col 6 -> d_solid
# fila 1, col 7 -> vacío en la referencia recortada
# fila 2, col 1 -> torus_outline
# fila 2, col 2 -> torus_shaded
# fila 2, col 3 -> torus_solid
# fila 2, col 4 -> sp3_shaded
# fila 2, col 5 -> sp3_solid
# fila 2, col 6 -> dz2_shaded
# fila 2, col 7 -> dz2_solid
# fila 3, col 1 -> sp_lobe_outline
# fila 3, col 2 -> sp_lobe_shaded
# fila 3, col 3 -> sp_lobe_solid
# fila 3, col 4 -> sigma_bonding_shaded
# fila 3, col 5 -> sigma_bonding_solid
# fila 3, col 6 -> pi_bonding_shaded
# fila 3, col 7 -> pi_bonding_solid
REFERENCE_PALETTE_CELLS = (
    ("s_outline", "s_shaded", "p_shaded", "p_solid", "d_shaded", "d_solid", None),
    ("torus_outline", "torus_shaded", "torus_solid", "sp3_shaded", "sp3_solid", "dz2_shaded", "dz2_solid"),
    ("sp_lobe_outline", "sp_lobe_shaded", "sp_lobe_solid", "sigma_bonding_shaded", "sigma_bonding_solid", "pi_bonding_shaded", "pi_bonding_solid"),
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


def _ellipse_path(x: float, y: float, w: float, h: float) -> QPainterPath:
    path = QPainterPath()
    path.addEllipse(QRectF(x, y, w, h))
    return path


def _ring_path(x: float, y: float, w: float, h: float, inner_scale_x: float, inner_scale_y: float) -> QPainterPath:
    path = QPainterPath()
    path.setFillRule(Qt.FillRule.OddEvenFill)
    path.addEllipse(QRectF(x, y, w, h))
    inset_x = w * (1.0 - inner_scale_x) * 0.5
    inset_y = h * (1.0 - inner_scale_y) * 0.5
    path.addEllipse(QRectF(x + inset_x, y + inset_y, w * inner_scale_x, h * inner_scale_y))
    return path


def _closed_cubic_path(
    start: tuple[float, float],
    segments: tuple[tuple[tuple[float, float], tuple[float, float], tuple[float, float]], ...],
) -> QPainterPath:
    path = QPainterPath(QPointF(*start))
    for c1, c2, end in segments:
        path.cubicTo(QPointF(*c1), QPointF(*c2), QPointF(*end))
    path.closeSubpath()
    return path


def _transform_path(path: QPainterPath, *, translate: tuple[float, float] = (0.0, 0.0), rotate: float = 0.0, scale_x: float = 1.0, scale_y: float = 1.0) -> QPainterPath:
    transform = QTransform()
    transform.translate(translate[0], translate[1])
    if rotate:
        transform.rotate(rotate)
    transform.scale(scale_x, scale_y)
    return transform.map(path)


_LOBE_PROFILES: dict[str, LobeProfile] = {
    "large": LobeProfile(0.625, 0.96875, 1.0, 0.5625, 0.7916666667, 0.2916666667, 0.625, 0.09375, -0.0104166667, -0.025, neck_inner_x=0.3541666667),
    "small": LobeProfile(0.6666666667, 0.9285714286, 1.0, 0.5714285714, 0.8333333333, 0.3214285714, 0.6666666667, 0.1428571429, 0.0357142857, 0.0, neck_inner_x=0.3333333333),
    "micro": LobeProfile(0.7142857143, 0.9230769231, 1.0, 0.6153846154, 0.8571428571, 0.3076923077, 0.7142857143, 0.1384615385, 0.0153846154, 0.0, neck_inner_x=0.3714285714),
    "horizontal_petal": LobeProfile(0.5625, 0.9, 1.0, 0.525, 0.78125, 0.225, 0.5625, 0.05, 0.0, 0.01, neck_inner_x=0.2625),
    "diag_petal": LobeProfile(0.6666666667, 0.9545454545, 1.0, 0.6363636364, 0.8888888889, 0.3636363636, 0.6666666667, 0.1363636364, 0.0136363636, 0.0, neck_inner_x=0.2777777778),
    "p_outline": LobeProfile(0.65, 0.9473684211, 1.0, 0.5, 0.8, 0.1842105263, 0.65, -0.0526315789, -0.1973684211, -0.2368421053, neck_inner_x=0.35),
    "p_refined": LobeProfile(0.7105263158, 0.9634146341, 1.0, 0.7073170732, 0.8947368421, 0.4756097561, 0.7368421053, 0.2195121951, 0.0634146341, 0.0, neck_inner_x=0.3552631579),
    "d_outline_vertical": LobeProfile(0.6111111111, 0.9428571429, 1.0, 0.5142857143, 0.7777777778, 0.2285714286, 0.6666666667, 0.0285714286, -0.1, -0.1428571429, neck_inner_x=0.3333333333),
    "d_outline_horizontal": LobeProfile(0.6153846154, 0.9375, 1.0, 0.5625, 0.7692307692, 0.21875, 0.6153846154, 0.03125, -0.03125, -0.03125, neck_inner_x=0.3076923077),
    "d_refined_vertical": LobeProfile(0.6756756757, 0.9526315789, 1.0, 0.6710526316, 0.8783783784, 0.4105263158, 0.7297297297, 0.1894736842, 0.0473684211, -0.0105263158, neck_inner_x=0.3513513514),
    "d_refined_horizontal": LobeProfile(0.6078431373, 0.9305555556, 1.0, 0.5833333333, 0.8039215686, 0.2444444444, 0.6470588235, 0.0777777778, 0.0166666667, 0.0, neck_inner_x=0.2745098039),
    "sp_outline": LobeProfile(0.7647058824, 0.9210526316, 1.0, 0.4736842105, 0.8823529412, 0.0789473684, 0.7647058824, -0.3421052632, -0.6842105263, -1.0, neck_inner_x=0.4411764706),
    "sp_refined_small": LobeProfile(0.7192982456, 0.9743589744, 1.0, 0.8358974359, 0.8421052632, 0.6871794872, 0.6666666667, 0.5538461538, 0.4410256410, 0.3384615385, neck_inner_x=0.2982456140),
    "sp_refined_large": LobeProfile(0.3513513514, 0.9045454545, 0.7162162162, 0.75, 0.8513513514, 0.4886363636, 1.0, 0.1318181818, -0.2045454545, -0.2454545455, neck_inner_x=0.7837837838),
    "sigma_large": LobeProfile(0.75, 0.9210526316, 1.0, 0.5263157895, 0.875, 0.2105263158, 0.75, -0.0789473684, -0.3947368421, -0.7368421053, neck_inner_x=0.4375),
    "sigma_small": LobeProfile(0.3846153846, 0.9166666667, 0.7692307692, 0.75, 0.9230769231, 0.5, 1.0, 0.1666666667, -0.0833333333, 0.0, neck_inner_x=0.6923076923),
    "dz2_outline": LobeProfile(0.6923076923, 0.9459459459, 1.0, 0.6216216216, 0.8461538462, 0.3783783784, 0.6923076923, 0.1621621622, 0.0270270270, 0.0, neck_inner_x=0.3076923077),
    "dz2_refined": LobeProfile(0.6721311475, 0.9589743590, 1.0, 0.6666666667, 0.8524590164, 0.4358974359, 0.7213114754, 0.2564102564, 0.1589743590, 0.1435897436, neck_inner_x=0.3442622951),
}


def _scale_dimensions(width: float, height: float, scale: tuple[float, float]) -> tuple[float, float]:
    return width * scale[0], height * scale[1]


def _profile(name: str) -> LobeProfile:
    return _LOBE_PROFILES[name]


def _vertical_lobe_path(width: float, height: float, neck_ratio: float, profile_name: str) -> QPainterPath:
    profile = _profile(profile_name)
    half_width = width * 0.5
    inner_x_ratio = profile.neck_inner_x if profile.neck_inner_x is not None else neck_ratio
    return _closed_cubic_path(
        (0.0, -height),
        (
            ((half_width * profile.tip_ctrl_x, -height * profile.tip_ctrl_y), (half_width * profile.shoulder_ctrl_x, -height * profile.shoulder_ctrl_y), (half_width * profile.waist_end_x, -height * profile.waist_end_y)),
            ((half_width * profile.neck_ctrl_x, -height * profile.neck_ctrl_y), (half_width * inner_x_ratio, -height * profile.neck_inner_y), (0.0, -height * profile.end_y)),
            ((-half_width * inner_x_ratio, -height * profile.neck_inner_y), (-half_width * profile.neck_ctrl_x, -height * profile.neck_ctrl_y), (-half_width * profile.waist_end_x, -height * profile.waist_end_y)),
            ((-half_width * profile.shoulder_ctrl_x, -height * profile.shoulder_ctrl_y), (-half_width * profile.tip_ctrl_x, -height * profile.tip_ctrl_y), (0.0, -height)),
        ),
    )


def _oriented_lobe_path(
    width: float,
    height: float,
    neck_ratio: float,
    profile_name: str,
    *,
    orientation: str = "up",
    offset: float = 0.0,
) -> QPainterPath:
    base = _vertical_lobe_path(width, height, neck_ratio, profile_name)
    if orientation == "up":
        return _transform_path(base, translate=(0.0, offset))
    if orientation == "down":
        return _transform_path(base, scale_y=-1.0, translate=(0.0, offset))
    if orientation == "left":
        return _transform_path(base, rotate=-90.0, translate=(offset, 0.0))
    if orientation == "right":
        return _transform_path(base, rotate=90.0, translate=(offset, 0.0))
    raise ValueError(f"Orientación de lóbulo no soportada: {orientation}")


def _ring_from_metrics(width: float, height: float, inner_scale: tuple[float, float], *, offset_y: float = 0.0) -> QPainterPath:
    return _ring_path(-width * 0.5, offset_y - height * 0.5, width, height, inner_scale[0], inner_scale[1])


def _vertical_lobe_large() -> QPainterPath:
    return _vertical_lobe_path(48.0, 48.0, 0.3541666667, "large")


def _vertical_lobe_small() -> QPainterPath:
    return _vertical_lobe_path(24.0, 28.0, 0.3333333333, "small")


def _vertical_lobe_micro() -> QPainterPath:
    return _vertical_lobe_path(14.0, 13.0, 0.3714285714, "micro")


def _horizontal_petal() -> QPainterPath:
    return _oriented_lobe_path(32.0, 40.0, 0.2625, "horizontal_petal", orientation="left")


def _diag_petal_upper_left() -> QPainterPath:
    base = _vertical_lobe_path(18.0, 22.0, 0.2777777778, "diag_petal")
    return _transform_path(base, rotate=-38.0, translate=(-20.0, -9.0))


def _diag_petal_lower_right() -> QPainterPath:
    return _transform_path(_diag_petal_upper_left(), scale_x=-1.0, scale_y=-1.0)


def _bounds_for_layers(*groups: tuple[GlyphLayer, ...]) -> QRectF:
    bounds: QRectF | None = None
    for group in groups:
        for layer in group:
            rect = layer.path.boundingRect()
            bounds = rect if bounds is None else bounds.united(rect)
    return bounds if bounds is not None else QRectF(-1.0, -1.0, 2.0, 2.0)


def _glyph_from_layers(
    glyph_id: str,
    preset: OrbitalGeometryPreset,
    *,
    paths_outline: tuple[GlyphLayer, ...],
    paths_shaded: tuple[GlyphLayer, ...],
    paths_solid: tuple[GlyphLayer, ...],
) -> GlyphDefinition:
    default_bounds = (
        QRectF(*preset.default_bounds_override)
        if preset.default_bounds_override is not None
        else _bounds_for_layers(paths_outline, paths_shaded, paths_solid)
    )
    return GlyphDefinition(
        id=glyph_id,
        paths_outline=paths_outline,
        paths_shaded=paths_shaded,
        paths_solid=paths_solid,
        default_bounds=default_bounds,
        anchor_center=QPointF(0.0, 0.0),
        anchor_bias=QPointF(*preset.anchor_bias),
        visual_padding=preset.visual_padding,
    )


def _preset_ring_path(preset: OrbitalGeometryPreset, *, outline_variant: bool = False) -> QPainterPath:
    width = preset.torus_outer_width
    height = preset.torus_outer_height
    inner_scale = preset.torus_inner_scale
    if outline_variant:
        width, height = _scale_dimensions(width, height, preset.outline_torus_scale)
        inner_scale = preset.outline_torus_inner_scale or inner_scale
    return _ring_from_metrics(width, height, inner_scale)


def _phase_shaded_layer(
    path: QPainterPath,
    *,
    phase: str,
    light_dir: str,
    gradient: str = "linear",
    stroke: bool = False,
    light_variant: str | None = None,
) -> GlyphLayer:
    return GlyphLayer(
        path,
        "shaded",
        gradient=gradient,
        stroke=stroke,
        phase=phase,
        light_dir=light_dir,
        light_variant=light_variant,
    )


def _phase_solid_layer(
    path: QPainterPath,
    *,
    phase: str,
    light_dir: str | None = None,
    stroke: bool = False,
) -> GlyphLayer:
    return GlyphLayer(
        path,
        "solid",
        stroke=stroke,
        phase=phase,
        light_dir=light_dir,
    )


def _build_s_round_glyph() -> GlyphDefinition:
    preset = ORBITAL_GEOMETRY_PRESETS["s_round"]
    sphere = _ellipse_path(
        -preset.main_lobe_width * 0.5,
        -preset.main_lobe_height * 0.5,
        preset.main_lobe_width,
        preset.main_lobe_height,
    )
    return _glyph_from_layers(
        "s_round",
        preset,
        paths_outline=(GlyphLayer(sphere, "outline"),),
        paths_shaded=(
            _phase_shaded_layer(
                sphere,
                phase="positive",
                light_dir="southeast",
                gradient="radial",
                stroke=False,
            ),
        ),
        paths_solid=(
            _phase_solid_layer(
                sphere,
                phase="positive",
                stroke=True,
            ),
        ),
    )


def _build_torus_flat_glyph() -> GlyphDefinition:
    preset = ORBITAL_GEOMETRY_PRESETS["torus_flat"]
    torus = _preset_ring_path(preset)
    return _glyph_from_layers(
        "torus_flat",
        preset,
        paths_outline=(GlyphLayer(torus, "outline"),),
        paths_shaded=(
            _phase_shaded_layer(
                torus,
                phase="positive",
                light_dir="south",
                gradient="elliptical",
            ),
        ),
        paths_solid=(
            _phase_solid_layer(
                torus,
                phase="positive",
            ),
        ),
    )


def _build_p_vertical_glyph() -> GlyphDefinition:
    preset = ORBITAL_GEOMETRY_PRESETS["p_vertical"]
    outline_w, outline_h = _scale_dimensions(preset.main_lobe_width, preset.main_lobe_height, preset.outline_main_scale)
    outline_profile = preset.outline_main_profile or preset.main_profile
    top_outline = _oriented_lobe_path(outline_w, outline_h, preset.neck_ratio, outline_profile, orientation="up", offset=preset.upper_lobe_offset)
    bottom_outline = _oriented_lobe_path(outline_w, outline_h, preset.neck_ratio, outline_profile, orientation="down", offset=preset.lower_lobe_offset)
    top_fill = _oriented_lobe_path(preset.main_lobe_width, preset.main_lobe_height, preset.neck_ratio, preset.main_profile, orientation="up", offset=preset.upper_lobe_offset)
    bottom_fill = _oriented_lobe_path(preset.secondary_lobe_width, preset.secondary_lobe_height, preset.neck_ratio, preset.secondary_profile or preset.main_profile, orientation="down", offset=preset.lower_lobe_offset)
    return _glyph_from_layers(
        "p_vertical",
        preset,
        paths_outline=(
            GlyphLayer(top_outline, "outline"),
            GlyphLayer(bottom_outline, "outline"),
        ),
        paths_shaded=(
            _phase_shaded_layer(top_fill, phase="positive", light_dir="south"),
            _phase_shaded_layer(bottom_fill, phase="negative", light_dir="north"),
        ),
        paths_solid=(
            _phase_solid_layer(top_fill, phase="positive"),
            _phase_solid_layer(bottom_fill, phase="negative"),
        ),
    )


def _build_d_clover_glyph() -> GlyphDefinition:
    preset = ORBITAL_GEOMETRY_PRESETS["d_clover"]
    outline_main_w, outline_main_h = _scale_dimensions(preset.main_lobe_width, preset.main_lobe_height, preset.outline_main_scale)
    outline_secondary_w, outline_secondary_h = _scale_dimensions(preset.secondary_lobe_width, preset.secondary_lobe_height, preset.outline_secondary_scale)
    outline_top = _oriented_lobe_path(outline_main_w, outline_main_h, preset.neck_ratio, preset.outline_main_profile or preset.main_profile, orientation="up")
    outline_bottom = _oriented_lobe_path(outline_main_w, outline_main_h, preset.neck_ratio, preset.outline_main_profile or preset.main_profile, orientation="down")
    outline_left = _oriented_lobe_path(outline_secondary_w, outline_secondary_h, preset.neck_ratio, preset.outline_secondary_profile or preset.secondary_profile or preset.main_profile, orientation="left")
    outline_right = _oriented_lobe_path(outline_secondary_w, outline_secondary_h, preset.neck_ratio, preset.outline_secondary_profile or preset.secondary_profile or preset.main_profile, orientation="right")
    fill_top = _oriented_lobe_path(preset.main_lobe_width, preset.main_lobe_height, preset.neck_ratio, preset.main_profile, orientation="up")
    fill_bottom = _oriented_lobe_path(preset.main_lobe_width, preset.main_lobe_height, preset.neck_ratio, preset.main_profile, orientation="down")
    fill_left = _oriented_lobe_path(preset.secondary_lobe_width, preset.secondary_lobe_height, preset.neck_ratio, preset.secondary_profile or preset.main_profile, orientation="left")
    fill_right = _oriented_lobe_path(preset.secondary_lobe_width, preset.secondary_lobe_height, preset.neck_ratio, preset.secondary_profile or preset.main_profile, orientation="right")
    return _glyph_from_layers(
        "d_clover",
        preset,
        paths_outline=(
            GlyphLayer(outline_left, "outline"),
            GlyphLayer(outline_top, "outline"),
            GlyphLayer(outline_right, "outline"),
            GlyphLayer(outline_bottom, "outline"),
        ),
        paths_shaded=(
            _phase_shaded_layer(fill_left, phase="positive", light_dir="east"),
            _phase_shaded_layer(fill_top, phase="negative", light_dir="south"),
            _phase_shaded_layer(fill_right, phase="positive", light_dir="west"),
            _phase_shaded_layer(fill_bottom, phase="negative", light_dir="north"),
        ),
        paths_solid=(
            _phase_solid_layer(fill_left, phase="positive"),
            _phase_solid_layer(fill_top, phase="negative"),
            _phase_solid_layer(fill_right, phase="positive"),
            _phase_solid_layer(fill_bottom, phase="negative"),
        ),
    )


def _build_sp3_hybrid_glyph() -> GlyphDefinition:
    preset = ORBITAL_GEOMETRY_PRESETS["sp3_hybrid"]
    top = _oriented_lobe_path(preset.secondary_lobe_width, preset.secondary_lobe_height, preset.neck_ratio, preset.secondary_profile or preset.main_profile, orientation="up", offset=preset.upper_lobe_offset)
    bottom = _oriented_lobe_path(preset.main_lobe_width, preset.main_lobe_height, preset.neck_ratio, preset.main_profile, orientation="down", offset=preset.lower_lobe_offset)
    return _glyph_from_layers(
        "sp3_hybrid",
        preset,
        paths_outline=(
            GlyphLayer(top, "outline"),
            GlyphLayer(bottom, "outline"),
        ),
        paths_shaded=(
            _phase_shaded_layer(top, phase="negative", light_dir="south"),
            _phase_shaded_layer(bottom, phase="positive", light_dir="north"),
        ),
        paths_solid=(
            _phase_solid_layer(top, phase="negative"),
            _phase_solid_layer(bottom, phase="positive"),
        ),
    )


def _build_sp_lobe_glyph() -> GlyphDefinition:
    preset = ORBITAL_GEOMETRY_PRESETS["sp_lobe"]
    top = _oriented_lobe_path(preset.secondary_lobe_width, preset.secondary_lobe_height, preset.neck_ratio, preset.secondary_profile or preset.main_profile, orientation="up", offset=preset.upper_lobe_offset)
    bottom = _oriented_lobe_path(preset.main_lobe_width, preset.main_lobe_height, preset.neck_ratio, preset.main_profile, orientation="down", offset=preset.lower_lobe_offset)
    if preset.outline_main_profile:
        outline_w, outline_h = _scale_dimensions(preset.main_lobe_width, preset.main_lobe_height, preset.outline_main_scale)
        outline = _vertical_lobe_path(outline_w, outline_h, preset.neck_ratio, preset.outline_main_profile)
    else:
        outline = top.united(bottom)
    return _glyph_from_layers(
        "sp_lobe",
        preset,
        paths_outline=(GlyphLayer(outline, "outline"),),
        paths_shaded=(
            _phase_shaded_layer(top, phase="negative", light_dir="south"),
            _phase_shaded_layer(bottom, phase="positive", light_dir="north"),
        ),
        paths_solid=(
            _phase_solid_layer(top, phase="negative"),
            _phase_solid_layer(bottom, phase="positive"),
        ),
    )


def _build_sigma_bonding_glyph() -> GlyphDefinition:
    preset = ORBITAL_GEOMETRY_PRESETS["sigma_bonding"]
    large = _oriented_lobe_path(preset.main_lobe_width, preset.main_lobe_height, preset.neck_ratio, preset.main_profile, orientation="up", offset=preset.upper_lobe_offset)
    small = _oriented_lobe_path(preset.secondary_lobe_width, preset.secondary_lobe_height, preset.neck_ratio, preset.secondary_profile or preset.main_profile, orientation="down", offset=preset.lower_lobe_offset)
    return _glyph_from_layers(
        "sigma_bonding",
        preset,
        paths_outline=(
            GlyphLayer(large, "outline"),
            GlyphLayer(small, "outline"),
        ),
        paths_shaded=(
            _phase_shaded_layer(large, phase="positive", light_dir="south"),
            _phase_shaded_layer(small, phase="negative", light_dir="north"),
        ),
        paths_solid=(
            _phase_solid_layer(large, phase="positive"),
            _phase_solid_layer(small, phase="negative"),
        ),
    )


def _build_dz2_axial_glyph() -> GlyphDefinition:
    preset = ORBITAL_GEOMETRY_PRESETS["dz2_axial"]
    outline_w, outline_h = _scale_dimensions(preset.main_lobe_width, preset.main_lobe_height, preset.outline_main_scale)
    outline_profile = preset.outline_main_profile or preset.main_profile
    top_outline = _oriented_lobe_path(outline_w, outline_h, preset.neck_ratio, outline_profile, orientation="up")
    bottom_outline = _oriented_lobe_path(outline_w, outline_h, preset.neck_ratio, outline_profile, orientation="down")
    outline_ring = _preset_ring_path(preset, outline_variant=True)
    top_fill = _oriented_lobe_path(preset.main_lobe_width, preset.main_lobe_height, preset.neck_ratio, preset.main_profile, orientation="up")
    bottom_fill = _oriented_lobe_path(preset.secondary_lobe_width, preset.secondary_lobe_height, preset.neck_ratio, preset.secondary_profile or preset.main_profile, orientation="down")
    fill_ring = _preset_ring_path(preset)
    return _glyph_from_layers(
        "dz2_axial",
        preset,
        paths_outline=(
            GlyphLayer(top_outline, "outline"),
            GlyphLayer(outline_ring, "outline"),
            GlyphLayer(bottom_outline, "outline"),
        ),
        paths_shaded=(
            _phase_shaded_layer(top_fill, phase="positive", light_dir="south"),
            _phase_shaded_layer(fill_ring, phase="negative", light_dir="south", gradient="elliptical"),
            _phase_shaded_layer(bottom_fill, phase="positive", light_dir="north"),
        ),
        paths_solid=(
            _phase_solid_layer(top_fill, phase="positive"),
            _phase_solid_layer(fill_ring, phase="negative"),
            _phase_solid_layer(bottom_fill, phase="positive"),
        ),
    )


def _build_pi_bonding_glyph() -> GlyphDefinition:
    preset = ORBITAL_GEOMETRY_PRESETS["pi_bonding"]
    top = _oriented_lobe_path(preset.main_lobe_width, preset.main_lobe_height, preset.neck_ratio, preset.main_profile, orientation="up", offset=preset.upper_lobe_offset)
    bottom = _oriented_lobe_path(preset.secondary_lobe_width, preset.secondary_lobe_height, preset.neck_ratio, preset.secondary_profile or preset.main_profile, orientation="down", offset=preset.lower_lobe_offset)
    ring = _preset_ring_path(preset)
    return _glyph_from_layers(
        "pi_bonding",
        preset,
        paths_outline=(
            GlyphLayer(top, "outline"),
            GlyphLayer(ring, "outline"),
            GlyphLayer(bottom, "outline"),
        ),
        paths_shaded=(
            _phase_shaded_layer(top, phase="positive", light_dir="south"),
            _phase_shaded_layer(ring, phase="neutral", light_dir="south", gradient="elliptical"),
            _phase_shaded_layer(bottom, phase="negative", light_dir="north"),
        ),
        paths_solid=(
            _phase_solid_layer(top, phase="positive"),
            _phase_solid_layer(ring, phase="neutral"),
            _phase_solid_layer(bottom, phase="negative"),
        ),
    )


def _build_f_clover_glyph() -> GlyphDefinition:
    f_top = _transform_path(_vertical_lobe_small(), translate=(0.0, -20.0))
    f_bottom = _transform_path(f_top, scale_x=1.0, scale_y=-1.0)
    f_left = _transform_path(_horizontal_petal(), translate=(-4.0, 0.0), scale_x=0.64, scale_y=0.64)
    f_right = _transform_path(f_left, scale_x=-1.0, scale_y=1.0)
    f_diag_a = _diag_petal_upper_left()
    f_diag_b = _diag_petal_lower_right()
    outline = (
        GlyphLayer(f_left, "outline"),
        GlyphLayer(f_top, "outline"),
        GlyphLayer(f_diag_a, "outline"),
        GlyphLayer(f_right, "outline"),
        GlyphLayer(f_bottom, "outline"),
        GlyphLayer(f_diag_b, "outline"),
    )
    shaded = (
        GlyphLayer(f_left, "shaded", stroke=False),
        GlyphLayer(f_top, "shaded", stroke=False),
        GlyphLayer(f_diag_a, "shaded", stroke=False),
        GlyphLayer(f_right, "shaded", stroke=False),
        GlyphLayer(f_bottom, "shaded", stroke=False),
        GlyphLayer(f_diag_b, "shaded", stroke=False),
    )
    solid = (
        GlyphLayer(f_left, "solid", stroke=False),
        GlyphLayer(f_top, "solid", stroke=False),
        GlyphLayer(f_diag_a, "solid", stroke=False),
        GlyphLayer(f_right, "solid", stroke=False),
        GlyphLayer(f_bottom, "solid", stroke=False),
        GlyphLayer(f_diag_b, "solid", stroke=False),
    )
    return GlyphDefinition(
        id="f_clover",
        paths_outline=outline,
        paths_shaded=shaded,
        paths_solid=solid,
        default_bounds=_bounds_for_layers(outline, shaded, solid),
        anchor_center=QPointF(0.0, 0.0),
        anchor_bias=QPointF(0.0, 0.0),
        visual_padding=10.0,
    )


def _build_fz3_axial_glyph() -> GlyphDefinition:
    fz3_top = _transform_path(_vertical_lobe_small(), translate=(0.0, -28.0))
    fz3_mid = _transform_path(_vertical_lobe_micro(), translate=(0.0, 0.0))
    fz3_bottom = _transform_path(_vertical_lobe_small(), scale_x=1.0, scale_y=-1.0, translate=(0.0, 28.0))
    fz3_ring = _ring_path(-23.0, -6.0, 46.0, 12.0, 0.48, 0.46)
    outline = (
        GlyphLayer(fz3_top, "outline"),
        GlyphLayer(fz3_mid, "outline"),
        GlyphLayer(fz3_ring, "outline"),
        GlyphLayer(fz3_bottom, "outline"),
    )
    shaded = (
        GlyphLayer(fz3_top, "shaded", stroke=False),
        GlyphLayer(fz3_mid, "shaded", stroke=False),
        GlyphLayer(fz3_ring, "shaded", gradient="radial", stroke=False),
        GlyphLayer(fz3_bottom, "shaded", stroke=False),
    )
    solid = (
        GlyphLayer(fz3_top, "solid", stroke=False),
        GlyphLayer(fz3_mid, "solid", stroke=False),
        GlyphLayer(fz3_ring, "solid", stroke=False),
        GlyphLayer(fz3_bottom, "solid", stroke=False),
    )
    return GlyphDefinition(
        id="fz3_axial",
        paths_outline=outline,
        paths_shaded=shaded,
        paths_solid=solid,
        default_bounds=_bounds_for_layers(outline, shaded, solid),
        anchor_center=QPointF(0.0, 0.0),
        anchor_bias=QPointF(0.0, 0.0),
        visual_padding=10.0,
    )


def _build_glyph_library() -> dict[str, GlyphDefinition]:
    return {
        "s_round": _build_s_round_glyph(),
        "torus_flat": _build_torus_flat_glyph(),
        "p_vertical": _build_p_vertical_glyph(),
        "d_clover": _build_d_clover_glyph(),
        "sp3_hybrid": _build_sp3_hybrid_glyph(),
        "sp_lobe": _build_sp_lobe_glyph(),
        "sigma_bonding": _build_sigma_bonding_glyph(),
        "dz2_axial": _build_dz2_axial_glyph(),
        "pi_bonding": _build_pi_bonding_glyph(),
        "f_clover": _build_f_clover_glyph(),
        "fz3_axial": _build_fz3_axial_glyph(),
    }


ORBITAL_GLYPH_LIBRARY = _build_glyph_library()


class OrbitalRenderer:
    """Renderer compartido por toolbar y documento."""

    def spec_for_kind(self, kind: str) -> OrbitalGlyphSpec:
        return ORBITAL_SPECS.get(kind) or ORBITAL_SPECS[DEFAULT_ORBITAL_KIND]

    def glyph_for_spec(self, spec: OrbitalGlyphSpec) -> GlyphDefinition:
        return ORBITAL_GLYPH_LIBRARY[spec.glyph_id]

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
    ) -> tuple[GlyphLayer, ...]:
        glyph = self.glyph_for_spec(spec)
        transform = self._anchor_transform(glyph, anchor0, anchor1)
        transformed: list[GlyphLayer] = []
        for layer in self._layers_for_style(glyph, spec.style):
            transformed.append(
                GlyphLayer(
                    transform.map(layer.path),
                    layer.paint,
                    gradient=layer.gradient,
                    stroke=layer.stroke,
                    phase=layer.phase,
                    light_dir=layer.light_dir,
                    light_variant=layer.light_variant,
                )
            )
        return tuple(transformed)

    def combined_path(self, spec: OrbitalGlyphSpec, anchor0: QPointF, anchor1: QPointF) -> QPainterPath:
        path = QPainterPath()
        for layer in self.transformed_layers(spec, anchor0, anchor1):
            path.addPath(layer.path)
        return path

    def bounding_rect(self, spec: OrbitalGlyphSpec, anchor0: QPointF, anchor1: QPointF) -> QRectF:
        return self.combined_path(spec, anchor0, anchor1).boundingRect()

    def _stroke_width(self, glyph: GlyphDefinition, anchor0: QPointF, anchor1: QPointF) -> float:
        extent = math.hypot(anchor1.x() - anchor0.x(), anchor1.y() - anchor0.y())
        scale = extent / self._reference_extent(glyph)
        return max(1.1, min(2.2, scale * 1.9))

    @staticmethod
    def _gradient_stops_for_layer(layer: GlyphLayer) -> tuple[tuple[float, str], ...]:
        if layer.light_variant and layer.light_variant in _GRADIENT_STOP_PRESETS:
            return _GRADIENT_STOP_PRESETS[layer.light_variant]
        if layer.phase == "negative":
            return _GRADIENT_STOP_PRESETS["dual_negative"]
        if layer.phase == "positive":
            return _GRADIENT_STOP_PRESETS["dual_positive"]
        return _GRADIENT_STOP_PRESETS["default"]

    @staticmethod
    def _apply_gradient_stops(gradient: QLinearGradient | QRadialGradient, layer: GlyphLayer) -> None:
        for stop, color in OrbitalRenderer._gradient_stops_for_layer(layer):
            gradient.setColorAt(stop, QColor(color))

    @staticmethod
    def _light_vector(light_dir: str | None, rect: QRectF) -> tuple[float, float]:
        if light_dir and light_dir in _LIGHT_DIRECTION_VECTORS:
            return _LIGHT_DIRECTION_VECTORS[light_dir]
        if rect.height() >= rect.width():
            return _LIGHT_DIRECTION_VECTORS["south"]
        return _LIGHT_DIRECTION_VECTORS["east"]

    @staticmethod
    def _gradient_profile(rect: QRectF) -> dict[str, float | bool]:
        major = max(rect.width(), rect.height(), 1.0)
        minor = max(min(rect.width(), rect.height()), 1.0)
        is_small = major <= 24.0 or rect.width() * rect.height() <= 420.0 or minor <= 10.0
        if is_small:
            return {
                "small": True,
                "focus_bias": 0.16,
                "radius_scale": 0.84,
                "linear_span": 0.92,
            }
        return {
            "small": False,
            "focus_bias": 0.30,
            "radius_scale": 0.96,
            "linear_span": 1.10,
        }

    @staticmethod
    def _elliptical_brush_transform(rect: QRectF) -> QTransform:
        center = rect.center()
        width = max(rect.width(), 1.0)
        height = max(rect.height(), 1.0)
        transform = QTransform()
        transform.translate(center.x(), center.y())
        if width >= height:
            transform.scale(1.0, width / height)
        else:
            transform.scale(height / width, 1.0)
        transform.translate(-center.x(), -center.y())
        return transform

    @classmethod
    def _gradient_for_path(cls, layer: GlyphLayer) -> QBrush:
        rect = layer.path.boundingRect()
        profile = cls._gradient_profile(rect)
        vector_x, vector_y = cls._light_vector(layer.light_dir, rect)
        center = rect.center()

        if layer.gradient in {"radial", "elliptical"}:
            major = max(rect.width(), rect.height(), 1.0)
            radius = major * float(profile["radius_scale"])
            focus_bias = major * float(profile["focus_bias"])
            focus = QPointF(
                center.x() - vector_x * focus_bias,
                center.y() - vector_y * focus_bias,
            )
            gradient = QRadialGradient(center, radius, focus)
            cls._apply_gradient_stops(gradient, layer)
            brush = QBrush(gradient)
            if layer.gradient == "elliptical":
                brush.setTransform(cls._elliptical_brush_transform(rect))
            return brush

        span = float(profile["linear_span"])
        start = QPointF(
            center.x() - vector_x * rect.width() * 0.5 * span,
            center.y() - vector_y * rect.height() * 0.5 * span,
        )
        end = QPointF(
            center.x() + vector_x * rect.width() * 0.5 * span,
            center.y() + vector_y * rect.height() * 0.5 * span,
        )
        gradient = QLinearGradient(start, end)
        cls._apply_gradient_stops(gradient, layer)
        return QBrush(gradient)

    @staticmethod
    def _outline_pen(width: float) -> QPen:
        return QPen(
            _STROKE_COLOR,
            width,
            Qt.PenStyle.SolidLine,
            Qt.PenCapStyle.RoundCap,
            Qt.PenJoinStyle.RoundJoin,
        )

    def _paint_layer(
        self,
        painter: QPainter,
        layer: GlyphLayer,
        stroke_width: float,
        *,
        stroke_shaded_lobes: bool,
    ) -> None:
        if layer.paint == "outline":
            painter.setPen(self._outline_pen(stroke_width))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawPath(layer.path)
            return

        if layer.paint == "shaded":
            painter.setBrush(self._gradient_for_path(layer))
            if layer.stroke and stroke_shaded_lobes:
                painter.setPen(self._outline_pen(stroke_width))
            elif layer.stroke:
                painter.setPen(self._outline_pen(max(1.0, stroke_width * 0.72)))
            else:
                painter.setPen(Qt.PenStyle.NoPen)
            painter.drawPath(layer.path)
            return

        painter.setBrush(QBrush(_SOLID_FILL))
        if layer.stroke:
            painter.setPen(self._outline_pen(max(1.0, stroke_width * 0.92)))
        else:
            painter.setPen(Qt.PenStyle.NoPen)
        painter.drawPath(layer.path)

    def paint_glyph(
        self,
        painter: QPainter,
        spec: OrbitalGlyphSpec,
        anchor0: QPointF,
        anchor1: QPointF,
        *,
        stroke_shaded_lobes: bool | None = None,
    ) -> None:
        glyph = self.glyph_for_spec(spec)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
        stroke_width = self._stroke_width(glyph, anchor0, anchor1)
        stroke_flag = spec.stroke_shaded_lobes if stroke_shaded_lobes is None else bool(stroke_shaded_lobes)
        for layer in self.transformed_layers(spec, anchor0, anchor1):
            self._paint_layer(painter, layer, stroke_width, stroke_shaded_lobes=stroke_flag)

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
        target_center = inner.center() + QPointF(
            glyph.anchor_bias.x() * scale,
            glyph.anchor_bias.y() * scale,
        )
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


_RENDERER = OrbitalRenderer()


def draw_orbital_icon(kind: str, size: int = ORBITAL_ICON_SIZE) -> QIcon:
    """Icono Qt para toolbar/paleta."""
    return _RENDERER.draw_icon(kind, size=size)


def render_orbital_palette_image(model: OrbitalPaletteModel = ORBITAL_PALETTE_MODEL) -> QImage:
    """Renderiza el grid orbital completo."""
    return _RENDERER.render_palette_image(model)


def orbital_renderer() -> OrbitalRenderer:
    """Devuelve la instancia compartida del renderer orbital."""
    return _RENDERER
