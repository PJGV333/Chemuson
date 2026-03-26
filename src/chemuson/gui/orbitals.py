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


@dataclass(frozen=True)
class GlyphDefinition:
    """Glifo manual trazado en coordenadas normalizadas."""

    id: str
    paths_outline: tuple[GlyphLayer, ...]
    paths_shaded: tuple[GlyphLayer, ...]
    paths_solid: tuple[GlyphLayer, ...]
    default_bounds: QRectF
    anchor_center: QPointF
    visual_padding: float


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


def _label_for(orbital_type: OrbitalType, style: OrbitalStyle) -> str:
    return f"{_TYPE_LABELS[orbital_type]} ({_STYLE_LABELS[style]})"


def _build_all_specs() -> dict[str, OrbitalGlyphSpec]:
    specs: dict[str, OrbitalGlyphSpec] = {}
    for orbital_type in OrbitalType:
        for style in OrbitalStyle:
            kind = f"{orbital_type.value}_{style.value}"
            extent = 0.92
            if orbital_type in {OrbitalType.S, OrbitalType.TORUS}:
                extent = 0.84
            if orbital_type in {OrbitalType.DZ2, OrbitalType.FZ3, OrbitalType.PI_BONDING}:
                extent = 0.98
            specs[kind] = OrbitalGlyphSpec(
                kind=kind,
                orbital_type=orbital_type,
                glyph_id=_ORBITAL_GLYPH_IDS[orbital_type],
                style=style,
                label=_label_for(orbital_type, style),
                canvas_extent_scale=extent,
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


def _vertical_lobe_large() -> QPainterPath:
    return _closed_cubic_path(
        (0.0, -48.0),
        (
            ((15.0, -46.5), (24.0, -27.0), (19.0, -14.0)),
            ((15.0, -4.5), (8.5, 0.5), (0.0, 1.2)),
            ((-8.5, 0.5), (-15.0, -4.5), (-19.0, -14.0)),
            ((-24.0, -27.0), (-15.0, -46.5), (0.0, -48.0)),
        ),
    )


def _vertical_lobe_small() -> QPainterPath:
    return _closed_cubic_path(
        (0.0, -28.0),
        (
            ((8.0, -26.0), (12.0, -16.0), (10.0, -9.0)),
            ((8.0, -4.0), (4.0, -1.0), (0.0, 0.0)),
            ((-4.0, -1.0), (-8.0, -4.0), (-10.0, -9.0)),
            ((-12.0, -16.0), (-8.0, -26.0), (0.0, -28.0)),
        ),
    )


def _vertical_lobe_micro() -> QPainterPath:
    return _closed_cubic_path(
        (0.0, -13.0),
        (
            ((5.0, -12.0), (7.0, -8.0), (6.0, -4.0)),
            ((5.0, -1.8), (2.6, -0.2), (0.0, 0.0)),
            ((-2.6, -0.2), (-5.0, -1.8), (-6.0, -4.0)),
            ((-7.0, -8.0), (-5.0, -12.0), (0.0, -13.0)),
        ),
    )


def _horizontal_petal() -> QPainterPath:
    return _closed_cubic_path(
        (-40.0, 0.0),
        (
            ((-36.0, -9.0), (-21.0, -16.0), (-9.0, -12.5)),
            ((-2.0, -9.0), (0.0, -4.2), (-0.4, 0.0)),
            ((0.0, 4.2), (-2.0, 9.0), (-9.0, 12.5)),
            ((-21.0, 16.0), (-36.0, 9.0), (-40.0, 0.0)),
        ),
    )


def _diag_petal_upper_left() -> QPainterPath:
    base = _closed_cubic_path(
        (0.0, -22.0),
        (
            ((6.0, -21.0), (9.0, -14.0), (8.0, -8.0)),
            ((6.0, -3.0), (2.5, -0.3), (0.0, 0.0)),
            ((-2.5, -0.3), (-6.0, -3.0), (-8.0, -8.0)),
            ((-9.0, -14.0), (-6.0, -21.0), (0.0, -22.0)),
        ),
    )
    return _transform_path(base, rotate=-38.0, translate=(-20.0, -9.0))


def _diag_petal_lower_right() -> QPainterPath:
    return _transform_path(_diag_petal_upper_left(), scale_x=-1.0, scale_y=-1.0)


def _build_glyph_library() -> dict[str, GlyphDefinition]:
    circle_bounds = QRectF(-34.0, -34.0, 68.0, 68.0)
    sphere = _ellipse_path(-34.0, -34.0, 68.0, 68.0)
    torus = _ring_path(-36.0, -10.0, 72.0, 20.0, 0.40, 0.34)

    p_top = _closed_cubic_path(
        (0.0, -38.0),
        (
            ((13.0, -36.0), (20.0, -19.0), (16.0, -7.0)),
            ((13.0, 2.0), (7.0, 7.5), (0.0, 9.0)),
            ((-7.0, 7.5), (-13.0, 2.0), (-16.0, -7.0)),
            ((-20.0, -19.0), (-13.0, -36.0), (0.0, -38.0)),
        ),
    )
    p_bottom = _transform_path(p_top, scale_x=1.0, scale_y=-1.0)
    p_refined_top = _closed_cubic_path(
        (0.0, -41.0),
        (
            ((10.8, -39.5), (15.2, -29.0), (13.6, -19.5)),
            ((11.2, -9.0), (5.4, -2.6), (0.0, 0.0)),
            ((-5.4, -2.6), (-11.2, -9.0), (-13.6, -19.5)),
            ((-15.2, -29.0), (-10.8, -39.5), (0.0, -41.0)),
        ),
    )
    p_refined_bottom = _transform_path(p_refined_top, scale_x=1.0, scale_y=-1.0)

    d_top = _closed_cubic_path(
        (0.0, -35.0),
        (
            ((11.0, -33.0), (18.0, -18.0), (14.0, -8.0)),
            ((12.0, -1.0), (6.0, 3.5), (0.0, 5.0)),
            ((-6.0, 3.5), (-12.0, -1.0), (-14.0, -8.0)),
            ((-18.0, -18.0), (-11.0, -33.0), (0.0, -35.0)),
        ),
    )
    d_bottom = _transform_path(d_top, scale_x=1.0, scale_y=-1.0)
    d_left = _closed_cubic_path(
        (-32.0, 0.0),
        (
            ((-30.0, -8.0), (-18.0, -13.0), (-7.0, -10.0)),
            ((-1.0, -8.0), (1.0, -4.0), (1.0, 0.0)),
            ((1.0, 4.0), (-1.0, 8.0), (-7.0, 10.0)),
            ((-18.0, 13.0), (-30.0, 8.0), (-32.0, 0.0)),
        ),
    )
    d_right = _transform_path(d_left, scale_x=-1.0, scale_y=1.0)
    d_refined_top = _closed_cubic_path(
        (0.0, -38.0),
        (
            ((10.0, -36.2), (14.8, -25.5), (13.0, -15.6)),
            ((10.8, -7.2), (5.2, -1.8), (0.0, 0.4)),
            ((-5.2, -1.8), (-10.8, -7.2), (-13.0, -15.6)),
            ((-14.8, -25.5), (-10.0, -36.2), (0.0, -38.0)),
        ),
    )
    d_refined_bottom = _transform_path(d_refined_top, scale_x=1.0, scale_y=-1.0)
    d_refined_left = _closed_cubic_path(
        (-36.0, 0.0),
        (
            ((-33.5, -6.2), (-21.0, -10.2), (-8.8, -8.2)),
            ((-2.8, -6.6), (-0.6, -2.8), (0.0, 0.0)),
            ((-0.6, 2.8), (-2.8, 6.6), (-8.8, 8.2)),
            ((-21.0, 10.2), (-33.5, 6.2), (-36.0, 0.0)),
        ),
    )
    d_refined_right = _transform_path(d_refined_left, scale_x=-1.0, scale_y=1.0)

    sp_top = _transform_path(_vertical_lobe_small(), translate=(0.0, -22.0))
    sp_bottom = _transform_path(_vertical_lobe_large(), scale_x=1.0, scale_y=-1.0, translate=(0.0, 1.5))

    sp_single = _closed_cubic_path(
        (0.0, -38.0),
        (
            ((13.0, -35.0), (17.0, -18.0), (15.0, -3.0)),
            ((13.0, 13.0), (7.5, 26.0), (0.0, 38.0)),
            ((-7.5, 26.0), (-13.0, 13.0), (-15.0, -3.0)),
            ((-17.0, -18.0), (-13.0, -35.0), (0.0, -38.0)),
        ),
    )
    sp_refined_small = _closed_cubic_path(
        (0.0, -39.0),
        (
            ((4.1, -38.0), (5.7, -32.6), (4.8, -26.8)),
            ((3.8, -21.6), (1.7, -17.2), (0.0, -13.2)),
            ((-1.7, -17.2), (-3.8, -21.6), (-4.8, -26.8)),
            ((-5.7, -32.6), (-4.1, -38.0), (0.0, -39.0)),
        ),
    )
    sp_refined_large = _closed_cubic_path(
        (0.0, -10.8),
        (
            ((11.6, -9.0), (14.8, 5.8), (12.6, 21.5)),
            ((10.6, 33.0), (5.2, 39.8), (0.0, 44.0)),
            ((-5.2, 39.8), (-10.6, 33.0), (-12.6, 21.5)),
            ((-14.8, 5.8), (-11.6, -9.0), (0.0, -10.8)),
        ),
    )

    sigma_large = _closed_cubic_path(
        (0.0, -38.0),
        (
            ((12.0, -35.0), (16.0, -20.0), (14.0, -8.0)),
            ((12.0, 3.0), (7.0, 15.0), (0.0, 28.0)),
            ((-7.0, 15.0), (-12.0, 3.0), (-14.0, -8.0)),
            ((-16.0, -20.0), (-12.0, -35.0), (0.0, -38.0)),
        ),
    )
    sigma_small = _closed_cubic_path(
        (0.0, 22.0),
        (
            ((4.5, 21.0), (6.5, 24.0), (6.0, 28.0)),
            ((5.0, 31.0), (2.5, 33.0), (0.0, 34.0)),
            ((-2.5, 33.0), (-5.0, 31.0), (-6.0, 28.0)),
            ((-6.5, 24.0), (-4.5, 21.0), (0.0, 22.0)),
        ),
    )

    dz2_top = _closed_cubic_path(
        (0.0, -37.0),
        (
            ((9.0, -35.0), (13.0, -23.0), (11.0, -14.0)),
            ((9.0, -6.0), (4.0, -1.0), (0.0, 0.0)),
            ((-4.0, -1.0), (-9.0, -6.0), (-11.0, -14.0)),
            ((-13.0, -23.0), (-9.0, -35.0), (0.0, -37.0)),
        ),
    )
    dz2_bottom = _transform_path(dz2_top, scale_x=1.0, scale_y=-1.0)
    dz2_ring = _ring_path(-33.0, -7.0, 66.0, 14.0, 0.44, 0.36)
    dz2_refined_top = _closed_cubic_path(
        (0.0, -39.0),
        (
            ((8.2, -37.4), (12.2, -26.0), (10.4, -17.0)),
            ((8.8, -10.0), (4.2, -6.2), (0.0, -5.6)),
            ((-4.2, -6.2), (-8.8, -10.0), (-10.4, -17.0)),
            ((-12.2, -26.0), (-8.2, -37.4), (0.0, -39.0)),
        ),
    )
    dz2_refined_bottom = _transform_path(dz2_refined_top, scale_x=1.0, scale_y=-1.0)
    dz2_refined_ring = _ring_path(-31.0, -4.8, 62.0, 9.6, 0.46, 0.26)

    pi_top = _transform_path(dz2_top, translate=(0.0, -2.0))
    pi_bottom = _transform_path(pi_top, scale_x=1.0, scale_y=-1.0)
    pi_ring = _ring_path(-24.0, -5.5, 48.0, 11.0, 0.44, 0.40)

    f_top = _transform_path(_vertical_lobe_small(), translate=(0.0, -20.0))
    f_bottom = _transform_path(f_top, scale_x=1.0, scale_y=-1.0)
    f_left = _transform_path(_horizontal_petal(), translate=(-4.0, 0.0), scale_x=0.64, scale_y=0.64)
    f_right = _transform_path(f_left, scale_x=-1.0, scale_y=1.0)
    f_diag_a = _diag_petal_upper_left()
    f_diag_b = _diag_petal_lower_right()

    fz3_top = _transform_path(_vertical_lobe_small(), translate=(0.0, -28.0))
    fz3_mid = _transform_path(_vertical_lobe_micro(), translate=(0.0, 0.0))
    fz3_bottom = _transform_path(_vertical_lobe_small(), scale_x=1.0, scale_y=-1.0, translate=(0.0, 28.0))
    fz3_ring = _ring_path(-23.0, -6.0, 46.0, 12.0, 0.48, 0.46)

    library = {
        "s_round": GlyphDefinition(
            id="s_round",
            paths_outline=(GlyphLayer(sphere, "outline"),),
            paths_shaded=(GlyphLayer(sphere, "shaded", gradient="radial"),),
            paths_solid=(GlyphLayer(sphere, "solid", stroke=True),),
            default_bounds=circle_bounds,
            anchor_center=QPointF(0.0, 0.0),
            visual_padding=6.0,
        ),
        "torus_flat": GlyphDefinition(
            id="torus_flat",
            paths_outline=(GlyphLayer(torus, "outline"),),
            paths_shaded=(GlyphLayer(torus, "shaded", gradient="radial"),),
            paths_solid=(GlyphLayer(torus, "solid", stroke=False),),
            default_bounds=QRectF(-36.0, -10.0, 72.0, 20.0),
            anchor_center=QPointF(0.0, 0.0),
            visual_padding=6.0,
        ),
        "p_vertical": GlyphDefinition(
            id="p_vertical",
            paths_outline=(
                GlyphLayer(p_top, "outline"),
                GlyphLayer(p_bottom, "outline"),
            ),
            paths_shaded=(
                GlyphLayer(p_refined_top, "shaded", stroke=False),
                GlyphLayer(p_refined_bottom, "outline"),
            ),
            paths_solid=(
                GlyphLayer(p_refined_top, "solid", stroke=False),
                GlyphLayer(p_refined_bottom, "outline"),
            ),
            default_bounds=QRectF(-16.0, -41.0, 32.0, 82.0),
            anchor_center=QPointF(0.0, 0.0),
            visual_padding=4.5,
        ),
        "d_clover": GlyphDefinition(
            id="d_clover",
            paths_outline=(
                GlyphLayer(d_left, "outline"),
                GlyphLayer(d_top, "outline"),
                GlyphLayer(d_right, "outline"),
                GlyphLayer(d_bottom, "outline"),
            ),
            paths_shaded=(
                GlyphLayer(d_refined_left, "outline"),
                GlyphLayer(d_refined_top, "shaded", stroke=False),
                GlyphLayer(d_refined_right, "outline"),
                GlyphLayer(d_refined_bottom, "shaded", stroke=False),
            ),
            paths_solid=(
                GlyphLayer(d_refined_left, "outline"),
                GlyphLayer(d_refined_top, "solid", stroke=False),
                GlyphLayer(d_refined_right, "outline"),
                GlyphLayer(d_refined_bottom, "solid", stroke=False),
            ),
            default_bounds=QRectF(-36.0, -38.0, 72.0, 76.0),
            anchor_center=QPointF(0.0, 0.0),
            visual_padding=5.0,
        ),
        "sp3_hybrid": GlyphDefinition(
            id="sp3_hybrid",
            paths_outline=(
                GlyphLayer(sp_top, "outline"),
                GlyphLayer(sp_bottom, "outline"),
            ),
            paths_shaded=(
                GlyphLayer(sp_top, "shaded"),
                GlyphLayer(sp_bottom, "outline"),
            ),
            paths_solid=(
                GlyphLayer(sp_top, "solid", stroke=False),
                GlyphLayer(sp_bottom, "outline"),
            ),
            default_bounds=QRectF(-21.0, -50.0, 42.0, 101.0),
            anchor_center=QPointF(0.0, 0.0),
            visual_padding=6.0,
        ),
        "sp_lobe": GlyphDefinition(
            id="sp_lobe",
            paths_outline=(GlyphLayer(sp_single, "outline"),),
            paths_shaded=(
                GlyphLayer(sp_refined_small, "shaded", stroke=False),
                GlyphLayer(sp_refined_large, "shaded", stroke=False),
            ),
            paths_solid=(
                GlyphLayer(sp_refined_small, "solid", stroke=False),
                GlyphLayer(sp_refined_large, "solid", stroke=False),
            ),
            default_bounds=QRectF(-17.0, -38.0, 34.0, 76.0),
            anchor_center=QPointF(0.0, 0.0),
            visual_padding=5.0,
        ),
        "sigma_bonding": GlyphDefinition(
            id="sigma_bonding",
            paths_outline=(
                GlyphLayer(sigma_large, "outline"),
                GlyphLayer(sigma_small, "outline"),
            ),
            paths_shaded=(
                GlyphLayer(sigma_large, "outline"),
                GlyphLayer(sigma_small, "shaded"),
            ),
            paths_solid=(
                GlyphLayer(sigma_large, "outline"),
                GlyphLayer(sigma_small, "solid", stroke=False),
            ),
            default_bounds=QRectF(-16.0, -38.0, 32.0, 72.0),
            anchor_center=QPointF(0.0, 0.0),
            visual_padding=5.0,
        ),
        "dz2_axial": GlyphDefinition(
            id="dz2_axial",
            paths_outline=(
                GlyphLayer(dz2_top, "outline"),
                GlyphLayer(dz2_ring, "outline"),
                GlyphLayer(dz2_bottom, "outline"),
            ),
            paths_shaded=(
                GlyphLayer(dz2_refined_top, "shaded", stroke=False),
                GlyphLayer(dz2_refined_ring, "shaded", gradient="radial", stroke=False),
                GlyphLayer(dz2_refined_bottom, "shaded", stroke=False),
            ),
            paths_solid=(
                GlyphLayer(dz2_top, "solid", stroke=False),
                GlyphLayer(dz2_ring, "solid", stroke=False),
                GlyphLayer(dz2_bottom, "solid", stroke=False),
            ),
            default_bounds=QRectF(-31.0, -39.0, 62.0, 78.0),
            anchor_center=QPointF(0.0, 0.0),
            visual_padding=5.0,
        ),
        "pi_bonding": GlyphDefinition(
            id="pi_bonding",
            paths_outline=(
                GlyphLayer(pi_top, "outline"),
                GlyphLayer(pi_ring, "outline"),
                GlyphLayer(pi_bottom, "outline"),
            ),
            paths_shaded=(
                GlyphLayer(pi_top, "outline"),
                GlyphLayer(pi_ring, "shaded", gradient="radial"),
                GlyphLayer(pi_bottom, "outline"),
            ),
            paths_solid=(
                GlyphLayer(pi_top, "outline"),
                GlyphLayer(pi_ring, "solid", stroke=False),
                GlyphLayer(pi_bottom, "outline"),
            ),
            default_bounds=QRectF(-24.0, -39.0, 48.0, 78.0),
            anchor_center=QPointF(0.0, 0.0),
            visual_padding=5.0,
        ),
        "f_clover": GlyphDefinition(
            id="f_clover",
            paths_outline=(
                GlyphLayer(f_left, "outline"),
                GlyphLayer(f_top, "outline"),
                GlyphLayer(f_diag_a, "outline"),
                GlyphLayer(f_right, "outline"),
                GlyphLayer(f_bottom, "outline"),
                GlyphLayer(f_diag_b, "outline"),
            ),
            paths_shaded=(
                GlyphLayer(f_left, "outline"),
                GlyphLayer(f_top, "shaded"),
                GlyphLayer(f_diag_a, "shaded"),
                GlyphLayer(f_right, "outline"),
                GlyphLayer(f_bottom, "shaded"),
                GlyphLayer(f_diag_b, "outline"),
            ),
            paths_solid=(
                GlyphLayer(f_left, "outline"),
                GlyphLayer(f_top, "solid", stroke=False),
                GlyphLayer(f_diag_a, "solid", stroke=False),
                GlyphLayer(f_right, "outline"),
                GlyphLayer(f_bottom, "solid", stroke=False),
                GlyphLayer(f_diag_b, "outline"),
            ),
            default_bounds=QRectF(-42.0, -42.0, 84.0, 84.0),
            anchor_center=QPointF(0.0, 0.0),
            visual_padding=10.0,
        ),
        "fz3_axial": GlyphDefinition(
            id="fz3_axial",
            paths_outline=(
                GlyphLayer(fz3_top, "outline"),
                GlyphLayer(fz3_mid, "outline"),
                GlyphLayer(fz3_ring, "outline"),
                GlyphLayer(fz3_bottom, "outline"),
            ),
            paths_shaded=(
                GlyphLayer(fz3_top, "outline"),
                GlyphLayer(fz3_mid, "shaded"),
                GlyphLayer(fz3_ring, "shaded", gradient="radial"),
                GlyphLayer(fz3_bottom, "outline"),
            ),
            paths_solid=(
                GlyphLayer(fz3_top, "outline"),
                GlyphLayer(fz3_mid, "solid", stroke=False),
                GlyphLayer(fz3_ring, "solid", stroke=False),
                GlyphLayer(fz3_bottom, "outline"),
            ),
            default_bounds=QRectF(-23.0, -50.0, 46.0, 100.0),
            anchor_center=QPointF(0.0, 0.0),
            visual_padding=10.0,
        ),
    }
    return library


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
    def _gradient_for_path(path: QPainterPath, mode: str) -> QBrush:
        rect = path.boundingRect()
        if mode == "radial":
            radius = max(rect.width(), rect.height()) * 0.88
            focus = QPointF(rect.left() + rect.width() * 0.34, rect.top() + rect.height() * 0.28)
            gradient = QRadialGradient(focus, radius, focus)
            gradient.setColorAt(0.0, QColor("#FFFFFF"))
            gradient.setColorAt(0.38, QColor("#D9D9D9"))
            gradient.setColorAt(0.72, QColor("#676767"))
            gradient.setColorAt(1.0, QColor("#101010"))
            return QBrush(gradient)

        gradient = QLinearGradient(rect.left(), rect.top(), rect.right(), rect.bottom())
        gradient.setColorAt(0.0, QColor("#FFFFFF"))
        gradient.setColorAt(0.40, QColor("#D8D8D8"))
        gradient.setColorAt(0.72, QColor("#727272"))
        gradient.setColorAt(1.0, QColor("#151515"))
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
            painter.setBrush(self._gradient_for_path(layer.path, layer.gradient))
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
        delta = inner.center() - bbox.center()
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
