"""Generadores de SVG dinámico para iconos parametrizados (sin Qt).

Fase 2 del plan de modernización de la UI (ver OpenSpec
``2026-09-24-modernize-ui-svg-icons``). Estos builders producen cadenas SVG
24×24 en el mismo lenguaje que el set estático
(``theme/icons/i-*.svg``): trazo ``currentColor`` 1.75, esquinas
redondeadas. Los colores literales solo aparecen cuando el caller pasa un
color de dominio (CPK/semántico): fichas de átomo, esfera de coordinación y
relleno de diagramas de energía.

El :class:`~chemuson.gui.theme.icon_provider.IconProvider` resuelve
``icon_dynamic(key, ...)``, con ``key`` en :data:`BUILDERS`, a la función
correspondiente y pasa el resultado por su pipeline de tinte/caché/HiDPI.

Nota: funciones puras (str → str) para poder testearlas sin QApplication.
"""

from __future__ import annotations

import html
import math
from typing import Callable

__all__ = [
    "BUILDERS",
    "glyph_svg",
    "atom_svg",
    "sphere_svg",
    "charge_svg",
    "electrons_svg",
    "radical_svg",
    "ring_svg",
    "ring_template_svg",
    "energy_boxes_svg",
]

_VIEWBOX = "0 0 24 24"
_STROKE_W = "1.75"


def _wrap(body: str) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" '
        f'viewBox="{_VIEWBOX}">'
        f'<g fill="none" stroke="currentColor" stroke-width="{_STROKE_W}" '
        f'stroke-linecap="round" stroke-linejoin="round">{body}</g></svg>'
    )


def _esc(text: str) -> str:
    return html.escape(str(text), quote=True)


# ---------------------------------------------------------------------------
# Utilidades de color (hex puro, sin Qt)
# ---------------------------------------------------------------------------

def _parse_hex(value: str) -> tuple[int, int, int] | None:
    """Parsea ``#RGB``/``#RRGGBB`` a ``(r, g, b)`` o ``None`` si no es válido."""
    s = str(value).strip().lstrip("#")
    if len(s) == 3:
        s = "".join(ch * 2 for ch in s)
    if len(s) != 6:
        return None
    try:
        return (int(s[0:2], 16), int(s[2:4], 16), int(s[4:6], 16))
    except ValueError:
        return None


def _rgb_hex(r: int, g: int, b: int) -> str:
    return "#{:02X}{:02X}{:02X}".format(
        max(0, min(255, int(round(r)))),
        max(0, min(255, int(round(g)))),
        max(0, min(255, int(round(b)))),
    )


def _scale_hex(color: str, factor: float, default: str = "#8D99A6") -> str:
    rgb = _parse_hex(color)
    if rgb is None:
        return default
    return _rgb_hex(rgb[0] * factor, rgb[1] * factor, rgb[2] * factor)


def _luminance(color: str) -> float:
    """Luminancia BT.601 (0–255); 0 si el color no es válido."""
    rgb = _parse_hex(color)
    if rgb is None:
        return 0.0
    return (rgb[0] * 299 + rgb[1] * 587 + rgb[2] * 114) / 1000


def _contrast_text(color: str) -> str:
    """Texto blanco sobre fondos oscuros y negro sobre claros (umbral 128)."""
    return "#FFFFFF" if _luminance(color) < 128 else "#000000"


# ---------------------------------------------------------------------------
# Builders (cada uno recibe params y devuelve el SVG completo)
# ---------------------------------------------------------------------------

def glyph_svg(
    label: str = "?",
    *,
    font_size: float = 12.5,
    weight: int = 700,
    shape: str = "none",
    **_ignored: object,
) -> str:
    """Glifo tipográfico (letras de átomo, corchetes, signos, etc.).

    ``shape="circle"`` dibuja un círculo de fondo (fichas de elemento).
    """
    bg = ""
    if shape == "circle":
        bg = (
            f'<circle cx="12" cy="12" r="9" fill="none" '
            f'stroke="currentColor" stroke-width="{_STROKE_W}"/>'
        )
    return _wrap(
        f"{bg}"
        f'<text x="12" y="16.6" text-anchor="middle" font-family="sans-serif" '
        f'font-size="{float(font_size)}" font-weight="{int(weight)}" '
        f'fill="currentColor" stroke="none">{_esc(label)}</text>'
    )


def atom_svg(text: str = "C", *, fill: str = "#333333", **_ignored: object) -> str:
    """Ficha de átomo: círculo CPK relleno + símbolo centrado.

    ``fill`` es un color de dominio (CPK), no del tema; el color del texto
    se decide por luminancia (umbral 128), como en la versión QPainter.
    """
    fill = fill if _parse_hex(fill) is not None else "#333333"
    n = max(1, len(str(text)))
    font_size = 12.0 if n == 1 else (9.5 if n == 2 else 8.0)
    stroke = _scale_hex(fill, 0.78)
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" '
        f'viewBox="{_VIEWBOX}">'
        f'<circle cx="12" cy="12" r="9.5" fill="{fill}" stroke="{stroke}" '
        f'stroke-width="1.25"/>'
        f'<text x="12" y="16.1" text-anchor="middle" font-family="sans-serif" '
        f'font-size="{font_size}" font-weight="700" '
        f'fill="{_contrast_text(fill)}" stroke="none">{_esc(text)}</text></svg>'
    )


def sphere_svg(fill: str = "#8D99A6", **_ignored: object) -> str:
    """Esfera de coordinación genérica con gradiente radial (color de dominio)."""
    base = fill if _parse_hex(fill) is not None else "#8D99A6"
    highlight = _scale_hex(base, 1.70, default=base)
    shadow = _scale_hex(base, 1.0 / 1.65, default=base)
    border = _scale_hex(base, 0.55, default=base)
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" '
        f'viewBox="{_VIEWBOX}">'
        f'<defs><radialGradient id="sph" gradientUnits="userSpaceOnUse" '
        f'cx="9.3" cy="9.3" r="16.4" fx="8.1" fy="8.1">'
        f'<stop offset="0" stop-color="{highlight}"/>'
        f'<stop offset="0.55" stop-color="{base}"/>'
        f'<stop offset="1" stop-color="{shadow}"/>'
        f"</radialGradient></defs>"
        f'<circle cx="12" cy="12" r="8.2" fill="url(#sph)" stroke="{border}" '
        f'stroke-width="1.25"/></svg>'
    )


def charge_svg(sign: str = "+", **_ignored: object) -> str:
    """Carga circular con signo ``+`` o ``-``."""
    plus = '<line x1="12" y1="7.75" x2="12" y2="16.25"/>' if str(sign) == "+" else ""
    return _wrap(
        f'<circle cx="12" cy="12" r="8"/>'
        f'<line x1="7.75" y1="12" x2="16.25" y2="12"/>{plus}'
    )


def electrons_svg(
    count: int = 1,
    *,
    spread: float = 6.0,
    **_ignored: object,
) -> str:
    """Puntos de electrones (desapareado=1, par solitario=2, ...)."""
    n = max(1, int(count))
    # Escala la separación de la rejilla 32 (valor histórico) a la de 24 px.
    gap = max(2.0, min(14.0, float(spread) * 24.0 / 32.0))
    radius = 1.7
    dots = []
    for i in range(n):
        x = 12.0 + (i - (n - 1) / 2.0) * gap
        dots.append(
            f'<circle cx="{x:.2f}" cy="12" r="{radius}" '
            f'fill="currentColor" stroke="none"/>'
        )
    return _wrap("".join(dots))


def radical_svg(sign: str = "+", **_ignored: object) -> str:
    """Radical (punto) con signo de carga opcional en la esquina superior."""
    plus = '<line x1="15.4" y1="5.9" x2="15.4" y2="12.9"/>' if str(sign) == "+" else ""
    return _wrap(
        f'<circle cx="9.8" cy="14.2" r="1.8" fill="currentColor" stroke="none"/>'
        f'<line x1="11.9" y1="9.4" x2="18.9" y2="9.4"/>{plus}'
    )


def _ring_points(sides: int, radius: float) -> list[tuple[float, float]]:
    """Vértices del polígono con vértice arriba (orientación ChemDraw)."""
    n = max(3, int(sides))
    points: list[tuple[float, float]] = []
    for k in range(n):
        angle = math.radians(90.0 - k * 360.0 / n)
        points.append((12.0 + radius * math.cos(angle), 12.0 - radius * math.sin(angle)))
    return points


def ring_svg(
    sides: int = 6,
    *,
    aromatic: bool = True,
    **_ignored: object,
) -> str:
    """Anillo de ``sides`` lados (3+) con círculo interno si es aromático."""
    points = _ring_points(sides, 8.0)
    pts = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
    inner = '<circle cx="12" cy="12" r="4.25"/>' if aromatic else ""
    return _wrap(f'<polygon points="{pts}"/>{inner}')


def ring_template_svg(
    label: str = "",
    *,
    sides: int = 6,
    **_ignored: object,
) -> str:
    """Anillo con borde inferior en negrita (sugerencia Haworth) y etiqueta."""
    points = _ring_points(sides, 7.6)
    pts = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
    lowest = sorted(points, key=lambda p: p[1], reverse=True)[:2]
    bold_edge = ""
    if len(lowest) == 2:
        (ax, ay), (bx, by) = lowest
        bold_edge = (
            f'<line x1="{ax:.2f}" y1="{ay:.2f}" x2="{bx:.2f}" y2="{by:.2f}" '
            f'stroke-width="3.2"/>'
        )
    label_svg = ""
    if label:
        font_size = 8.5 if len(str(label)) <= 2 else 7.5
        label_svg = (
            f'<text x="12" y="15.6" text-anchor="middle" font-family="sans-serif" '
            f'font-size="{font_size}" font-weight="700" '
            f'fill="currentColor" stroke="none">{_esc(label)}</text>'
        )
    return _wrap(f'<polygon points="{pts}"/>{bold_edge}{label_svg}')


def energy_boxes_svg(
    boxes: int = 1,
    *,
    label: str = "",
    side: str = "left",
    fill: str = "#FFFFFF",
    stroke: bool = True,
    **_ignored: object,
) -> str:
    """Cajas de configuración electrónica (1..N) con flecha en la caja central.

    ``fill`` es un color de dominio que pasa el caller (preset del diagrama);
    el trazo usa ``currentColor``. Para N grande las cajas se fusionan en una
    banda (equivalente visual al pixmap subpíxel de la versión QPainter).
    """
    n = max(1, min(64, int(boxes)))
    margin_x, margin_y, height = 3.0, 5.5, 13.0
    bottom = margin_y + height
    label = str(label or "").strip()
    side = side if side in ("left", "right") else "left"
    label_w = 6.5 if (label and side in ("left", "right")) else 0.0
    gap_label = 1.5 if label_w else 0.0
    x0 = margin_x + (label_w + gap_label if side == "left" else 0.0)
    x1 = 24.0 - margin_x - (label_w + gap_label if side == "right" else 0.0)
    width = max(2.0, x1 - x0)

    body: list[str] = []
    if label:
        cx = (margin_x + (margin_x + label_w)) / 2.0 if side == "left" else (
            (x1 - label_w + x1) / 2.0
        )
        body.append(
            f'<text x="{cx:.2f}" y="12.3" text-anchor="middle" '
            f'font-family="sans-serif" font-size="7" font-weight="700" '
            f'fill="currentColor" stroke="none">{_esc(label)}</text>'
        )

    fill_attr = f' fill="{fill}"' if fill else ' fill="none"'
    stroke_attr = f' stroke="currentColor"' if stroke else ' stroke="none"'

    def box_attrs() -> str:
        return f"{fill_attr}{stroke_attr} stroke-width='1'"

    gap = 1.1 if n <= 14 else 0.6
    box_w = (width - gap * (n - 1)) / n if n > 1 else width
    if box_w < 0.7:
        # Banda fusionada (N muy grande): misma lectura que antes.
        body.append(f"<rect x='{x0:.2f}' y='{margin_y}' width='{width:.2f}' height='{height}'{box_attrs()} />")
    else:
        for i in range(n):
            bx = x0 + i * (box_w + gap)
            body.append(
                f"<rect x='{bx:.2f}' y='{margin_y}' width='{box_w:.2f}' "
                f"height='{height}'{box_attrs()} />"
            )
        if box_w >= 2.4:
            idx = min(n - 1, n // 2)
            cx = x0 + idx * (box_w + gap) + box_w / 2.0
            body.append(
                f"<line x1='{cx:.2f}' y1='{bottom - 1.5}' x2='{cx:.2f}' y2='8.5' "
                f"stroke='currentColor' stroke-width='1.3' fill='none'/>"
                f"<path d='M{cx - 1.7:.2f} 11.1L{cx:.2f} 8.5L{cx + 1.7:.2f} 11.1' "
                f"stroke='currentColor' stroke-width='1.3' fill='none'/>"
            )
    return _wrap("".join(body))


#: Registro resuelto por ``IconProvider.icon_dynamic`` (key → builder).
BUILDERS: dict[str, Callable[..., str]] = {
    "glyph": glyph_svg,
    "atom": atom_svg,
    "sphere": sphere_svg,
    "charge": charge_svg,
    "electrons": electrons_svg,
    "radical": radical_svg,
    "ring": ring_svg,
    "ring-template": ring_template_svg,
    "energy-boxes": energy_boxes_svg,
}
