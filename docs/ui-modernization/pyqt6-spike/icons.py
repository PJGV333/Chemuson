"""IconProvider para el spike: SVG → QIcon teñido y con caché.

- Los iconos son archivos SVG (carpeta `icons/`) con trazo `currentColor`
  (rejilla 24 px, trazo 1.75, esquinas redondeadas — misma métrica que el
  mockup y que el set de la Fase 2 del PLAN.md).
- El tinte se hace sustituyendo `currentColor` en el SVG y renderizando con
  QSvgRenderer (sin iterar píxeles).
- Caché por (nombre, color, tamaño).
"""
from __future__ import annotations

import html
from pathlib import Path

from PyQt6.QtCore import QByteArray, Qt, QRectF
from PyQt6.QtGui import QColor, QIcon, QPainter, QPixmap, QBrush
from PyQt6.QtSvg import QSvgRenderer

ICONS_DIR = Path(__file__).resolve().parent / "icons"
VIEWBOX = "0 0 24 24"


def _render_svg(data: bytes, size: int, dpr: float) -> QPixmap:
    r = QSvgRenderer(QByteArray(data))
    px = max(1, int(size * dpr))
    pm = QPixmap(px, px)
    pm.setDevicePixelRatio(dpr)
    pm.fill(Qt.GlobalColor.transparent)
    if r.isValid():
        p = QPainter(pm)
        r.render(p)
        p.end()
    return pm


def _glyph_svg(label: str, *, font_size: float = 12.5, weight: int = 700,
               shape: str = "none") -> bytes:
    """SVG mínimo para glifos de texto (letras de átomo, cargas, etc.).

    `shape` permite dibujar la forma de fondo (círculo/rectángulo) cuando el
    glifo va sobre una forma (p. ej. fichas de elemento).
    """
    esc = html.escape(label, quote=True)
    bg = ""
    if shape == "circle":
        bg = '<circle cx="12" cy="12" r="9" fill="none" stroke="currentColor" stroke-width="1.75"/>'
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="{VIEWBOX}">'
        f"{bg}"
        f'<text x="12" y="16.6" text-anchor="middle" font-family="sans-serif" '
        f'font-size="{font_size}" font-weight="{weight}" '
        f'fill="currentColor" stroke="none">{esc}</text></svg>'
    ).encode()


class IconProvider:
    """Carga SVG, tiñe y cachea QIcons.

    Los iconos "estáticos" viven en `icons/<name>.svg` (name sin prefijo,
    p. ej. `pointer` → `icons/i-pointer.svg`). Los "glifo" se generan al vuelo
    (letras de átomo, cargas) y también se cachean.
    """

    def __init__(self, dpr: float):
        self._dpr = dpr
        self._cache: dict[tuple[str, str, int], QIcon] = {}
        self._svg_cache: dict[str, bytes] = {}

    # -- carga -------------------------------------------------------------
    def _svg(self, name: str) -> bytes | None:
        if name in self._svg_cache:
            return self._svg_cache[name]
        if name.startswith("glyph:"):
            # formato: glyph:<label>|<shape>|<size>
            _, rest = name.split(":", 1)
            label, shape, fsize = (rest + "||12.5|none").split("|", 3)[:3]
            data = _glyph_svg(label, font_size=float(fsize), shape=shape)
        else:
            p = ICONS_DIR / f"i-{name}.svg"
            if not p.exists():
                return None
            data = p.read_bytes()
        self._svg_cache[name] = data
        return data

    # -- API ---------------------------------------------------------------
    def icon(self, name: str, color: str, size: int = 21) -> QIcon:
        key = (name, color, size)
        ic = self._cache.get(key)
        if ic is not None:
            return ic
        data = self._svg(name)
        if data is None:
            ic = QIcon()  # vacío: fail visible
        else:
            tinted = data.replace(b"currentColor", color.encode())
            ic = QIcon(_render_svg(tinted, size, self._dpr))
        self._cache[key] = ic
        return ic

    def pixmap(self, name: str, color: str, size: int = 21) -> QPixmap:
        """Pixmap suelto (para pintarlo a mano en paintEvent)."""
        data = self._svg(name)
        if data is None:
            return QPixmap()
        tinted = data.replace(b"currentColor", color.encode())
        return _render_svg(tinted, size, self._dpr)

    def dot_pixmap(self, color: str, size: int = 8, radius: float | None = None) -> QPixmap:
        """Punto circular de tamaño arbitrario (indicadores, severidades)."""
        from PyQt6.QtGui import QPainter, QBrush
        px = max(1, int(size * self._dpr))
        pm = QPixmap(px, px)
        pm.setDevicePixelRatio(self._dpr)
        pm.fill(Qt.GlobalColor.transparent)
        p = QPainter(pm)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(QBrush(QColor(color)))
        p.drawEllipse(QRectF(0.5, 0.5, size - 1, size - 1))
        p.end()
        return pm
