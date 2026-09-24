"""IconProvider: SVG → QIcon/QPixmap con tinte, caché y HiDPI.

Fase 1 (infraestructura) del plan de modernización de la UI; la migración
de los iconos de ``gui/icons.py`` a este provider es la Fase 2 (contrato
documentado abajo). Modelo del provider del spike aprobado
(``docs/ui-modernization/pyqt6-spike/icons.py``), adaptado al paquete.

Contrato (para la Fase 2):
- Iconos estáticos: ``src/chemuson/gui/theme/icons/i-<name>.svg`` (rejilla
  24 px, trazo ``currentColor``); la carpeta se puebla en la Fase 2.
- Tinte: sustitución de ``currentColor`` en el SVG + ``QSvgRenderer``
  (sin iterar píxeles). El color se pasa explícitamente: la Fase 2 lo
  alimentará con tokens (``icon``/``iconHover``/``iconActive``) del tema
  activo, y al ser parte de la clave de caché el provider es theme-aware.
- ``gui/icons.py`` NO se modifica en esta fase; sus funciones seguirán
  pintando con ``QPainter`` hasta la delegación de la Fase 2.

Notas de uso (pitfalls del spike):
- ``icon()`` devuelve ``QIcon`` (para ``QToolButton``); ``QLabel`` necesita
  ``pixmap()`` (pasar un ``QIcon`` donde toca un ``QPixmap`` provoca
  ``qFatal`` en QtWidgets).
- Fallo visible: nombre ausente → ``QIcon()``/``QPixmap()`` vacíos, sin
  excepciones.
"""

from __future__ import annotations

import html
from pathlib import Path

from PyQt6.QtCore import QByteArray, Qt
from PyQt6.QtGui import QGuiApplication, QIcon, QPainter, QPixmap
from PyQt6.QtSvg import QSvgRenderer

__all__ = [
    "DEFAULT_ICONS_DIR",
    "IconProvider",
]

#: Carpeta de SVGs estáticos (la Fase 2 la puebla con el set curado).
DEFAULT_ICONS_DIR: Path = Path(__file__).resolve().parent / "icons"

_VIEWBOX = "0 0 24 24"


def _render_svg(data: bytes, size: int, dpr: float) -> QPixmap:
    """Renderiza un SVG a ``QPixmap`` con ratio de dispositivo ``dpr``."""
    renderer = QSvgRenderer(QByteArray(data))
    px = max(1, int(round(size * dpr)))
    pixmap = QPixmap(px, px)
    pixmap.setDevicePixelRatio(dpr)
    pixmap.fill(Qt.GlobalColor.transparent)
    if renderer.isValid():
        painter = QPainter(pixmap)
        renderer.render(painter)
        painter.end()
    return pixmap


def _glyph_svg(
    label: str,
    *,
    font_size: float = 12.5,
    weight: int = 700,
    shape: str = "none",
) -> bytes:
    """SVG mínimo para glifos de texto (letras de átomo, cargas, etc.).

    ``shape`` permite dibujar la forma de fondo (``"circle"``) cuando el
    glifo va sobre una forma (p. ej. fichas de elemento).
    """
    escaped = html.escape(label, quote=True)
    bg = ""
    if shape == "circle":
        bg = (
            '<circle cx="12" cy="12" r="9" fill="none" '
            'stroke="currentColor" stroke-width="1.75"/>'
        )
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" '
        f'viewBox="{_VIEWBOX}">'
        f"{bg}"
        f'<text x="12" y="16.6" text-anchor="middle" font-family="sans-serif" '
        f'font-size="{font_size}" font-weight="{weight}" '
        f'fill="currentColor" stroke="none">{escaped}</text></svg>'
    ).encode()


class IconProvider:
    """Carga SVG, tiñe y cachea ``QIcon``/``QPixmap``.

    Args:
        dpr: Ratio de dispositivo para HiDPI. ``None`` usa el ratio de la
            :class:`QGuiApplication` activa (o 1.0 sin aplicación).
        icons_dir: Carpeta de SVGs estáticos (por defecto
            ``gui/theme/icons/``).
    """

    def __init__(
        self,
        dpr: float | None = None,
        icons_dir: Path | str | None = None,
    ) -> None:
        if dpr is None:
            app = QGuiApplication.instance()
            dpr = app.devicePixelRatio() if app is not None else 1.0
        self._dpr = float(dpr)
        self._icons_dir = Path(icons_dir) if icons_dir is not None else DEFAULT_ICONS_DIR
        self._icon_cache: dict[tuple[str, str, int], QIcon] = {}
        self._pixmap_cache: dict[tuple[str, str, int], QPixmap] = {}
        self._svg_cache: dict[str, bytes] = {}

    # ------------------------------------------------------------------
    # Propiedades
    # ------------------------------------------------------------------
    @property
    def dpr(self) -> float:
        """Ratio de dispositivo con el que se renderizan los iconos."""
        return self._dpr

    @property
    def icons_dir(self) -> Path:
        """Carpeta de SVGs estáticos."""
        return self._icons_dir

    def clear_cache(self) -> None:
        """Vacía las cachés de SVG/iconos/pixmaps (p. ej. al cargar assets)."""
        self._icon_cache.clear()
        self._pixmap_cache.clear()
        self._svg_cache.clear()

    # ------------------------------------------------------------------
    # Carga de SVG
    # ------------------------------------------------------------------
    def _svg(self, name: str) -> bytes | None:
        """Devuelve los bytes SVG del nombre (estático o glifo) o ``None``."""
        if name in self._svg_cache:
            return self._svg_cache[name]
        data: bytes | None
        if name.startswith("glyph:"):
            # Formato: glyph:<label>|<shape>|<size>
            _, rest = name.split(":", 1)
            label, shape, fsize = (rest + "||12.5|none").split("|", 3)[:3]
            data = _glyph_svg(label, font_size=float(fsize), shape=shape)
        else:
            path = self._icons_dir / f"i-{name}.svg"
            data = path.read_bytes() if path.exists() else None
        if data is not None:
            self._svg_cache[name] = data
        return data

    # ------------------------------------------------------------------
    # API pública
    # ------------------------------------------------------------------
    def icon(self, name: str, color: str, size: int = 20) -> QIcon:
        """``QIcon`` teñido para ``name`` (o vacío si no existe).

        Caché por ``(name, color, size)``: la misma clave devuelve la misma
        instancia, de modo que cambiar de tema (cambiar ``color``) produce
        iconos distintos cacheados por tema.
        """
        key = (name, color, size)
        cached = self._icon_cache.get(key)
        if cached is not None:
            return cached
        icon = self._to_icon(name, color, size)
        self._icon_cache[key] = icon
        return icon

    def pixmap(self, name: str, color: str, size: int = 20) -> QPixmap:
        """``QPixmap`` suelto teñido (para ``QLabel``/``paintEvent``).

        Vacío (``isNull()``) si el icono no existe.
        """
        key = (name, color, size)
        cached = self._pixmap_cache.get(key)
        if cached is not None:
            return cached
        data = self._svg(name)
        pixmap = QPixmap() if data is None else self._tinted_pixmap(data, color, size)
        self._pixmap_cache[key] = pixmap
        return pixmap

    # ------------------------------------------------------------------
    # Internos
    # ------------------------------------------------------------------
    def _tinted_pixmap(self, data: bytes, color: str, size: int) -> QPixmap:
        tinted = data.replace(b"currentColor", color.encode())
        return _render_svg(tinted, size, self._dpr)

    def _to_icon(self, name: str, color: str, size: int) -> QIcon:
        data = self._svg(name)
        if data is None:
            return QIcon()  # vacío: fallo visible
        return QIcon(self._tinted_pixmap(data, color, size))

    # ------------------------------------------------------------------
    # Utilidades
    # ------------------------------------------------------------------
    @staticmethod
    def theme_color(theme_name: str, role: str = "icon") -> str:
        """Color de tinte para un estado de icono desde los tokens del tema.

        ``role``: ``"icon"``, ``"iconHover"`` o ``"iconActive"``.
        """
        from chemuson.gui.theme.tokens import theme_color

        return theme_color(theme_name, role).name()
