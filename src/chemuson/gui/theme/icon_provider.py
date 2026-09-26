"""IconProvider: SVG → QIcon/QPixmap con tinte, caché y HiDPI.

Fase 1 (infraestructura) del plan de modernización de la UI; la Fase 2
(``openspec/changes/2026-09-24-modernize-ui-svg-icons``) puebla el set de
SVG estáticos, añade la API dinámica (``icon_dynamic``/``pixmap_dynamic``)
y convierte ``gui/icons.py`` en fachada sobre este provider. Modelo del
provider del spike aprobado (``docs/ui-modernization/pyqt6-spike/icons.py``).

Contrato:
- Iconos estáticos: ``src/chemuson/gui/theme/icons/i-<name>.svg`` (rejilla
  24 px, trazo ``currentColor``).
- Iconos dinámicos: ``icon_dynamic(key, color, size, **params)`` resuelve
  el SVG con el builder ``icon_svg.BUILDERS[key]`` (anillos parametrizados,
  átomos CPK, cargas, electrones, diagramas de energía...), pasa por el
  mismo render/tinte/HiDPI y se cachea con clave que incluye el key y los
  params canónicos.
- Tinte: sustitución de ``currentColor`` en el SVG + ``QSvgRenderer``
  (sin iterar píxeles). El color se pasa explícitamente y es parte de la
  clave de caché: el provider es theme-aware (el cambio light/dark produce
  iconos correctos por tema sin contaminación entre cachés).

Notas de uso (pitfalls del spike):
- ``icon()`` devuelve ``QIcon`` (para ``QToolButton``); ``QLabel`` necesita
  ``pixmap()`` (pasar un ``QIcon`` donde toca un ``QPixmap`` provoca
  ``qFatal`` en QtWidgets).
- Fallo visible: nombre ausente → ``QIcon()``/``QPixmap()`` vacíos, sin
  excepciones (igual para ``key`` desconocido en la API dinámica).
"""

from __future__ import annotations

from pathlib import Path

from PyQt6.QtCore import QByteArray, QRectF, Qt
from PyQt6.QtGui import QGuiApplication, QIcon, QPainter, QPixmap
from PyQt6.QtSvg import QSvgRenderer

from chemuson.gui.theme import icon_svg

__all__ = [
    "DEFAULT_ICONS_DIR",
    "IconProvider",
]

#: Carpeta de SVGs estáticos (poblada por la Fase 2).
DEFAULT_ICONS_DIR: Path = Path(__file__).resolve().parent / "icons"

_VIEWBOX = "0 0 24 24"


def _render_svg(data: bytes, size: int, dpr: float) -> QPixmap:
    """Renderiza un SVG a ``QPixmap`` con ratio de dispositivo ``dpr``."""
    renderer = QSvgRenderer(QByteArray(data))
    px = max(1, int(round(size * dpr)))
    pixmap = QPixmap(px, px)
    pixmap.fill(Qt.GlobalColor.transparent)
    if renderer.isValid():
        # Mantener DPR=1 durante el pintado y dibujar explícitamente en
        # coordenadas físicas evita el escalado/clipping del backend HiDPI.
        painter = QPainter(pixmap)
        renderer.render(painter, QRectF(0.0, 0.0, float(px), float(px)))
        painter.end()
    # Etiquetar el backing sólo después del render; su tamaño lógico queda
    # como px / dpr (21×21 a DPR 2 -> 42×42 físicos).
    pixmap.setDevicePixelRatio(dpr)
    return pixmap


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
        self._dpr = float(dpr) if dpr is not None else self.application_device_pixel_ratio()
        self._icons_dir = Path(icons_dir) if icons_dir is not None else DEFAULT_ICONS_DIR
        self._icon_cache: dict[tuple[str, str, int, float], QIcon] = {}
        self._pixmap_cache: dict[tuple[str, str, int, float], QPixmap] = {}
        self._svg_cache: dict[str, bytes] = {}

    @staticmethod
    def application_device_pixel_ratio() -> float:
        app = QGuiApplication.instance()
        if app is None:
            return 1.0
        screen = app.primaryScreen()
        return float(screen.devicePixelRatio()) if screen is not None else 1.0

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
            # Formato: glyph:<label>|<shape>|<font_size>
            _, rest = name.split(":", 1)
            label, shape, fsize = (rest + "||12.5|none").split("|", 3)[:3]
            data = icon_svg.glyph_svg(
                label, font_size=float(fsize), shape=shape
            ).encode()
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

        Caché por ``(name, color, size, dpr)``: un cambio de DPR no reutiliza
        un backing store renderizado para otra pantalla.
        """
        key = (name, color, size, self._dpr)
        cached = self._icon_cache.get(key)
        if cached is not None:
            return cached
        icon = self._to_icon(name, color, size)
        self._icon_cache[key] = icon
        return icon

    def set_device_pixel_ratio(self, dpr: float) -> None:
        """Actualiza el DPR objetivo e invalida sólo representaciones Qt.

        SVG ya cargados permanecen en caché; QIcon/QPixmap se vuelven a
        rasterizar para el nuevo monitor/DPR.
        """
        resolved = float(dpr)
        if resolved <= 0:
            raise ValueError("dpr must be positive")
        if resolved == self._dpr:
            return
        self._dpr = resolved
        self._icon_cache.clear()
        self._pixmap_cache.clear()

    def pixmap(self, name: str, color: str, size: int = 20) -> QPixmap:
        """``QPixmap`` suelto teñido (para ``QLabel``/``paintEvent``).

        Vacío (``isNull()``) si el icono no existe.
        """
        key = (name, color, size, self._dpr)
        cached = self._pixmap_cache.get(key)
        if cached is not None:
            return cached
        data = self._svg(name)
        pixmap = QPixmap() if data is None else self._tinted_pixmap(data, color, size)
        self._pixmap_cache[key] = pixmap
        return pixmap

    # ------------------------------------------------------------------
    # API dinámica (Fase 2): iconos parametrizados vía icon_svg.BUILDERS
    # ------------------------------------------------------------------
    @staticmethod
    def _dynamic_name(key: str, params: dict[str, object]) -> str | None:
        """Clave canónica de SVG para un icono dinámico (o ``None`` si el
        ``key`` no tiene builder registrado: fallo visible, sin excepción).
        """
        builder = icon_svg.BUILDERS.get(key)
        if builder is None:
            return None
        canonical = "|".join(
            f"{name}={params[name]!r}" for name in sorted(params)
        )
        return f"dyn:{key}:{canonical}"

    def icon_dynamic(
        self, key: str, color: str, size: int = 20, **params: object
    ) -> QIcon:
        """``QIcon`` dinámico teñido para ``key`` (p. ej. ``"ring"``,
        ``"atom"``, ``"energy-boxes"``) o vacío si no hay builder.

        La clave de caché incluye ``key``, los ``params`` canónicos, el
        ``color`` y el ``size``: cambiar parámetros o tema produce iconos
        distintos y el uso repetido devuelve la misma instancia.
        """
        name = self._dynamic_name(key, params)
        if name is None:
            return QIcon()
        cache_key = (name, color, size, self._dpr)
        cached = self._icon_cache.get(cache_key)
        if cached is not None:
            return cached
        icon = QIcon(self._dynamic_pixmap(name, color, size, key=key, params=params))
        self._icon_cache[cache_key] = icon
        return icon

    def pixmap_dynamic(
        self, key: str, color: str, size: int = 20, **params: object
    ) -> QPixmap:
        """``QPixmap`` dinámico teñido; vacío si ``key`` no existe."""
        name = self._dynamic_name(key, params)
        if name is None:
            return QPixmap()
        cache_key = (name, color, size, self._dpr)
        cached = self._pixmap_cache.get(cache_key)
        if cached is not None:
            return cached
        pixmap = self._dynamic_pixmap(name, color, size, key=key, params=params)
        self._pixmap_cache[cache_key] = pixmap
        return pixmap

    def _dynamic_pixmap(
        self, name: str, color: str, size: int, *, key: str, params: dict[str, object]
    ) -> QPixmap:
        builder = icon_svg.BUILDERS[key]
        data = builder(**params).encode()
        self._svg_cache[name] = data
        return self._tinted_pixmap(data, color, size)

    # ------------------------------------------------------------------
    # Internos
    # ------------------------------------------------------------------
    def _tinted_pixmap(self, data: bytes, color: str, size: int) -> QPixmap:
        tinted = data.replace(b"currentColor", color.encode())
        return _render_svg(tinted, size, self._dpr)

    def _to_icon(self, name: str, color: str, size: int) -> QIcon:
        if name.startswith("dyn:"):
            # Los dinámicos se resuelven en ``icon_dynamic`` (donde hay
            # ``key``/``params``); aquí solo por si se llama directamente.
            return QIcon()
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
