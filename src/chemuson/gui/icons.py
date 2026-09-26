"""Biblioteca de iconos de Chemuson (fachada de compatibilidad sobre SVG).

Historial: hasta la Fase 2 de la modernización de la UI este módulo
generaba todos los iconos con ``QPainter`` (~1320 líneas). Ahora es una
**fachada 1:1** que conserva la API pública exacta (mismos nombres, firmas
y retorno ``QIcon``) pero delega en el ``IconProvider`` de
``chemuson.gui.theme``:

- Iconos estáticos: set SVG propio en ``chemuson/gui/theme/icons/``
  (24×24, trazo ``currentColor`` 1.75; ver ``icons/LICENSE.txt``).
- Iconos parametrizados (anillos, átomos CPK, cargas, electrones,
  diagramas de energía, esfera de coordinación): SVG dinámico generado por
  ``chemuson.gui.theme.icon_svg``.
- Tinte: el color de cada icono es parte de la clave de caché del provider,
  por lo que el cambio light→dark→light produce iconos correctos por tema
  sin pixel-loops ni ``QIcon.fromTheme`` (contrato documentado en el
  OpenSpec ``2026-09-24-modernize-ui-svg-icons``).

Los colores por defecto de los iconos monocromos provienen de los design
tokens de la Fase 1 (``icon``/``text3``/``surface3``/``surface``). Los
colores químicos/CPK (:data:`ATOM_COLORS`, relleno de átomos, esfera de
coordinación, ``fill_color`` de diagramas de energía) son colores de
dominio y no cambian con el tema.
"""

from __future__ import annotations

from PyQt6.QtGui import QIcon, QPixmap
from PyQt6.QtCore import Qt

from chemuson.gui.theme.icon_provider import IconProvider
from chemuson.gui.theme.tokens import get_tokens

__all__ = [
    "ATOM_COLORS",
    "ARROW_ICON_MAP",
    "BOND_ICON_MAP",
    "GENERIC_ICON_MAP",
    "ICON_SIZE",
    "draw_arrow_icon",
    "draw_atom_icon",
    "draw_bond_icon",
    "draw_charge_icon",
    "draw_coordination_sphere_icon",
    "draw_electron_icon",
    "draw_energy_diagram_icon",
    "draw_energy_levels_icon",
    "draw_generic_icon",
    "draw_glyph_icon",
    "draw_molecular_orbital_icon",
    "draw_radical_charge_icon",
    "draw_ring_icon",
    "draw_ring_template_icon",
    "draw_wavy_anchor_icon",
    "get_benzene_icon",
    "get_double_bond_icon",
    "get_eraser_icon",
    "get_pointer_icon",
    "get_single_bond_icon",
    "icon_fill_color",
    "icon_foreground_color",
    "icon_muted_color",
    "icon_paper_color",
    "set_icon_theme",
]

# Tamaño estándar de iconos (px lógicos; el provider renderiza con HiDPI).
ICON_SIZE = 32

#: Estado del tema activo para iconos monocromos (contrato histórico).
_ICON_THEME = "light"

#: Paleta de colores para símbolos atómicos (CPK; colores de dominio).
ATOM_COLORS = {
    'C': '#333333',   # Carbon - dark gray
    'N': '#3050F8',   # Nitrogen - blue
    'O': '#FF0D0D',   # Oxygen - red
    'S': '#FFFF30',   # Sulfur - yellow
    'P': '#FF8000',   # Phosphorus - orange
    'F': '#90E050',   # Fluorine - light green
    'Cl': '#1FF01F',  # Chlorine - green
    'Br': '#A62929',  # Bromine - dark red
    'H': '#FFFFFF',   # Hydrogen - white
}

# ---------------------------------------------------------------------------
# Mapas 1:1 API histórica → SVG (inventario de la Fase 2)
# ---------------------------------------------------------------------------

#: Formas de ``draw_generic_icon`` → nombre de SVG estático.
GENERIC_ICON_MAP: dict[str, str] = {
    "pointer": "pointer",
    "eraser": "eraser",
    "pan": "pan",
    "sliders": "sliders",
    "rotate_left": "rotate-left",
    "rotate_right": "rotate-right",
    "flip_horizontal": "flip-horizontal",
    "flip_vertical": "flip-vertical",
    "zoom_in": "zoom-in",
    "zoom_out": "zoom-out",
    "chain": "chain",
    "lasso": "lasso",
    "corner": "corner",
    "frame": "frame",
    "rounded_frame": "rounded-frame",
    "tlc": "tlc",
    "electrophoresis": "electrophoresis",
    "document_new": "doc-new",
    "document_open": "doc-open",
    "document_save": "doc-save",
    "undo": "undo",
    "redo": "redo",
    "copy": "copy",
    "paste": "paste",
    "clean": "clean",
}

#: Tipos de ``draw_bond_icon`` → nombre de SVG estático.
BOND_ICON_MAP: dict[str, str] = {
    "single": "bond-single",
    "bold": "bond-bold",
    "double": "bond-double",
    "triple": "bond-triple",
    "aromatic": "bond-aromatic",
    "interaction": "bond-interaction",
    "coordination": "bond-coordination",
    "wedge": "bond-wedge",
    "hashed": "bond-hashed",
    "wavy": "bond-wavy",
    "flex": "bond-flex",
}

#: Kinds de ``draw_arrow_icon`` → nombre de SVG estático.
ARROW_ICON_MAP: dict[str, str] = {
    "forward": "arrow-forward",
    "retro": "arrow-retro",
    "both": "arrow-both",
    "equilibrium": "arrow-equilibrium",
    "forward_open": "arrow-forward-open",
    "retro_open": "arrow-retro-open",
    "both_open": "arrow-both-open",
    "equilibrium_open": "arrow-equilibrium-open",
    "forward_dashed": "arrow-forward-dashed",
    "retro_dashed": "arrow-retro-dashed",
    "both_dashed": "arrow-both-dashed",
    "equilibrium_dashed": "arrow-equilibrium-dashed",
    "line": "arrow-line",
    "line_dashed": "arrow-line-dashed",
    "retrosynthetic": "arrow-retrosynthetic",
    "curved": "arrow-curved",
    "curved_fishhook": "arrow-curved-fishhook",
}

_provider: IconProvider | None = None


def _get_provider() -> IconProvider:
    """Provider module-level (lazy: se crea al primer uso, con la app viva)."""
    global _provider
    if _provider is None:
        _provider = IconProvider()
    return _provider


# ---------------------------------------------------------------------------
# Tema de iconos y colores (contrato histórico → tokens de la Fase 1)
# ---------------------------------------------------------------------------

def set_icon_theme(theme_name: str) -> None:
    """Actualiza el tema activo usado por los iconos monocromos."""
    global _ICON_THEME
    _ICON_THEME = "dark" if theme_name == "dark" else "light"
    # La fachada conserva un provider global; sincroniza el raster cacheado
    # al DPR del monitor actual si Qt ya dispone de pantalla primaria.
    if _provider is not None:
        _provider.set_device_pixel_ratio(
            _provider.application_device_pixel_ratio()
        )


def _active_theme() -> str:
    return "dark" if _ICON_THEME == "dark" else "light"


def icon_foreground_color() -> str:
    """Color de trazo principal según el tema (token ``icon``)."""
    return str(get_tokens(_active_theme())["icon"])


def icon_muted_color() -> str:
    """Color secundario según el tema (token ``text3``)."""
    return str(get_tokens(_active_theme())["text3"])


def icon_fill_color() -> str:
    """Relleno neutro según el tema (token ``surface3``)."""
    return str(get_tokens(_active_theme())["surface3"])


def icon_paper_color() -> str:
    """Relleno claro/oscuro para detalles internos (token ``surface``)."""
    return str(get_tokens(_active_theme())["surface"])


def _blank_icon() -> QIcon:
    """``QIcon`` no nulo pero en blanco (fallback, como el pixmap vacío
    de la versión QPainter)."""
    pixmap = QPixmap(ICON_SIZE, ICON_SIZE)
    pixmap.fill(Qt.GlobalColor.transparent)
    return QIcon(pixmap)


def _static(name: str) -> QIcon:
    """Icono estático del set, teñido con el color del tema activo."""
    return _get_provider().icon(name, icon_foreground_color(), ICON_SIZE)


# ---------------------------------------------------------------------------
# API pública (mismas firmas que la versión QPainter)
# ---------------------------------------------------------------------------

def draw_generic_icon(shape: str) -> QIcon:
    """Icono de herramienta genérica (puntero, borrador, etc.).

    Args:
        shape: 'pointer', 'eraser', 'pan', 'zoom_in', 'zoom_out', 'chain',
            'lasso', 'rotate_left', 'rotate_right', 'flip_horizontal',
            'flip_vertical', 'corner', 'frame', 'rounded_frame', 'tlc',
            'electrophoresis', 'document_new', 'document_open',
            'document_save', 'undo', 'redo', 'copy', 'paste', 'clean'.

    Returns:
        ``QIcon`` (en blanco si la forma no existe en el inventario, como
        la versión anterior).
    """
    name = GENERIC_ICON_MAP.get(shape)
    if name is None:
        return _blank_icon()
    return _static(name)


def draw_bond_icon(bond_type: str = 'single') -> QIcon:
    """Icono de enlace (simple, doble, cuña, etc.).

    Args:
        bond_type: Tipo de enlace ('single', 'double', 'bold', etc.).

    Returns:
        ``QIcon`` (en blanco si el tipo no existe en el inventario).
    """
    name = BOND_ICON_MAP.get(bond_type)
    if name is None:
        return _blank_icon()
    return _static(name)


def draw_arrow_icon(kind: str = "forward") -> QIcon:
    """Flechas usadas en la paleta de anotaciones.

    Args:
        kind: 'forward', 'retro', 'both', 'equilibrium', variantes
            ``_open``/``_dashed``, 'line', 'line_dashed', 'retrosynthetic',
            'curved', 'curved_fishhook'.

    Returns:
        ``QIcon`` (línea simple si el kind no existe, como la versión
        anterior).
    """
    name = ARROW_ICON_MAP.get(kind, "arrow-line")
    return _static(name)


def draw_atom_icon(text: str, color: str = None) -> QIcon:
    """Icono de átomo con el símbolo centrado (SVG dinámico CPK).

    Args:
        text: Símbolo del elemento (p. ej., "C", "N", "Cl").
        color: Color hex; si es ``None``, usa :data:`ATOM_COLORS`.

    Returns:
        ``QIcon`` con el símbolo del elemento.
    """
    if color is None:
        color = ATOM_COLORS.get(text, '#333333')
    return _get_provider().icon_dynamic(
        "atom", icon_foreground_color(), ICON_SIZE,
        text=str(text), fill=str(color),
    )


def draw_coordination_sphere_icon(color: str = "#8D99A6") -> QIcon:
    """Icono de esfera de coordinación genérica (gradiente radial)."""
    return _get_provider().icon_dynamic(
        "sphere", icon_foreground_color(), ICON_SIZE, fill=str(color),
    )


def draw_glyph_icon(text: str, color: str = None) -> QIcon:
    """Icono tipográfico minimalista (letras, corchetes, símbolos).

    El tinte por defecto sigue el tema activo; ``color`` explícito lo
    anula (p. ej. muestras de color de la barra de texto).
    """
    tint = str(color) if color is not None else icon_foreground_color()
    font_size = 11.5 if len(str(text)) == 1 else 8.5
    return _get_provider().icon_dynamic(
        "glyph", tint, ICON_SIZE,
        label=str(text), font_size=font_size,
    )


def draw_charge_icon(sign: str) -> QIcon:
    """Icono de carga circular con signo + o -."""
    return _get_provider().icon_dynamic(
        "charge", icon_foreground_color(), ICON_SIZE, sign=str(sign),
    )


def draw_electron_icon(count: int = 1, spread: float = 6.0) -> QIcon:
    """Puntos de electrones (simple, par, etc.)."""
    return _get_provider().icon_dynamic(
        "electrons", icon_foreground_color(), ICON_SIZE,
        count=int(count), spread=float(spread),
    )


def draw_energy_diagram_icon(
    box_count: int,
    *,
    label_text: str = "",
    label_side: str = "left",
    fill_color: str = "#FFFFFF",
    stroke_visible: bool = True,
) -> QIcon:
    """Icono de cajas de configuración electrónica (1..N).

    ``fill_color`` es un color de dominio que pasa el caller (preset del
    diagrama); el trazo sigue el tema activo.
    """
    return _get_provider().icon_dynamic(
        "energy-boxes", icon_foreground_color(), ICON_SIZE,
        boxes=int(box_count),
        label=str(label_text),
        side=str(label_side),
        fill=str(fill_color),
        stroke=bool(stroke_visible),
    )


def draw_energy_levels_icon() -> QIcon:
    """Icono abreviado de escalera de niveles de energía."""
    return _static("energy-levels")


def draw_molecular_orbital_icon() -> QIcon:
    """Icono esquemático de orbital molecular."""
    return _static("molecular-orbital")


def draw_radical_charge_icon(sign: str) -> QIcon:
    """Radical (punto) con un pequeño signo de carga."""
    return _get_provider().icon_dynamic(
        "radical", icon_foreground_color(), ICON_SIZE, sign=str(sign),
    )


def draw_wavy_anchor_icon() -> QIcon:
    """Icono de ancla ondulada."""
    return _static("wavy-anchor")


def draw_ring_icon(size: int = 6, aromatic: bool = True) -> QIcon:
    """Icono de anillo con el número de lados indicado (SVG dinámico).

    Args:
        size: Número de lados del polígono.
        aromatic: Si se dibuja el círculo interno aromático.

    Returns:
        ``QIcon`` con el polígono del anillo.
    """
    return _get_provider().icon_dynamic(
        "ring", icon_foreground_color(), ICON_SIZE,
        sides=int(size), aromatic=bool(aromatic),
    )


def draw_ring_template_icon(label: str, size: int = 6) -> QIcon:
    """Anillo con etiqueta para presets de plantillas (SVG dinámico)."""
    return _get_provider().icon_dynamic(
        "ring-template", icon_foreground_color(), ICON_SIZE,
        label=str(label), sides=int(size),
    )


# Convenience functions for common icons
def get_pointer_icon() -> QIcon:
    """Atajo para icono de puntero."""
    return draw_generic_icon('pointer')

def get_eraser_icon() -> QIcon:
    """Atajo para icono de borrador."""
    return draw_generic_icon('eraser')

def get_single_bond_icon() -> QIcon:
    """Atajo para icono de enlace simple."""
    return draw_bond_icon('single')

def get_double_bond_icon() -> QIcon:
    """Atajo para icono de enlace doble."""
    return draw_bond_icon('double')

def get_benzene_icon() -> QIcon:
    """Atajo para icono de anillo bencénico."""
    return draw_ring_icon()
