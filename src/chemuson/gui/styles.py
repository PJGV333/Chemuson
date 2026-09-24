"""
Estilos visuales de la aplicación Chemuson (fachada de compatibilidad).

Este módulo conserva la API pública histórica para los callers existentes
(``main_window.py`` y ``toolbar.py``), pero la fuente de verdad visual es el
sistema de design tokens de :mod:`chemuson.gui.theme` (Fase 1 del plan de
modernización de la UI; OpenSpec
``2026-09-24-modernize-ui-theme-foundation``).

Los generadores de QSS son los del nuevo sistema; ``LIGHT_COLORS`` y
``DARK_COLORS`` se conservan como alias *deprecated* de los tokens (mapeo
1:1 de clave legada → token). Para código nuevo, usa
``chemuson.gui.theme`` directamente.
"""

from chemuson.gui.theme import (
    DEFAULT_THEME_NAME,
    get_main_stylesheet,
    get_tokens,
    get_tool_palette_stylesheet,
)

__all__ = [
    "DARK_COLORS",
    "DEFAULT_THEME",
    "LIGHT_COLORS",
    "MAIN_STYLESHEET",
    "TOOL_PALETTE_STYLESHEET",
    "get_main_stylesheet",
    "get_tool_palette_stylesheet",
]

#: Tema por defecto (API histórica; igual a ``theme.DEFAULT_THEME_NAME``).
DEFAULT_THEME = DEFAULT_THEME_NAME

# Mapeo de claves de la paleta legada a tokens actuales.
_LEGACY_KEY_TO_TOKEN = {
    "primary_dark": "surface",
    "primary_medium": "borderStrong",
    "accent_primary": "accent",
    "accent_hover": "accentHover",
    "accent_pressed": "accentStrong",
    "bg_main": "bg",
    "bg_elevated": "surface",
    "bg_toolbar": "surface2",
    "bg_dock": "surface2",
    "border_light": "border",
    "border_medium": "borderStrong",
    "border_dark": "borderStrong",
    "text_primary": "text1",
    "text_secondary": "text2",
    "text_muted": "text3",
    "text_inverse": "onAccent",
    "palette_bg": "surface2",
    "palette_border": "border",
    "palette_button_bg": "surface",
    "palette_button_border": "borderStrong",
    "palette_button_hover": "surface2",
    "palette_selected_bg": "accentSoft",
    "palette_selected_border": "accentBorder",
}


def _legacy_palette(theme_name: str) -> dict[str, str]:
    """Paleta legada como alias de los tokens del tema (deprecated)."""
    tokens = get_tokens(theme_name)
    return {legacy: str(tokens[token]) for legacy, token in _LEGACY_KEY_TO_TOKEN.items()}


#: Colores Modo Claro (deprecated: alias de los tokens ``light``).
LIGHT_COLORS = _legacy_palette("light")

#: Colores Modo Oscuro (deprecated: alias de los tokens ``dark``).
DARK_COLORS = _legacy_palette("dark")

# Hojas de estilo generadas en el import con el tema por defecto
# (comportamiento histórico de este módulo).
MAIN_STYLESHEET = get_main_stylesheet(DEFAULT_THEME)
TOOL_PALETTE_STYLESHEET = get_tool_palette_stylesheet(DEFAULT_THEME)
