"""Design tokens de la UI de Chemuson (fuente de verdad visual).

Fase 1 del plan de modernización de la UI (``docs/ui-modernization/PLAN.md``).

Los valores del tema copian la tabla de tokens del spike PyQt6 aprobado
(``docs/ui-modernization/pyqt6-spike/theme.py``, commit ``59e977d``,
``SMOKE: OK — 0 fallos``), que a su vez traduce el mockup
(``docs/ui-modernization/mockup-ui.html``) y la tabla §2.2 de PLAN.md.
Las métricas (spacing/radios/tipografía) completan lo que el spike dejaba
en constants sueltas, según la grilla de 8 px de PLAN.md §2.2.

Solo este módulo (y los generadores de ``qss.py``) conocen colores de UI;
el resto de la UI los consume a través de ``get_tokens``/``theme_color``.
Los colores químicos de dominio (CPK, hoja del canvas) viven en sus
subsistemas propios y no se tocan en esta fase.
"""

from __future__ import annotations

from PyQt6.QtGui import QColor

__all__ = [
    "DARK_TOKENS",
    "DEFAULT_THEME_NAME",
    "LIGHT_TOKENS",
    "METRICS",
    "SYSTEM_THEME_NAME",
    "THEME_NAMES",
    "get_tokens",
    "theme_color",
]

# ---------------------------------------------------------------------------
# Nombres de tema
# ---------------------------------------------------------------------------

#: Temas resueltos disponibles (la API acepta además ``system``).
THEME_NAMES: tuple[str, ...] = ("light", "dark")

#: Valor por defecto cuando no hay preferencia persistida.
DEFAULT_THEME_NAME: str = "light"

#: Nombre lógico "seguir sistema": se resuelve a light/dark en
#: ``chemuson.gui.theme.resolve_theme_name`` (preparado, aún sin opción en la
#: UI; ver OpenSpec ``2026-09-24-modernize-ui-theme-foundation``).
SYSTEM_THEME_NAME: str = "system"

# ---------------------------------------------------------------------------
# Tokens de color (misma tabla que el spike aprobado)
# ---------------------------------------------------------------------------

LIGHT_TOKENS: dict[str, object] = {
    "bg": "#F1F5F9",
    "surface": "#FFFFFF",
    "surface2": "#F8FAFC",
    "surface3": "#EEF2F7",
    "border": "#E2E8F0",
    "borderStrong": "#CBD5E1",
    "text1": "#0F172A",
    "text2": "#475569",
    "text3": "#94A3B8",
    "accent": "#0E7490",
    "accentStrong": "#155E75",
    "accentHover": "#0891B2",
    "accentSoft": "#ECFEFF",
    "accentBorder": "#A5F3FC",
    "onAccent": "#FFFFFF",
    "danger": "#DC2626",
    "dangerSoft": "#FEF2F2",
    "warn": "#D97706",
    "warnSoft": "#FFFBEB",
    "ok": "#059669",
    "okSoft": "#ECFDF5",
    # Tokens de dominio (canvas/hoja): se incluyen para que las fases de
    # canvas/pulido los consuman; esta fase no los aplica al canvas.
    #
    # Nota (fidelidad vs spike): el spike escribía los rellenos suaves del
    # tema oscuro como ``rgba(r,g,b,a)``; aquí se usan hex de 8 dígitos
    # ``#AARRGGBB`` (mismo color exacto) porque Qt los acepta tanto en QSS
    # como en el constructor de ``QColor`` (``QColor`` no parsea ``rgba()``),
    # lo que mantiene la invarianta "todo token es un color válido".
    "sheet": "#FFFFFF",
    "sheetGrid": "#E9EEF4",
    "canvasBg": "#E9EEF4",
    # Sombras (QGraphicsDropShadowEffect; QSS no tiene box-shadow).
    "shadow1": QColor(15, 23, 42, 25),
    "shadow2": QColor(15, 23, 42, 38),
    "sheetShadow": QColor(15, 23, 42, 31),
    # Color de tinte de iconos por estado (lo consumirá la Fase 2).
    "icon": "#475569",
    "iconHover": "#0F172A",
    "iconActive": "#155E75",
}

DARK_TOKENS: dict[str, object] = {
    "bg": "#0B1120",
    "surface": "#0F172A",
    "surface2": "#16213A",
    "surface3": "#1E293B",
    "border": "#263349",
    "borderStrong": "#334155",
    "text1": "#F1F5F9",
    "text2": "#C3CEDF",
    "text3": "#64748B",
    "accent": "#22D3EE",
    "accentStrong": "#67E8F9",
    "accentHover": "#06B6D4",
    "accentSoft": "#1F22D3EE",
    "accentBorder": "#7322D3EE",
    "onAccent": "#083344",
    "danger": "#F87171",
    "dangerSoft": "#21F87171",
    "warn": "#FBBF24",
    "warnSoft": "#21FBBF24",
    "ok": "#34D399",
    "okSoft": "#2134D399",
    "sheet": "#FFFFFF",
    "sheetGrid": "#E9EEF4",
    "canvasBg": "#0A0F1A",
    "shadow1": QColor(0, 0, 0, 115),
    "shadow2": QColor(0, 0, 0, 140),
    "sheetShadow": QColor(0, 0, 0, 140),
    "icon": "#C3CEDF",
    "iconHover": "#F1F5F9",
    "iconActive": "#67E8F9",
}

_TOKEN_TABLES: dict[str, dict[str, object]] = {
    "light": LIGHT_TOKENS,
    "dark": DARK_TOKENS,
}

# ---------------------------------------------------------------------------
# Métricas UI (grilla 8 px, PLAN.md §2.2; valores en px)
# ---------------------------------------------------------------------------

METRICS: dict[str, int] = {
    # Espaciado
    "spacingXs": 4,
    "spacingSm": 8,
    "spacingMd": 12,
    "spacingLg": 16,
    # Radios
    "radiusSurface": 10,
    "radiusBtn": 9,
    "radiusChip": 8,
    # Tipografía
    "fontBase": 13,
    "fontSmall": 12,
    "fontTiny": 11,
    # Fase 3 (shell superior)
    "appbarH": 54,
    # Fase 4 (rail de herramientas + flyouts)
    "railW": 58,
    "railBtn": 44,
    "flyoutW": 244,
}


def get_tokens(theme_name: str) -> dict[str, object]:
    """Devuelve la tabla de tokens del tema solicitado.

    Args:
        theme_name: ``"light"`` o ``"dark"``. Cualquier otro valor resuelve
            a :data:`DEFAULT_THEME_NAME` (sin excepciones).

    Returns:
        Tabla de tokens (inmutable por convención: no mutar in-place).
    """
    return _TOKEN_TABLES.get(theme_name, _TOKEN_TABLES[DEFAULT_THEME_NAME])


def theme_color(theme_name: str, key: str) -> QColor:
    """Devuelve el token ``key`` del tema como :class:`QColor`.

    Los tokens que ya son ``QColor`` (sombras) se devuelven como copia.
    """
    value = get_tokens(theme_name)[key]
    if isinstance(value, QColor):
        return QColor(value)
    return QColor(str(value))
