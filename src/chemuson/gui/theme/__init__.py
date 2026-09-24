"""Sistema de temas de la UI de Chemuson (fuente de verdad visual).

Fase 1 del plan de modernización de la UI
(``docs/ui-modernization/PLAN.md``). Este paquete centraliza:

- :mod:`~chemuson.gui.theme.tokens` — design tokens light/dark + métricas.
- :mod:`~chemuson.gui.theme.qss` — generadores de hojas de estilo.
- :mod:`~chemuson.gui.theme.palette` — ``QPalette`` y fuente base.
- :mod:`~chemuson.gui.theme.icon_provider` — provider SVG→QIcon (Fase 2).

La fachada de compatibilidad histórica es ``chemuson.gui.styles``.
"""

from __future__ import annotations

from chemuson.gui.theme.palette import build_qpalette, theme_font
from chemuson.gui.theme.qss import get_main_stylesheet, get_tool_palette_stylesheet
from chemuson.gui.theme.tokens import (
    DARK_TOKENS,
    DEFAULT_THEME_NAME,
    LIGHT_TOKENS,
    METRICS,
    SYSTEM_THEME_NAME,
    THEME_NAMES,
    get_tokens,
    theme_color,
)

__all__ = [
    "DARK_TOKENS",
    "DEFAULT_THEME_NAME",
    "LIGHT_TOKENS",
    "METRICS",
    "SYSTEM_THEME_NAME",
    "THEME_NAMES",
    "apply_theme",
    "build_qpalette",
    "get_main_stylesheet",
    "get_tokens",
    "get_tool_palette_stylesheet",
    "resolve_theme_name",
    "set_theme_from_system",
    "system_theme_name",
    "theme_color",
    "theme_font",
]


def system_theme_name() -> str:
    """Devuelve ``"light"`` o ``"dark"`` según el esquema del sistema.

    Usa ``QStyleHints.colorScheme()`` (Qt >= 6.5) cuando hay una
    ``QApplication`` y el esquema es explícito; si no, cae a
    :data:`DEFAULT_THEME_NAME`.
    """
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import QApplication

    app = QApplication.instance()
    if app is not None:
        hints = app.styleHints()
        if hints is not None and hasattr(hints, "colorScheme"):
            scheme = hints.colorScheme()
            if scheme == Qt.ColorScheme.Dark:
                return "dark"
            if scheme == Qt.ColorScheme.Light:
                return "light"
    return DEFAULT_THEME_NAME


def resolve_theme_name(name: str) -> str:
    """Normaliza un nombre de tema a ``"light"`` o ``"dark"``.

    - ``"light"``/``"dark"`` (mayúsculas o espacios): resuelven a sí mismos.
    - ``"system"``: resuelve al esquema actual del sistema operativo.
    - Cualquier otro valor: :data:`DEFAULT_THEME_NAME` (sin excepciones).
    """
    if not isinstance(name, str):
        return DEFAULT_THEME_NAME
    normalized = name.strip().lower()
    if normalized in THEME_NAMES:
        return normalized
    if normalized == SYSTEM_THEME_NAME:
        return system_theme_name()
    return DEFAULT_THEME_NAME


def apply_theme(target: object, theme_name: str) -> str:
    """Aplica el tema (fuente + ``QPalette`` + QSS principal) a un destino.

    Args:
        target: ``QApplication`` o ``QWidget`` (p. ej. la ventana principal).
        theme_name: ``"light"``, ``"dark"`` o ``"system"`` (ver
            :func:`resolve_theme_name`).

    Returns:
        El nombre de tema resuelto que se aplicó (``"light"`` o ``"dark"``).
    """
    resolved = resolve_theme_name(theme_name)
    target.setFont(theme_font())
    target.setPalette(build_qpalette(resolved))
    target.setStyleSheet(get_main_stylesheet(resolved))
    return resolved


def set_theme_from_system(target: object) -> str:
    """Aplica a ``target`` el tema que coincide con el esquema del sistema.

    Resolución en el momento de la llamada (una-shot); Qt no expone una
    señal portable de cambio de esquema, por lo que un re-apply futuro se
    integrará con la fase de ajustes.
    """
    return apply_theme(target, SYSTEM_THEME_NAME)
