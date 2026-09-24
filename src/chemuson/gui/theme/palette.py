"""Paleta Qt y fuente base construidas desde los design tokens.

Fase 1 del plan de modernización de la UI. El QPalette alinea los widgets
que heredan del sistema (menús nativos, popups, etc.) con los tokens, igual
que hace el spike aprobado; la fuente aplica la familia del sistema con
fallback (Inter/Segoe UI/SF Pro Text/Roboto/DejaVu Sans) y la base de 13 px
de PLAN.md §2.2.
"""

from __future__ import annotations

from PyQt6.QtGui import QColor, QFont, QFontDatabase, QPalette

from chemuson.gui.theme.tokens import METRICS, get_tokens

__all__ = [
    "build_qpalette",
    "pick_font_family",
    "theme_font",
]

#: Familias candidatas a la fuente base (orden de preferencia del spike).
_FONT_CANDIDATES: tuple[str, ...] = (
    "Inter",
    "Segoe UI",
    "SF Pro Text",
    "Roboto",
    "DejaVu Sans",
)


def pick_font_family(candidates: tuple[str, ...] = _FONT_CANDIDATES) -> str:
    """Devuelve la primera familia candidata instalada (o la primera del SO)."""
    available = set(QFontDatabase.families())
    for family in candidates:
        if family in available:
            return family
    return QFontDatabase().families()[0]


def build_qpalette(theme_name: str) -> QPalette:
    """Construye un :class:`QPalette` completo desde los tokens del tema."""
    t = get_tokens(theme_name)

    def c(key: str) -> QColor:
        return QColor(str(t[key]))

    palette = QPalette()
    palette.setColor(QPalette.ColorRole.Window, c("bg"))
    palette.setColor(QPalette.ColorRole.WindowText, c("text1"))
    palette.setColor(QPalette.ColorRole.Base, c("surface2"))
    palette.setColor(QPalette.ColorRole.AlternateBase, c("surface3"))
    palette.setColor(QPalette.ColorRole.Text, c("text1"))
    palette.setColor(QPalette.ColorRole.Button, c("surface"))
    palette.setColor(QPalette.ColorRole.ButtonText, c("text1"))
    palette.setColor(QPalette.ColorRole.BrightText, c("text1"))
    palette.setColor(QPalette.ColorRole.Highlight, c("accent"))
    palette.setColor(QPalette.ColorRole.HighlightedText, c("onAccent"))
    palette.setColor(QPalette.ColorRole.PlaceholderText, c("text3"))
    palette.setColor(QPalette.ColorRole.ToolTipBase, c("surface"))
    palette.setColor(QPalette.ColorRole.ToolTipText, c("text1"))
    palette.setColor(QPalette.ColorRole.Link, c("accent"))

    for group in (QPalette.ColorGroup.Disabled, QPalette.ColorGroup.Inactive):
        palette.setColor(group, QPalette.ColorRole.WindowText, c("text3"))
        palette.setColor(group, QPalette.ColorRole.Text, c("text3"))
        palette.setColor(group, QPalette.ColorRole.ButtonText, c("text3"))
        palette.setColor(group, QPalette.ColorRole.PlaceholderText, c("text3"))

    return palette


def theme_font(size_px: int | None = None, bold: bool = False) -> QFont:
    """Fuente del tema: familia con fallback y tamaño base de 13 px."""
    font = QFont()
    font.setFamily(pick_font_family())
    font.setPixelSize(METRICS["fontBase"] if size_px is None else size_px)
    font.setBold(bold)
    return font
