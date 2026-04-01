"""Proveedor de iconos con fallback y contraste consistente por tema."""

from __future__ import annotations

from PyQt6.QtGui import QColor, QIcon, QPainter, QPixmap

from chemuson.gui.icons import draw_generic_icon


def _tint_icon(icon: QIcon, color: str, size: int = 24) -> QIcon:
    base = icon.pixmap(size, size)
    if base.isNull():
        return icon
    tinted = QPixmap(base.size())
    tinted.fill(QColor("transparent"))
    painter = QPainter(tinted)
    painter.drawPixmap(0, 0, base)
    painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_SourceIn)
    painter.fillRect(tinted.rect(), QColor(color))
    painter.end()
    return QIcon(tinted)


def themed_icon(theme_name: str, fallback_shape: str, *, is_dark: bool) -> QIcon:
    icon = QIcon.fromTheme(theme_name)
    if icon.isNull():
        return draw_generic_icon(fallback_shape)
    tint = "#EAF2FF" if is_dark else "#1E293B"
    return _tint_icon(icon, tint)


def apply_main_action_icons(window, *, theme: str) -> None:
    is_dark = theme == "dark"
    window.action_new.setIcon(themed_icon("document-new", "chain", is_dark=is_dark))
    window.action_open.setIcon(themed_icon("document-open", "pointer", is_dark=is_dark))
    window.action_save.setIcon(themed_icon("document-save", "pointer", is_dark=is_dark))
    window.action_undo.setIcon(themed_icon("edit-undo", "rotate_left", is_dark=is_dark))
    window.action_redo.setIcon(themed_icon("edit-redo", "rotate_right", is_dark=is_dark))
    window.action_copy.setIcon(themed_icon("edit-copy", "pointer", is_dark=is_dark))
    window.action_paste.setIcon(themed_icon("edit-paste", "frame", is_dark=is_dark))
    window.action_branch_auto_arrange.setIcon(themed_icon("edit-clear", "eraser", is_dark=is_dark))
    window.action_clean_2d.setIcon(themed_icon("edit-clear", "eraser", is_dark=is_dark))
