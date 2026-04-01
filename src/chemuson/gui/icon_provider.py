"""Proveedor de iconos con fallback y contraste consistente por tema."""

from __future__ import annotations

from PyQt6.QtCore import QRect, Qt
from PyQt6.QtGui import QColor, QIcon, QPainter, QPixmap

from chemuson.gui.icons import draw_generic_icon


def _alpha_bounding_rect(pixmap: QPixmap) -> QRect:
    """Calcula el bounding rect de píxeles no transparentes."""
    image = pixmap.toImage()
    width = image.width()
    height = image.height()
    min_x, min_y = width, height
    max_x = max_y = -1
    for y in range(height):
        for x in range(width):
            if image.pixelColor(x, y).alpha() > 0:
                min_x = min(min_x, x)
                min_y = min(min_y, y)
                max_x = max(max_x, x)
                max_y = max(max_y, y)
    if max_x < min_x or max_y < min_y:
        return QRect(0, 0, width, height)
    return QRect(min_x, min_y, max_x - min_x + 1, max_y - min_y + 1)


def _tint_icon(icon: QIcon, color: str, size: int = 24) -> QIcon:
    request = max(48, size * 3)
    base = icon.pixmap(request, request)
    if base.isNull():
        return icon
    alpha_rect = _alpha_bounding_rect(base)
    glyph = base.copy(alpha_rect)
    target_px = int(size * 0.96)
    scaled = glyph.scaled(
        target_px,
        target_px,
        Qt.AspectRatioMode.KeepAspectRatio,
        Qt.TransformationMode.SmoothTransformation,
    )
    canvas = QPixmap(size, size)
    canvas.fill(QColor("transparent"))
    x = (size - scaled.width()) // 2
    y = (size - scaled.height()) // 2
    painter = QPainter(canvas)
    painter.drawPixmap(x, y, scaled)
    painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_SourceIn)
    painter.fillRect(canvas.rect(), QColor(color))
    painter.end()
    return QIcon(canvas)


def retint_icon(icon: QIcon, color: str, size: int) -> QIcon:
    """Retinta un icono existente respetando el tamaño objetivo."""
    return _tint_icon(icon, color, size=size)


def themed_icon(theme_name: str, fallback_shape: str, *, is_dark: bool, size: int = 24) -> QIcon:
    icon = QIcon.fromTheme(theme_name)
    if icon.isNull():
        return draw_generic_icon(fallback_shape)
    tint = "#EAF2FF" if is_dark else "#1E293B"
    return _tint_icon(icon, tint, size=size)


def apply_main_action_icons(window, *, theme: str) -> None:
    is_dark = theme == "dark"
    icon_size = 24
    if hasattr(window, "main_toolbar"):
        icon_size = max(24, int(window.main_toolbar.iconSize().width()))
    window.action_new.setIcon(themed_icon("document-new", "chain", is_dark=is_dark, size=icon_size))
    window.action_open.setIcon(themed_icon("document-open", "pointer", is_dark=is_dark, size=icon_size))
    window.action_save.setIcon(themed_icon("document-save", "pointer", is_dark=is_dark, size=icon_size))
    window.action_undo.setIcon(themed_icon("edit-undo", "rotate_left", is_dark=is_dark, size=icon_size))
    window.action_redo.setIcon(themed_icon("edit-redo", "rotate_right", is_dark=is_dark, size=icon_size))
    window.action_copy.setIcon(themed_icon("edit-copy", "pointer", is_dark=is_dark, size=icon_size))
    window.action_paste.setIcon(themed_icon("edit-paste", "frame", is_dark=is_dark, size=icon_size))
    window.action_branch_auto_arrange.setIcon(themed_icon("edit-clear", "eraser", is_dark=is_dark, size=icon_size))
    window.action_clean_2d.setIcon(themed_icon("edit-clear", "eraser", is_dark=is_dark, size=icon_size))
