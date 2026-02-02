"""
Text Formatting Toolbar
Provides controls for font settings, text styles, and colors.
"""
from __future__ import annotations

from PyQt6.QtWidgets import (
    QToolBar,
    QFontComboBox,
    QSpinBox,
    QToolButton,
    QWidget,
    QHBoxLayout,
    QColorDialog,
)
from PyQt6.QtGui import QAction, QIcon, QFont, QColor
from PyQt6.QtCore import pyqtSignal, Qt, QSize
from gui.icons import draw_glyph_icon

class TextFormatToolbar(QToolBar):
    """
    Toolbar for text formatting options (Font, Size, B/I/U, Color).
    """
    
    # Signal emitted when any font style changes
    # Args: font_family, font_size, is_bold, is_italic, is_underline, is_sub, is_sup, property_name
    format_changed = pyqtSignal(str, int, bool, bool, bool, bool, bool, str)
    
    # Signal emitted when color changes
    color_changed = pyqtSignal(QColor)

    def __init__(self, parent=None):
        super().__init__("Formato de Texto", parent)
        self.setIconSize(QSize(16, 16))
        self.setMovable(False)
        
        # --- Font Family ---
        self.font_combo = QFontComboBox()
        self.font_combo.setCurrentFont(QFont("Arial"))
        self.font_combo.currentFontChanged.connect(lambda: self._emit_change("family"))
        self.addWidget(self.font_combo)
        
        self.addSeparator()

        # --- Font Size ---
        self.size_spin = QSpinBox()
        self.size_spin.setRange(6, 144)
        self.size_spin.setValue(12)
        self.size_spin.valueChanged.connect(lambda: self._emit_change("size"))
        self.addWidget(self.size_spin)
        
        self.addSeparator()

        # --- B / I / U ---
        self.action_bold = self._add_toggle_action("B", "Negrita", "bold")
        self.action_bold.setFont(QFont("Times", 10, QFont.Weight.Bold))
        
        self.action_italic = self._add_toggle_action("I", "Cursiva", "italic")
        font_i = QFont("Times", 10)
        font_i.setItalic(True)
        self.action_italic.setFont(font_i)
        
        self.action_underline = self._add_toggle_action("U", "Subrayado", "underline")
        font_u = QFont("Times", 10)
        font_u.setUnderline(True)
        self.action_underline.setFont(font_u)
        
        self.addSeparator()

        # --- Sub / Sup ---
        self.action_sub = self._add_toggle_action("x₂", "Subíndice", "sub")
        self.action_sup = self._add_toggle_action("x²", "Superíndice", "sup")
        
        # Ensure sub/sup are exclusive-ish
        self.action_sub.triggered.connect(lambda c: self._handle_exclusive(self.action_sub, self.action_sup, "sub"))
        self.action_sup.triggered.connect(lambda c: self._handle_exclusive(self.action_sup, self.action_sub, "sup"))

        self.addSeparator()
        
        # --- Color ---
        self.color_btn = QToolButton()
        self.color_btn.setToolTip("Color de texto")
        self.color_btn.clicked.connect(self._pick_color)
        self._current_color = QColor(Qt.GlobalColor.black)
        self._update_color_icon()
        self.addWidget(self.color_btn)

    def _add_toggle_action(self, text: str, tooltip: str, prop_name: str) -> QAction:
        action = QAction(text, self)
        action.setToolTip(tooltip)
        action.setCheckable(True)
        action.triggered.connect(lambda: self._emit_change(prop_name))
        self.addAction(action)
        return action

    def _handle_exclusive(self, trigger_action: QAction, other_action: QAction, prop_name: str):
        if trigger_action.isChecked():
            other_action.setChecked(False)
        self._emit_change(prop_name)

    def _pick_color(self):
        color = QColorDialog.getColor(self._current_color, self, "Seleccionar color de texto")
        if color.isValid():
            self._current_color = color
            self._update_color_icon()
            self.color_changed.emit(color)

    def _update_color_icon(self):
        # Create a simple icon with the current color
        pixmap = draw_glyph_icon("■", color=self._current_color.name())
        self.color_btn.setIcon(pixmap)

    def _emit_change(self, property_name: str = "all"):
        # Gather all states and emit
        family = self.font_combo.currentFont().family()
        size = self.size_spin.value()
        b = self.action_bold.isChecked()
        i = self.action_italic.isChecked()
        u = self.action_underline.isChecked()
        sub = self.action_sub.isChecked()
        sup = self.action_sup.isChecked()
        
        self.format_changed.emit(family, size, b, i, u, sub, sup, property_name)
    
    def set_state(self, font: QFont, settings: dict):
        """Update toolbar to reflect current selection state."""
        self.blockSignals(True)
        self.font_combo.setCurrentFont(font)
        self.size_spin.setValue(max(6, int(font.pointSizeF())))
        self.action_bold.setChecked(font.bold())
        self.action_italic.setChecked(font.italic())
        self.action_underline.setChecked(font.underline())
        
        # Check specific flags/settings for sub/sup if available
        is_sub = settings.get("sub", False)
        is_sup = settings.get("sup", False)
        self.action_sub.setChecked(is_sub)
        self.action_sup.setChecked(is_sup)
        
        color = settings.get("color")
        if color and isinstance(color, QColor):
            self._current_color = color
            self._update_color_icon()
            
        self.blockSignals(False)

