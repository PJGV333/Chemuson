"""
Barra de formato de texto para Chemuson.

Provee controles de fuente, tamaño, estilos y color de etiquetas.
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
from PyQt6.QtGui import QAction, QIcon, QFont, QColor, QKeySequence
from PyQt6.QtCore import pyqtSignal, Qt, QSize
from chemuson.gui.icons import draw_glyph_icon

class TextFormatToolbar(QToolBar):
    """
    Barra de herramientas para formato de texto (fuente, tamaño, estilos, color).
    """
    
    # Signal emitted when any font style changes
    # Args: font_family, font_size, is_bold, is_italic, is_underline, is_sub, is_sup, property_name
    format_changed = pyqtSignal(str, int, bool, bool, bool, bool, bool, str)
    
    # Signal emitted when color changes
    color_changed = pyqtSignal(QColor)
    alignment_changed = pyqtSignal(Qt.AlignmentFlag)
    opacity_changed = pyqtSignal(int)

    def __init__(self, parent=None):
        """Inicializa la barra de formato de texto.

        Args:
            parent: Widget padre opcional.
        """
        super().__init__("Formato de Texto", parent)
        self.setIconSize(QSize(16, 16))
        self.setMovable(False)
        self._syncing_state = False
        self._alignment_buttons: list[tuple[QToolButton, str]] = []
        
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
        self.action_bold.setShortcut(QKeySequence("Ctrl+B"))
        self.action_bold.setShortcutContext(Qt.ShortcutContext.ApplicationShortcut)
        
        self.action_italic = self._add_toggle_action("I", "Cursiva", "italic")
        font_i = QFont("Times", 10)
        font_i.setItalic(True)
        self.action_italic.setFont(font_i)
        self.action_italic.setShortcut(QKeySequence("Ctrl+I"))
        self.action_italic.setShortcutContext(Qt.ShortcutContext.ApplicationShortcut)
        
        self.action_underline = self._add_toggle_action("U", "Subrayado", "underline")
        font_u = QFont("Times", 10)
        font_u.setUnderline(True)
        self.action_underline.setFont(font_u)
        self.action_underline.setShortcut(QKeySequence("Ctrl+U"))
        self.action_underline.setShortcutContext(Qt.ShortcutContext.ApplicationShortcut)
        
        self.addSeparator()

        # --- Sub / Sup ---
        self.action_sub = self._add_toggle_action("x₂", "Subíndice", "sub")
        self.action_sup = self._add_toggle_action("x²", "Superíndice", "sup")
        self.action_sub.setShortcut(QKeySequence("Ctrl+="))
        self.action_sub.setShortcutContext(Qt.ShortcutContext.ApplicationShortcut)
        self.action_sup.setShortcuts([QKeySequence("Ctrl+Shift+="), QKeySequence("Ctrl++")])
        self.action_sup.setShortcutContext(Qt.ShortcutContext.ApplicationShortcut)
        
        # Ensure sub/sup are exclusive-ish
        self.action_sub.triggered.connect(lambda c: self._handle_exclusive(self.action_sub, self.action_sup, "sub"))
        self.action_sup.triggered.connect(lambda c: self._handle_exclusive(self.action_sup, self.action_sub, "sup"))

        self.addSeparator()

        # --- Alignment ---
        self._add_align_button("format-justify-left", "Alinear a la izquierda", "≡", Qt.AlignmentFlag.AlignLeft)
        self._add_align_button("format-justify-center", "Centrar", "≣", Qt.AlignmentFlag.AlignHCenter)
        self._add_align_button("format-justify-fill", "Justificar", "☰", Qt.AlignmentFlag.AlignJustify)

        self.addSeparator()
        
        # --- Color ---
        self.color_btn = QToolButton()
        self.color_btn.setToolTip("Color de texto")
        self.color_btn.clicked.connect(self._pick_color)
        self._current_color = QColor(Qt.GlobalColor.black)
        self._update_color_icon()
        self.addWidget(self.color_btn)

        self.addSeparator()

        # --- Opacity ---
        self.opacity_spin = QSpinBox()
        self.opacity_spin.setRange(0, 100)
        self.opacity_spin.setSingleStep(5)
        self.opacity_spin.setValue(100)
        self.opacity_spin.setSuffix("%")
        self.opacity_spin.setPrefix("Op ")
        self.opacity_spin.setToolTip(
            "Opacidad. Sin selección aplica a todo el canvas; con selección, solo a los elementos seleccionados."
        )
        self.opacity_spin.valueChanged.connect(self._emit_opacity_change)
        self.addWidget(self.opacity_spin)

    def _add_toggle_action(self, text: str, tooltip: str, prop_name: str) -> QAction:
        """Crea una acción conmutadora (toggle) para formato.

        Args:
            text: Texto mostrado en el botón.
            tooltip: Texto de ayuda.
            prop_name: Nombre de la propiedad asociada.

        Returns:
            Acción Qt configurada.
        """
        action = QAction(text, self)
        action.setToolTip(tooltip)
        action.setCheckable(True)
        action.triggered.connect(lambda: self._emit_change(prop_name))
        self.addAction(action)
        return action

    def _add_align_button(
        self,
        theme_name: str,
        tooltip: str,
        fallback_glyph: str,
        alignment: Qt.AlignmentFlag,
    ) -> None:
        """Añade un botón de alineación con icono temático o de respaldo."""
        button = QToolButton(self)
        button.setIcon(draw_glyph_icon(fallback_glyph))
        button.setToolTip(tooltip)
        button.setAutoRaise(True)
        button.clicked.connect(lambda: self.alignment_changed.emit(alignment))
        self.addWidget(button)
        self._alignment_buttons.append((button, fallback_glyph))

    def _handle_exclusive(self, trigger_action: QAction, other_action: QAction, prop_name: str):
        """Fuerza exclusividad entre subíndice y superíndice."""
        if trigger_action.isChecked():
            other_action.setChecked(False)
        self._emit_change(prop_name)

    def _pick_color(self):
        """Abre un diálogo para seleccionar el color de texto."""
        color = QColorDialog.getColor(self._current_color, self, "Seleccionar color de texto")
        if color.isValid():
            self._current_color = color
            self._update_color_icon()
            self.color_changed.emit(color)

    def _update_color_icon(self):
        """Actualiza el icono de color actual."""
        # Crear un icono simple con el color actual
        pixmap = draw_glyph_icon("■", color=self._current_color.name())
        self.color_btn.setIcon(pixmap)

    def refresh_icons(self) -> None:
        """Regenera iconos dependientes del tema."""
        for button, glyph in self._alignment_buttons:
            button.setIcon(draw_glyph_icon(glyph))
        self._update_color_icon()

    def _emit_change(self, property_name: str = "all"):
        """Emite la señal con el estado actual de formato."""
        if self._syncing_state:
            return
        # Recolectar estado actual y emitir señal.
        family = self.font_combo.currentFont().family()
        size = self.size_spin.value()
        b = self.action_bold.isChecked()
        i = self.action_italic.isChecked()
        u = self.action_underline.isChecked()
        sub = self.action_sub.isChecked()
        sup = self.action_sup.isChecked()
        
        self.format_changed.emit(family, size, b, i, u, sub, sup, property_name)

    def _emit_opacity_change(self, value: int) -> None:
        """Emite cambios del control de opacidad evitando realimentación."""
        if self._syncing_state:
            return
        self.opacity_changed.emit(int(value))

    def set_state(self, font: QFont, settings: dict):
        """Actualiza la barra según el estado de selección actual."""
        self._syncing_state = True
        try:
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
        finally:
            self._syncing_state = False

    def set_opacity_percent(self, value: int) -> None:
        """Sincroniza el control de opacidad sin reemitir la señal."""
        self._syncing_state = True
        try:
            self.opacity_spin.setValue(max(0, min(100, int(value))))
        finally:
            self._syncing_state = False
