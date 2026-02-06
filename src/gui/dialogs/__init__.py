"""
Diálogos de Chemuson.

Incluye preferencias, estilos y ayudas de inicio rápido.
"""
from __future__ import annotations

from dataclasses import replace

from PyQt6.QtCore import pyqtSignal, Qt
from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QColorDialog,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QTabWidget,
    QTextBrowser,
    QVBoxLayout,
    QWidget,
)

from core.model import ChemState
from gui.style import DrawingStyle


class PreferencesDialog(QDialog):
    """Preferencias de la aplicación en pestañas."""

    preferences_changed = pyqtSignal(dict)

    def __init__(self, current_state: ChemState, current_style: DrawingStyle, parent=None) -> None:
        """Inicializa el diálogo.

        Args:
            current_state: Descripción del parámetro.
            current_style: Descripción del parámetro.
            parent: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado del diálogo.
        """
        super().__init__(parent)
        self.setWindowTitle("Preferencias")
        self.setMinimumWidth(420)
        self._current_style = current_style

        tabs = QTabWidget(self)
        tabs.addTab(self._build_general_tab(current_state), "General")
        tabs.addTab(self._build_appearance_tab(), "Apariencia")
        tabs.addTab(self._build_rdkit_tab(), "RDKit")

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self._on_accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addWidget(tabs)
        layout.addWidget(buttons)

        self._show_carbons = current_state.show_implicit_carbons
        self._show_hydrogens = current_state.show_implicit_hydrogens
        self._use_aromatic_circles = current_state.use_aromatic_circles

    def _build_general_tab(self, current_state: ChemState) -> QWidget:
        """Construye general tab.

        Args:
            current_state: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado del diálogo.
        """
        widget = QWidget()
        layout = QVBoxLayout(widget)

        self.carbons_checkbox = QCheckBox("Mostrar carbonos implícitos")
        self.carbons_checkbox.setChecked(current_state.show_implicit_carbons)
        layout.addWidget(self.carbons_checkbox)

        self.hydrogens_checkbox = QCheckBox("Mostrar hidrógenos implícitos")
        self.hydrogens_checkbox.setChecked(current_state.show_implicit_hydrogens)
        layout.addWidget(self.hydrogens_checkbox)

        self.aromatic_checkbox = QCheckBox("Aromáticos como círculos")
        self.aromatic_checkbox.setChecked(current_state.use_aromatic_circles)
        layout.addWidget(self.aromatic_checkbox)

        layout.addStretch()
        return widget

    def _build_appearance_tab(self) -> QWidget:
        """Construye appearance tab.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado del diálogo.
        """
        widget = QWidget()
        layout = QVBoxLayout(widget)
        form = QFormLayout()

        self.bond_cap_combo = QComboBox()
        self.bond_cap_combo.addItem("Redondeados", "round")
        self.bond_cap_combo.addItem("Rectos", "flat")
        current_cap = (
            "round"
            if self._current_style.cap_style == Qt.PenCapStyle.RoundCap
            else "flat"
        )
        index = self.bond_cap_combo.findData(current_cap)
        if index >= 0:
            self.bond_cap_combo.setCurrentIndex(index)
        form.addRow("Bordes de enlace", self.bond_cap_combo)

        layout.addLayout(form)
        layout.addStretch()
        return widget

    def _build_rdkit_tab(self) -> QWidget:
        """Construye rdkit tab.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado del diálogo.
        """
        widget = QWidget()
        layout = QVBoxLayout(widget)
        label = QLabel("Preferencias de RDKit (longitud de enlace, limpieza 2D) próximamente.")
        label.setStyleSheet("color: #666666;")
        layout.addWidget(label)
        layout.addStretch()
        return widget

    def _on_accept(self) -> None:
        """Maneja accept.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado del diálogo.
        """
        prefs = {
            "show_carbons": self.carbons_checkbox.isChecked(),
            "show_hydrogens": self.hydrogens_checkbox.isChecked(),
            "aromatic_circles": self.aromatic_checkbox.isChecked(),
            "bond_caps": self.bond_cap_combo.currentData(),
        }
        self.preferences_changed.emit(prefs)
        self.accept()


class StyleDialog(QDialog):
    """Diálogo para editar propiedades del estilo de dibujo."""

    def __init__(self, current_style: DrawingStyle, bond_length: float, parent=None) -> None:
        """Inicializa el diálogo.

        Args:
            current_style: Descripción del parámetro.
            bond_length: Descripción del parámetro.
            parent: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado del diálogo.
        """
        super().__init__(parent)
        self.setWindowTitle("Estilo de dibujo")
        self.setMinimumWidth(360)

        self._style = current_style
        self._bond_length = bond_length
        self._bond_color = QColor(current_style.bond_color)
        self._atom_fill_color = QColor(current_style.atom_fill_color)
        self._atom_stroke_color = QColor(current_style.atom_stroke_color)
        self._result_style: DrawingStyle | None = None
        self._result_bond_length: float | None = None

        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignmentFlag.AlignLeft)

        self.stroke_spin = QDoubleSpinBox()
        self.stroke_spin.setRange(0.5, 8.0)
        self.stroke_spin.setSingleStep(0.1)
        self.stroke_spin.setValue(current_style.stroke_px)
        form.addRow("Grosor de línea", self.stroke_spin)

        self.length_spin = QDoubleSpinBox()
        self.length_spin.setRange(10.0, 120.0)
        self.length_spin.setSingleStep(1.0)
        self.length_spin.setValue(bond_length)
        form.addRow("Longitud de enlace", self.length_spin)

        bond_row, self.bond_color_btn = self._build_color_row(self._bond_color)
        form.addRow("Color de enlaces", bond_row)

        atom_fill_row, self.atom_fill_btn = self._build_color_row(self._atom_fill_color)
        form.addRow("Relleno de vértices", atom_fill_row)

        atom_stroke_row, self.atom_stroke_btn = self._build_color_row(self._atom_stroke_color)
        form.addRow("Borde de vértices", atom_stroke_row)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self._on_accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(buttons)

    def _build_color_row(self, color: QColor) -> tuple[QWidget, QPushButton]:
        """Construye color row.

        Args:
            color: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado del diálogo.
        """
        button = QPushButton()
        button.setFixedWidth(60)
        self._set_button_color(button, color)
        button.clicked.connect(lambda: self._pick_color(button))

        row = QWidget()
        layout = QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(button)
        layout.addStretch()
        return row, button

    def _set_button_color(self, button: QPushButton, color: QColor) -> None:
        """Configura button color.

        Args:
            button: Descripción del parámetro.
            color: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado del diálogo.
        """
        button.setStyleSheet(
            f"background-color: {color.name()}; border: 1px solid #666666;")

    def _pick_color(self, button: QPushButton) -> None:
        """Selecciona color.

        Args:
            button: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado del diálogo.
        """
        if button is self.bond_color_btn:
            initial = QColor(self._bond_color)
        elif button is self.atom_fill_btn:
            initial = QColor(self._atom_fill_color)
        elif button is self.atom_stroke_btn:
            initial = QColor(self._atom_stroke_color)
        else:
            initial = QColor(button.palette().button().color())
        picked = QColorDialog.getColor(initial, self, "Seleccionar color")
        if not picked.isValid():
            return
        self._set_button_color(button, picked)
        if button is self.bond_color_btn:
            self._bond_color = picked
        elif button is self.atom_fill_btn:
            self._atom_fill_color = picked
        elif button is self.atom_stroke_btn:
            self._atom_stroke_color = picked

    def _on_accept(self) -> None:
        """Maneja accept.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado del diálogo.
        """
        stroke = self.stroke_spin.value()
        bond_length = self.length_spin.value()
        self._result_style = replace(
            self._style,
            bond_length_px=bond_length,
            stroke_px=stroke,
            bond_color=self._bond_color.name(),
            atom_fill_color=self._atom_fill_color.name(),
            atom_stroke_color=self._atom_stroke_color.name(),
        )
        self._result_bond_length = bond_length
        self.accept()

    def selected_style(self) -> tuple[DrawingStyle, float]:
        """Método auxiliar para selected style.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado del diálogo.
        """
        if self._result_style is None or self._result_bond_length is None:
            return self._style, self._bond_length
        return self._result_style, self._result_bond_length


class AtomLabelDialog(QDialog):
    """Dialog for selecting atom/group label and attachment anchor."""

    def __init__(
        self,
        label: str,
        anchor: str | None,
        items: list[str],
        element_symbols: set[str],
        parent=None,
    ) -> None:
        """Inicializa el diálogo.

        Args:
            label: Descripción del parámetro.
            anchor: Descripción del parámetro.
            items: Descripción del parámetro.
            element_symbols: Descripción del parámetro.
            parent: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado del diálogo.
        """
        super().__init__(parent)
        self.setWindowTitle("Elemento")
        self._element_symbols = set(element_symbols)
        self._anchor_default = (anchor or "").strip()

        self.label_combo = QComboBox()
        self.label_combo.setEditable(True)
        self.label_combo.addItems(items)
        self.label_combo.setCurrentText(label)

        self.anchor_combo = QComboBox()
        self.anchor_combo.setEditable(True)

        form = QFormLayout()
        form.addRow("Símbolo del elemento o grupo funcional:", self.label_combo)
        form.addRow("Átomo de unión:", self.anchor_combo)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(buttons)

        self._refresh_anchor_candidates()
        editor = self.label_combo.lineEdit()
        if isinstance(editor, QLineEdit):
            editor.textChanged.connect(self._refresh_anchor_candidates)

    def _refresh_anchor_candidates(self) -> None:
        """Método auxiliar para  refresh anchor candidates.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado del diálogo.
        """
        label = self.label_combo.currentText().strip()
        candidates = self._anchor_candidates_for_label(label)
        current = self.anchor_combo.currentText().strip() or self._anchor_default
        items = ["Auto"] + candidates
        if current and current not in items and current != "Auto":
            items.append(current)
        self.anchor_combo.blockSignals(True)
        self.anchor_combo.clear()
        self.anchor_combo.addItems(items)
        if current:
            self.anchor_combo.setCurrentText(current)
        else:
            self.anchor_combo.setCurrentIndex(0)
        self.anchor_combo.blockSignals(False)

    def _anchor_candidates_for_label(self, label: str) -> list[str]:
        """Método auxiliar para  anchor candidates for label.

        Args:
            label: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado del diálogo.
        """
        cleaned = label.strip()
        if not cleaned:
            return []
        ordered: list[str] = []
        seen: set[str] = set()
        i = 0
        while i < len(cleaned):
            ch = cleaned[i]
            if not ch.isalpha():
                i += 1
                continue
            if ch.isupper():
                symbol = None
                if i + 1 < len(cleaned) and cleaned[i + 1].islower():
                    candidate = cleaned[i : i + 2]
                    if candidate in self._element_symbols:
                        symbol = candidate
                        i += 2
                if symbol is None and ch in self._element_symbols:
                    symbol = ch
                    i += 1
                if symbol is None:
                    i += 1
                    continue
                if symbol not in seen:
                    seen.add(symbol)
                    ordered.append(symbol)
                continue
            i += 1
        return ordered

    def value(self) -> tuple[str, str | None]:
        """Método auxiliar para value.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado del diálogo.
        """
        label = self.label_combo.currentText().strip()
        anchor = self.anchor_combo.currentText().strip()
        if not anchor or anchor.lower() == "auto":
            anchor = None
        return label, anchor


class QuickStartDialog(QDialog):
    """First-launch tutorial dialog."""

    def __init__(self, parent=None) -> None:
        """Inicializa el diálogo.

        Args:
            parent: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado del diálogo.
        """
        super().__init__(parent)
        self.setWindowTitle("Guía Rápida")
        self.setMinimumWidth(520)

        text = QTextBrowser()
        text.setOpenExternalLinks(False)
        text.setHtml(
            "<h3>Guía Rápida de Chemuson</h3>"
            "<p>Bienvenido a Chemuson, su editor molecular libre.</p>"
            "<ul>"
            "<li><b>Dibujar Átomos:</b> Seleccione un elemento en el panel izquierdo y haga clic en el folio.</li>"
            "<li><b>Dibujar Enlaces:</b> Haga clic en un átomo y arrastre hacia otro o hacia un espacio vacío para crear un enlace."
            " Haga clic en un enlace existente para cambiar su orden de forma incremental o use la paleta de enlaces.</li>"
            "<li><b>Dibujar Anillos:</b> Seleccione una herramienta de anillo y haga clic o arrastre desde un enlace para fusionar anillos.</li>"
            "<li><b>Borrar:</b> Use la herramienta de borrado o seleccione elementos y presione Supr.</li>"
            "<li><b>Zoom:</b> Use la rueda del ratón o los botones Zoom+/Zoom- del menú Ver o la barra de herramientas.</li>"
            "<li><b>Aromáticos:</b> Puede alternar la visualización entre 'Enlaces dobles' y 'Círculo aromático' en el menú Ver.</li>"
            "</ul>"
            "<p><i>Tip: Use 'Limpiar 2D' en el menú Estructura para organizar automáticamente sus moléculas.</i></p>"
        )

        self.no_show_checkbox = QCheckBox("No mostrar de nuevo")

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok)
        buttons.accepted.connect(self.accept)

        layout = QVBoxLayout(self)
        layout.addWidget(text)
        layout.addWidget(self.no_show_checkbox, alignment=Qt.AlignmentFlag.AlignLeft)
        layout.addWidget(buttons)
