"""
Chemuson Dialogs
Preferences and quick start dialogs.
"""
from PyQt6.QtWidgets import (
    QDialog,
    QTabWidget,
    QWidget,
    QVBoxLayout,
    QLabel,
    QCheckBox,
    QDialogButtonBox,
    QTextBrowser,
)
from PyQt6.QtCore import pyqtSignal, Qt

from core.model import ChemState


class PreferencesDialog(QDialog):
    """Application preferences with tabbed layout."""

    preferences_changed = pyqtSignal(dict)

    def __init__(self, current_state: ChemState, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Preferencias")
        self.setMinimumWidth(420)

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
        widget = QWidget()
        layout = QVBoxLayout(widget)
        label = QLabel("Opciones de color, fuentes e iconos próximamente.")
        label.setStyleSheet("color: #666666;")
        layout.addWidget(label)
        layout.addStretch()
        return widget

    def _build_rdkit_tab(self) -> QWidget:
        widget = QWidget()
        layout = QVBoxLayout(widget)
        label = QLabel("Preferencias de RDKit (longitud de enlace, limpieza 2D) próximamente.")
        label.setStyleSheet("color: #666666;")
        layout.addWidget(label)
        layout.addStretch()
        return widget

    def _on_accept(self) -> None:
        prefs = {
            "show_carbons": self.carbons_checkbox.isChecked(),
            "show_hydrogens": self.hydrogens_checkbox.isChecked(),
            "aromatic_circles": self.aromatic_checkbox.isChecked(),
        }
        self.preferences_changed.emit(prefs)
        self.accept()


class QuickStartDialog(QDialog):
    """First-launch tutorial dialog."""

    def __init__(self, parent=None) -> None:
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
