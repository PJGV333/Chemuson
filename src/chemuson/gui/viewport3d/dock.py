"""Dock discreto para el visor molecular 3D."""

from __future__ import annotations

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (
    QCheckBox,
    QButtonGroup,
    QComboBox,
    QDockWidget,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from chemuson.geometry3d import SceneMolecule3D
from chemuson.gui.theme import build_qpalette, resolve_theme_name
from chemuson.gui.theme.tokens import get_tokens

from .widget import MolecularViewportWidget


class MolecularViewportDock(QDockWidget):
    """Contenedor integrable que no expone el modelo químico al renderer."""

    liveMMToggled = pyqtSignal(bool)
    geometryConstraintRequested = pyqtSignal(str, object, float)

    def __init__(self, parent=None) -> None:
        super().__init__("Visor molecular 3D", parent)
        self._theme_name = "light"
        self.setObjectName("molecularViewportDock")
        self.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea)
        self.setMinimumWidth(300)
        container = QWidget(self)
        layout = QVBoxLayout(container)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(4)
        self.viewport = MolecularViewportWidget(container)
        self.advanced_controls_button = QPushButton("Controles de edición ▸", container)
        self.advanced_controls_button.setObjectName("molecularViewportAdvancedControlsToggle")
        self.advanced_controls_button.setCheckable(True)
        self.advanced_controls_button.setToolTip("Mostrar u ocultar herramientas geométricas")
        self.advanced_controls_panel = QWidget(container)
        self.advanced_controls_panel.setObjectName("molecularViewportAdvancedControls")
        controls = QFormLayout(self.advanced_controls_panel)
        edit_row = QHBoxLayout()
        self.edit_mode_checkbox = QCheckBox("Editar geometría 3D")
        self.live_mm_checkbox = QCheckBox("Live MM (RDKit)")
        edit_row.addWidget(self.edit_mode_checkbox)
        edit_row.addWidget(self.live_mm_checkbox)
        self.tool_combo = QComboBox()
        self.tool_combo.addItem("Mover átomo/fragmento", "move")
        self.tool_combo.addItem("Construir molécula", "build")
        self.tool_combo.addItem("Seleccionar", "select")
        self.tool_combo.addItem("Enlazar átomos", "bond")
        self.tool_combo.addItem("Medir geometría", "measure")
        self.tool_buttons: dict[str, QToolButton] = {}
        self.tool_button_group = QButtonGroup(container)
        self.tool_button_group.setExclusive(True)
        tool_row = QHBoxLayout()
        for tool, label, shortcut, description in (
            ("select", "Sel", "S", "Seleccionar átomo o enlace"),
            ("build", "Dib", "D", "Dibujar átomos y enlaces"),
            ("move", "Mov", "M", "Mover átomo o fragmento"),
            ("bond", "Enl", "B", "Crear o cambiar enlaces"),
            ("measure", "Med", "R", "Medir distancia, ángulo o dihedro"),
        ):
            button = QToolButton(container)
            button.setText(label)
            button.setCheckable(True)
            button.setToolTip(f"{description} ({shortcut})")
            button.setAccessibleName(description)
            self.tool_button_group.addButton(button)
            button.clicked.connect(
                lambda _checked=False, selected=tool: self._activate_tool_button(selected)
            )
            self.tool_buttons[tool] = button
            tool_row.addWidget(button)
        self.tool_button_group.setExclusive(False)
        self.tool_buttons["move"].setChecked(False)
        self.tool_button_group.setExclusive(True)
        self.fit_view_button = QToolButton(container)
        self.fit_view_button.setText("Fit")
        self.fit_view_button.setToolTip("Ajustar molécula a la vista (F)")
        self.fit_view_button.clicked.connect(self.viewport.fit_scene)
        tool_row.addWidget(self.fit_view_button)
        self.element_combo = QComboBox()
        self.element_combo.addItems(("C", "N", "O", "F", "P", "S", "Cl", "Br", "I"))
        self.bond_order_combo = QComboBox()
        self.bond_order_combo.addItem("Sencillo", 1)
        self.bond_order_combo.addItem("Doble", 2)
        self.bond_order_combo.addItem("Triple", 3)
        self.bond_order_combo.addItem("Aromático", 4)
        self.representation_combo = QComboBox()
        self.representation_combo.addItem("Bolas y varillas", "ball_and_stick")
        self.representation_combo.addItem("Alambre", "wireframe")
        self.representation_combo.addItem("Espacio lleno", "space_filling")
        self.live_forcefield_combo = QComboBox()
        self.live_forcefield_combo.addItem("MMFF94", "MMFF94")
        self.live_forcefield_combo.addItem("MMFF94s", "MMFF94s")
        self.live_forcefield_combo.addItem("UFF", "UFF")
        controls.addRow("Herramienta", self.tool_combo)
        controls.addRow("Elemento nuevo", self.element_combo)
        controls.addRow("Enlace nuevo", self.bond_order_combo)
        controls.addRow("Representación", self.representation_combo)
        controls.addRow("Campo Live MM", self.live_forcefield_combo)
        self.constraint_combo = QComboBox()
        self.constraint_combo.addItem("Distancia (Å)", "bond_length")
        self.constraint_combo.addItem("Ángulo (°)", "bond_angle")
        self.constraint_combo.addItem("Dihedro (°)", "dihedral")
        self.atom_ids_edit = QLineEdit()
        self.atom_ids_edit.setPlaceholderText("IDs separados por coma")
        self.target_value = QDoubleSpinBox()
        self.target_value.setDecimals(2)
        self.target_value.setRange(0.1, 5.0)
        self.target_value.setValue(1.5)
        self.apply_constraint_button = QPushButton("Aplicar geometría")
        controls.addRow("Parámetro", self.constraint_combo)
        controls.addRow("Átomos", self.atom_ids_edit)
        controls.addRow("Objetivo", self.target_value)
        controls.addRow(self.apply_constraint_button)
        self.advanced_controls_panel.setVisible(False)
        self.advanced_controls_button.toggled.connect(self._set_advanced_controls_visible)
        self.constraint_combo.currentIndexChanged.connect(self._update_target_range)
        self.apply_constraint_button.clicked.connect(self._emit_geometry_constraint)
        self.edit_mode_checkbox.toggled.connect(self.viewport.set_edit_mode)
        self.viewport.editToolChanged.connect(self._sync_tool_controls)
        self.viewport.editModeChanged.connect(self._sync_edit_mode_controls)
        self.live_mm_checkbox.toggled.connect(self.liveMMToggled.emit)
        self.tool_combo.currentIndexChanged.connect(self._set_edit_tool)
        self.element_combo.currentTextChanged.connect(self.viewport.set_build_element)
        self.representation_combo.currentIndexChanged.connect(self._set_representation)
        self.status_label = QLabel("Sin conformero 3D")
        self.status_label.setObjectName("molecularViewportStatus")
        self.status_label.setWordWrap(True)
        layout.addLayout(edit_row)
        layout.addLayout(tool_row)
        layout.addWidget(self.advanced_controls_button)
        layout.addWidget(self.advanced_controls_panel)
        layout.addWidget(self.viewport, 1)
        layout.addWidget(self.status_label)
        self.setWidget(container)
        self.viewport.statusChanged.connect(self.set_status)
        self.topLevelChanged.connect(self._on_dock_top_level_changed)
        self.dockLocationChanged.connect(self._on_dock_location_changed)
        self.set_theme(self._theme_name)
        self._sync_edit_mode_controls(self.viewport.edit_mode)

    def set_theme(self, theme_name: str) -> None:
        """Aplica paleta/QSS del dock sin afectar otros paneles o KDE."""
        self._theme_name = resolve_theme_name(theme_name)
        tokens = get_tokens(self._theme_name)
        self.setPalette(build_qpalette(self._theme_name))
        self.setAutoFillBackground(True)
        self.setStyleSheet(
            """
            QDockWidget#molecularViewportDock {
                background-color: %(bg)s;
                color: %(text1)s;
            }
            QDockWidget#molecularViewportDock::title {
                background-color: %(surface)s;
                color: %(text1)s;
                border-bottom: 1px solid %(border)s;
                padding: 5px 8px;
            }
            QDockWidget#molecularViewportDock QLabel {
                color: %(text1)s;
                background: transparent;
            }
            QDockWidget#molecularViewportDock QCheckBox {
                color: %(text1)s;
                background: transparent;
                spacing: 6px;
            }
            QDockWidget#molecularViewportDock QCheckBox:disabled {
                color: %(text3)s;
            }
            QDockWidget#molecularViewportDock QCheckBox::indicator {
                width: 16px;
                height: 16px;
                background-color: %(surface2)s;
                border: 1px solid %(borderStrong)s;
                border-radius: 3px;
            }
            QDockWidget#molecularViewportDock QCheckBox::indicator:checked {
                background-color: %(accent)s;
                border-color: %(accentStrong)s;
            }
            QDockWidget#molecularViewportDock QPushButton,
            QDockWidget#molecularViewportDock QToolButton,
            QDockWidget#molecularViewportDock QComboBox,
            QDockWidget#molecularViewportDock QLineEdit,
            QDockWidget#molecularViewportDock QDoubleSpinBox {
                background-color: %(surface2)s;
                color: %(text1)s;
                border: 1px solid %(border)s;
                border-radius: 4px;
                padding: 4px 6px;
            }
            QDockWidget#molecularViewportDock QPushButton:hover,
            QDockWidget#molecularViewportDock QToolButton:hover,
            QDockWidget#molecularViewportDock QComboBox:hover {
                background-color: %(surface3)s;
                border-color: %(borderStrong)s;
            }
            QDockWidget#molecularViewportDock QPushButton:checked,
            QDockWidget#molecularViewportDock QToolButton:checked {
                background-color: %(accentSoft)s;
                border-color: %(accentBorder)s;
                color: %(accentStrong)s;
            }
            QDockWidget#molecularViewportDock QPushButton:disabled,
            QDockWidget#molecularViewportDock QToolButton:disabled,
            QDockWidget#molecularViewportDock QComboBox:disabled,
            QDockWidget#molecularViewportDock QLineEdit:disabled,
            QDockWidget#molecularViewportDock QDoubleSpinBox:disabled {
                background-color: %(surface3)s;
                color: %(text3)s;
                border-color: %(border)s;
            }
            QDockWidget#molecularViewportDock QComboBox QAbstractItemView {
                background-color: %(surface)s;
                color: %(text1)s;
                selection-background-color: %(accentSoft)s;
                selection-color: %(accentStrong)s;
            }
            QDockWidget#molecularViewportDock QToolTip {
                background-color: %(surface)s;
                color: %(text1)s;
                border: 1px solid %(border)s;
            }
            """
            % {key: str(value) for key, value in tokens.items() if isinstance(value, str)}
        )

    def _on_dock_top_level_changed(self, _floating: bool) -> None:
        self.set_theme(self._theme_name)

    def _on_dock_location_changed(self, _area: Qt.DockWidgetArea) -> None:
        self.set_theme(self._theme_name)

    def _clear_active_tool_button(self) -> None:
        self.tool_button_group.setExclusive(False)
        for button in self.tool_buttons.values():
            button.setChecked(False)
        self.tool_button_group.setExclusive(True)

    def _sync_edit_mode_controls(self, enabled: bool) -> None:
        self.edit_mode_checkbox.blockSignals(True)
        self.edit_mode_checkbox.setChecked(bool(enabled))
        self.edit_mode_checkbox.blockSignals(False)
        if enabled:
            self._sync_tool_controls(self.viewport.edit_tool)
        else:
            self._clear_active_tool_button()

    def _set_advanced_controls_visible(self, visible: bool) -> None:
        self.advanced_controls_panel.setVisible(visible)
        label = "Controles de edición ▾" if visible else "Controles de edición ▸"
        self.advanced_controls_button.setText(label)

    def _set_representation(self, _index: int) -> None:
        self.viewport.set_representation(str(self.representation_combo.currentData()))

    def _set_edit_tool(self, _index: int) -> None:
        tool = str(self.tool_combo.currentData())
        self.viewport.set_edit_tool(tool)
        if not self.edit_mode_checkbox.isChecked():
            self.edit_mode_checkbox.setChecked(True)

    def _activate_tool_button(self, tool: str) -> None:
        index = self.tool_combo.findData(tool)
        if index >= 0:
            self.tool_combo.setCurrentIndex(index)
        if not self.edit_mode_checkbox.isChecked():
            self.edit_mode_checkbox.setChecked(True)

    def _sync_tool_controls(self, tool: str) -> None:
        if not self.viewport.edit_mode:
            self._clear_active_tool_button()
            return
        index = self.tool_combo.findData(tool)
        if index >= 0 and self.tool_combo.currentIndex() != index:
            self.tool_combo.setCurrentIndex(index)
        button = self.tool_buttons.get(tool)
        if button is not None:
            button.setChecked(True)

    def _update_target_range(self, _index: int) -> None:
        operation = self.constraint_combo.currentData()
        if operation == "bond_length":
            self.target_value.setDecimals(2)
            self.target_value.setRange(0.1, 5.0)
            self.target_value.setSingleStep(0.1)
            self.target_value.setValue(1.5)
        elif operation == "bond_angle":
            self.target_value.setDecimals(1)
            self.target_value.setRange(1.0, 179.0)
            self.target_value.setSingleStep(1.0)
            self.target_value.setValue(109.5)
        else:
            self.target_value.setDecimals(1)
            self.target_value.setRange(-180.0, 180.0)
            self.target_value.setSingleStep(5.0)
            self.target_value.setValue(180.0)

    def _emit_geometry_constraint(self) -> None:
        try:
            atom_ids = tuple(
                int(value.strip())
                for value in self.atom_ids_edit.text().split(",")
                if value.strip()
            )
        except ValueError:
            self.set_status("Los IDs de átomos deben ser enteros separados por coma")
            return
        operation = str(self.constraint_combo.currentData())
        expected = {"bond_length": 2, "bond_angle": 3, "dihedral": 4}[operation]
        if len(atom_ids) != expected:
            self.set_status(f"{operation}: se requieren {expected} IDs de átomos")
            return
        self.geometryConstraintRequested.emit(
            operation, atom_ids, float(self.target_value.value())
        )

    def selected_live_forcefield(self) -> str:
        return str(self.live_forcefield_combo.currentData())

    def selected_bond_order(self) -> tuple[int, bool]:
        order = int(self.bond_order_combo.currentData() or 1)
        return (1 if order == 4 else order, order == 4)

    def set_scene(
        self,
        scene: SceneMolecule3D | None,
        *,
        preserve_interaction: bool = False,
    ) -> None:
        self.viewport.set_scene(scene, preserve_interaction=preserve_interaction)
        if scene is None or not scene.atoms:
            self.set_status("Sin conformero 3D")
        else:
            self.set_status(f"{len(scene.atoms)} átomos · {len(scene.bonds)} enlaces")

    def clear_scene(self) -> None:
        self.set_scene(None)

    def set_status(self, message: str) -> None:
        self.status_label.setText(str(message))
