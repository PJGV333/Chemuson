"""Dock discreto para el visor molecular 3D."""

from __future__ import annotations

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDockWidget,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from chemuson.geometry3d import SceneMolecule3D

from .widget import MolecularViewportWidget


class MolecularViewportDock(QDockWidget):
    """Contenedor integrable que no expone el modelo químico al renderer."""

    liveMMToggled = pyqtSignal(bool)
    geometryConstraintRequested = pyqtSignal(str, object, float)

    def __init__(self, parent=None) -> None:
        super().__init__("Visor molecular 3D", parent)
        self.setObjectName("molecularViewportDock")
        self.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea)
        self.setMinimumWidth(360)
        container = QWidget(self)
        layout = QVBoxLayout(container)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(4)
        self.viewport = MolecularViewportWidget(container)
        controls = QFormLayout()
        edit_row = QHBoxLayout()
        self.edit_mode_checkbox = QCheckBox("Editar geometría 3D")
        self.live_mm_checkbox = QCheckBox("Live MM (RDKit)")
        edit_row.addWidget(self.edit_mode_checkbox)
        edit_row.addWidget(self.live_mm_checkbox)
        self.tool_combo = QComboBox()
        self.tool_combo.addItem("Mover átomo/fragmento", "move")
        self.tool_combo.addItem("Construir molécula", "build")
        self.element_combo = QComboBox()
        self.element_combo.addItems(("C", "N", "O", "F", "P", "S", "Cl", "Br", "I"))
        self.bond_order_combo = QComboBox()
        self.bond_order_combo.addItem("Sencillo", 1)
        self.bond_order_combo.addItem("Doble", 2)
        self.bond_order_combo.addItem("Triple", 3)
        self.bond_order_combo.addItem("Aromático", 4)
        self.live_forcefield_combo = QComboBox()
        self.live_forcefield_combo.addItem("MMFF94", "MMFF94")
        self.live_forcefield_combo.addItem("MMFF94s", "MMFF94s")
        self.live_forcefield_combo.addItem("UFF", "UFF")
        controls.addRow("Herramienta", self.tool_combo)
        controls.addRow("Elemento nuevo", self.element_combo)
        controls.addRow("Enlace nuevo", self.bond_order_combo)
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
        self.constraint_combo.currentIndexChanged.connect(self._update_target_range)
        self.apply_constraint_button.clicked.connect(self._emit_geometry_constraint)
        self.edit_mode_checkbox.toggled.connect(self.viewport.set_edit_mode)
        self.live_mm_checkbox.toggled.connect(self.liveMMToggled.emit)
        self.tool_combo.currentIndexChanged.connect(self._set_edit_tool)
        self.element_combo.currentTextChanged.connect(self.viewport.set_build_element)
        self.status_label = QLabel("Sin conformero 3D")
        self.status_label.setObjectName("molecularViewportStatus")
        self.status_label.setWordWrap(True)
        layout.addLayout(edit_row)
        layout.addLayout(controls)
        layout.addWidget(self.viewport, 1)
        layout.addWidget(self.status_label)
        self.setWidget(container)
        self.viewport.statusChanged.connect(self.set_status)

    def _set_edit_tool(self, _index: int) -> None:
        tool = str(self.tool_combo.currentData())
        self.viewport.set_edit_tool(tool)
        if tool == "build" and not self.edit_mode_checkbox.isChecked():
            self.edit_mode_checkbox.setChecked(True)

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
