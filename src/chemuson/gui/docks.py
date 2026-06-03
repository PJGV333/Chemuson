"""
Docks de Chemuson.

Paneles laterales para plantillas, inspección de propiedades y apariencia.
"""
import json
from pathlib import Path

from PyQt6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QDockWidget,
    QFileDialog,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QFormLayout,
    QComboBox,
    QTreeWidget,
    QTreeWidgetItem,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QPushButton,
    QMenu,
    QTabWidget,
)
from PyQt6.QtCore import Qt, pyqtSignal
from chemuson.gui.style import DrawingStyle


class PlantillasDock(QDockWidget):
    """
    Dock que muestra listas de plantillas químicas.
    """
    template_selected = pyqtSignal(dict)
    new_category_requested = pyqtSignal()
    import_requested = pyqtSignal()
    export_requested = pyqtSignal()
    rename_category_requested = pyqtSignal(str)
    delete_category_requested = pyqtSignal(str)
    rename_template_requested = pyqtSignal(str)
    delete_template_requested = pyqtSignal(str)

    def __init__(self, parent=None):
        """Inicializa el dock de plantillas."""
        super().__init__("Plantillas", parent)
        self.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea)
        
        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)

        controls = QHBoxLayout()
        controls.setContentsMargins(6, 6, 6, 2)
        self.btn_new_category = QPushButton("Nueva categoría")
        self.btn_import = QPushButton("Importar")
        self.btn_export = QPushButton("Exportar")
        self.btn_new_category.clicked.connect(self.new_category_requested.emit)
        self.btn_import.clicked.connect(self.import_requested.emit)
        self.btn_export.clicked.connect(self.export_requested.emit)
        controls.addWidget(self.btn_new_category)
        controls.addWidget(self.btn_import)
        controls.addWidget(self.btn_export)
        layout.addLayout(controls)
        
        self.tree = QTreeWidget()
        self.tree.setHeaderHidden(True)
        self.tree.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.tree.customContextMenuRequested.connect(self._show_context_menu)
        self.tree.itemActivated.connect(self._emit_template)
        self.tree.itemDoubleClicked.connect(self._emit_template)
        layout.addWidget(self.tree)
        
        self.setWidget(container)
        self.set_templates([])

    def set_templates(self, grouped_templates: list[dict]) -> None:
        """Carga categorías y plantillas en el árbol del dock."""
        self.tree.clear()
        if not grouped_templates:
            placeholder = QTreeWidgetItem(["(Sin plantillas)"])
            placeholder.setFlags(placeholder.flags() & ~Qt.ItemFlag.ItemIsSelectable)
            self.tree.addTopLevelItem(placeholder)
            return
        for group in grouped_templates:
            category = str(group.get("name", "")).strip()
            if not category:
                continue
            group_item = QTreeWidgetItem([category])
            group_item.setData(
                0,
                Qt.ItemDataRole.UserRole,
                {"kind": "category", "name": category},
            )
            for template in group.get("templates", []):
                name = str(template.get("name", "")).strip()
                template_id = str(template.get("id", "")).strip()
                if not name or not template_id:
                    continue
                child = QTreeWidgetItem([name])
                icon = template.get("icon")
                if icon is not None:
                    child.setIcon(0, icon)
                child.setData(
                    0,
                    Qt.ItemDataRole.UserRole,
                    {
                        "kind": "template",
                        "id": template_id,
                        "name": name,
                        "category": category,
                    },
                )
                group_item.addChild(child)
            self.tree.addTopLevelItem(group_item)
            group_item.setExpanded(True)

    def _emit_template(self, item: QTreeWidgetItem, _column: int = 0) -> None:
        """Emite la plantilla seleccionada si el item tiene datos."""
        payload = item.data(0, Qt.ItemDataRole.UserRole)
        if isinstance(payload, dict) and payload.get("kind") == "template":
            self.template_selected.emit(payload)

    def _show_context_menu(self, pos) -> None:
        """Muestra menú contextual según el tipo de item."""
        item = self.tree.itemAt(pos)
        if item is None:
            return
        payload = item.data(0, Qt.ItemDataRole.UserRole)
        if not isinstance(payload, dict):
            return
        menu = QMenu(self.tree)
        kind = payload.get("kind")
        if kind == "template":
            insert_action = menu.addAction("Insertar")
            rename_action = menu.addAction("Renombrar plantilla...")
            delete_action = menu.addAction("Eliminar plantilla")
            chosen = menu.exec(self.tree.viewport().mapToGlobal(pos))
            if chosen == insert_action:
                self.template_selected.emit(payload)
            elif chosen == rename_action:
                self.rename_template_requested.emit(str(payload.get("id", "")))
            elif chosen == delete_action:
                self.delete_template_requested.emit(str(payload.get("id", "")))
            return
        if kind == "category":
            rename_action = menu.addAction("Renombrar categoría...")
            delete_action = menu.addAction("Eliminar categoría")
            chosen = menu.exec(self.tree.viewport().mapToGlobal(pos))
            if chosen == rename_action:
                self.rename_category_requested.emit(str(payload.get("name", "")))
            elif chosen == delete_action:
                self.delete_category_requested.emit(str(payload.get("name", "")))


class InspectorDock(QDockWidget):
    """
    Dock que muestra propiedades del objeto seleccionado.
    """
    property_changed = pyqtSignal(str, object)

    def __init__(self, parent=None):
        """Inicializa el dock de inspección."""
        super().__init__("Inspector", parent)
        self.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea)
        
        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)
        
        self.info_label = QLabel("Nada seleccionado")
        self.info_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.info_label.setStyleSheet("color: #666666; font-style: italic; padding: 10px;")
        layout.addWidget(self.info_label)
        
        self.prop_table = QTableWidget(0, 2)
        self.prop_table.setHorizontalHeaderLabels(["Propiedad", "Valor"])
        self.prop_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.prop_table.verticalHeader().setVisible(False)
        self.prop_table.setAlternatingRowColors(True)
        self.prop_table.setVisible(False)
        layout.addWidget(self.prop_table)
        
        layout.addStretch()
        self.setWidget(container)

    def set_atom(self, atom) -> None:
        """Carga detalles de un átomo en el inspector."""
        self.update_selection(1, 0, {
            "element": atom.element,
            "id": atom.id,
            "charge": atom.charge,
            "x": atom.x,
            "y": atom.y,
        })

    def set_bond(self, bond) -> None:
        """Carga detalles de un enlace en el inspector."""
        self.update_selection(0, 1, 0, {
            "order": bond.order,
            "style": bond.style,
            "aromatic": bond.is_aromatic,
        })

    def update_selection(self, num_atoms: int, num_bonds: int, num_text: int, details: dict):
        """Actualiza el inspector con información de selección."""
        if (
            num_atoms == 0
            and num_bonds == 0
            and num_text == 0
            and details.get("type") != "energy_diagram"
        ):
            self.info_label.setText("Nada seleccionado")
            self.info_label.setVisible(True)
            self.prop_table.setVisible(False)
            return
            
        self.info_label.setVisible(False)
        self.prop_table.setVisible(True)
        self.prop_table.setRowCount(0)
        
        data = []
        if details.get("type") == "energy_diagram":
            data = [
                ("Tipo", "Diagrama de energía"),
                ("Preset", str(details.get("kind", "?"))),
                ("Etiqueta", str(details.get("label", ""))),
                ("Cajas", str(details.get("boxes", "?"))),
                ("Ocupación", str(details.get("occupancies", ""))),
            ]
        elif num_atoms == 1 and num_bonds == 0 and num_text == 0:
            data = [
                ("Tipo", "Átomo"),
                ("Elemento", details.get("element", "?")),
                ("ID", str(details.get("id", "?"))),
                ("Carga", str(details.get("charge", 0))),
                ("X", f"{details.get('x', 0):.1f}"),
                ("Y", f"{details.get('y', 0):.1f}"),
            ]
        elif num_atoms == 0 and num_bonds == 1 and num_text == 0:
            data = [
                ("Tipo", "Enlace"),
                ("Orden", str(details.get("order", 1))),
                ("Estilo", str(details.get("style", "Plain"))),
                ("Aromático", "Sí" if details.get("aromatic") else "No"),
            ]
        elif num_text == 1 and num_atoms == 0 and num_bonds == 0:
            font = details.get("font")
            data = [
                ("Tipo", "Texto"),
                ("Fuente", font.family() if font else "?"),
                ("Tamaño", str(font.pointSize()) if font else "?"),
                ("Sub/Sup", "Sub" if details.get("sub") else ("Sup" if details.get("sup") else "Normal")),
            ]
        else:
            data = [
                ("Selección", "Múltiple"),
                ("Átomos", str(num_atoms)),
                ("Enlaces", str(num_bonds)),
                ("Texto", str(num_text)),
            ]
            
        self.prop_table.setRowCount(len(data))
        for i, (key, val) in enumerate(data):
            self.prop_table.setItem(i, 0, QTableWidgetItem(key))
            self.prop_table.setItem(i, 1, QTableWidgetItem(val))


class ChemicalPropertiesDock(QDockWidget):
    """Dock de propiedades químicas calculadas para el documento activo."""

    def __init__(self, parent=None):
        """Inicializa el dock de inteligencia química básica."""
        super().__init__("Propiedades químicas", parent)
        self.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea)

        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)

        self.info_label = QLabel("Sin estructura")
        self.info_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.info_label.setStyleSheet("color: #666666; font-style: italic; padding: 10px;")
        layout.addWidget(self.info_label)

        self.prop_table = QTableWidget(0, 2)
        self.prop_table.setHorizontalHeaderLabels(["Propiedad", "Valor"])
        self.prop_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.prop_table.verticalHeader().setVisible(False)
        self.prop_table.setAlternatingRowColors(True)
        layout.addWidget(self.prop_table)

        controls = QHBoxLayout()
        controls.setContentsMargins(6, 4, 6, 6)
        self.btn_copy_properties = QPushButton("Copiar propiedades")
        self.btn_copy_properties.clicked.connect(self.copy_properties_tsv)
        controls.addStretch()
        controls.addWidget(self.btn_copy_properties)
        layout.addLayout(controls)

        self.setWidget(container)
        self._rows: list[tuple[str, str]] = []

    def update_properties(self, rows: list[tuple[str, str]]) -> None:
        """Actualiza la tabla de propiedades calculadas."""
        self._rows = [(str(key), str(value)) for key, value in (rows or [])]
        self.prop_table.setRowCount(0)
        if not self._rows:
            self.info_label.setVisible(True)
            self.prop_table.setVisible(False)
            self.btn_copy_properties.setVisible(False)
            return

        self.info_label.setVisible(False)
        self.prop_table.setVisible(True)
        self.btn_copy_properties.setVisible(True)
        self.prop_table.setRowCount(len(self._rows))
        for row, (key, value) in enumerate(self._rows):
            self.prop_table.setItem(row, 0, QTableWidgetItem(str(key)))
            self.prop_table.setItem(row, 1, QTableWidgetItem(str(value)))

    def copy_properties_tsv(self) -> str:
        """Copia las propiedades visibles como TSV y devuelve el texto."""
        if not self._rows:
            return ""
        lines = ["Propiedad\tValor"]
        lines.extend(f"{key}\t{value}" for key, value in self._rows)
        text = "\n".join(lines)
        QApplication.clipboard().setText(text)
        return text


class SpectroscopyDock(QDockWidget):
    """Dock con predicción espectral y selección cruzada de átomos."""

    peak_atom_selected = pyqtSignal(int)

    def __init__(self, parent=None):
        super().__init__("Espectros", parent)
        self.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea)

        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)

        self.info_label = QLabel("Sin predicción espectral")
        self.info_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.info_label.setStyleSheet("color: #666666; font-style: italic; padding: 10px;")
        layout.addWidget(self.info_label)

        self.tabs = QTabWidget()
        self.proton_table = self._make_table(["δ ppm", "Int.", "Entorno", "Conf.", "Átomo"])
        self.carbon_table = self._make_table(["δ ppm", "Entorno", "Conf.", "Átomo"])
        self.mass_table = self._make_table(["m/z", "Int. %", "Pico", "Conf."])
        self.tabs.addTab(self.proton_table, "¹H NMR")
        self.tabs.addTab(self.carbon_table, "¹³C NMR")
        self.tabs.addTab(self.mass_table, "MS")
        layout.addWidget(self.tabs)

        controls = QHBoxLayout()
        controls.setContentsMargins(6, 4, 6, 6)
        self.btn_copy_table = QPushButton("Copiar tabla")
        self.btn_copy_table.clicked.connect(self.copy_current_table_tsv)
        self.btn_export_table = QPushButton("Exportar tabla...")
        self.btn_export_table.clicked.connect(self.export_current_table_tsv)
        controls.addStretch()
        controls.addWidget(self.btn_copy_table)
        controls.addWidget(self.btn_export_table)
        layout.addLayout(controls)

        self.setWidget(container)
        self.update_prediction(None)

    def update_prediction(self, prediction) -> None:
        """Carga una predicción espectral en las tablas del dock."""
        if prediction is not None:
            source = str(getattr(prediction, "source", "heuristic-v1"))
            confidence = float(getattr(prediction, "confidence", 0.0) or 0.0)
            self.info_label.setText(
                f"Predicción estimada ({source}); confianza global {confidence:.2f}. No sustituye datos experimentales."
            )
        else:
            self.info_label.setText("Sin predicción espectral")
        self._fill_proton(getattr(prediction, "proton_nmr", []) if prediction is not None else [])
        self._fill_carbon(getattr(prediction, "carbon_nmr", []) if prediction is not None else [])
        self._fill_mass(getattr(prediction, "mass_spectrum", []) if prediction is not None else [])
        has_rows = any(
            table.rowCount() > 0
            for table in (self.proton_table, self.carbon_table, self.mass_table)
        )
        self.info_label.setVisible(True)
        self.tabs.setVisible(has_rows)
        self.btn_copy_table.setVisible(has_rows)
        self.btn_export_table.setVisible(has_rows)

    def _make_table(self, headers: list[str]) -> QTableWidget:
        table = QTableWidget(0, len(headers))
        table.setHorizontalHeaderLabels(headers)
        table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        table.verticalHeader().setVisible(False)
        table.setAlternatingRowColors(True)
        table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        table.setSelectionMode(QTableWidget.SelectionMode.SingleSelection)
        table.itemSelectionChanged.connect(lambda table=table: self._emit_selected_atom(table))
        return table

    def _fill_proton(self, peaks) -> None:
        rows = [
            (
                f"{float(peak.shift_ppm):.2f}",
                f"{int(peak.hydrogens)}H",
                str(peak.environment),
                f"{float(getattr(peak, 'confidence', 0.0)):.2f}",
                int(peak.atom_id),
            )
            for peak in peaks
        ]
        self._fill_table(self.proton_table, rows, atom_column=4)

    def _fill_carbon(self, peaks) -> None:
        rows = [
            (
                f"{float(peak.shift_ppm):.1f}",
                str(peak.environment),
                f"{float(getattr(peak, 'confidence', 0.0)):.2f}",
                int(peak.atom_id),
            )
            for peak in peaks
        ]
        self._fill_table(self.carbon_table, rows, atom_column=3)

    def _fill_mass(self, peaks) -> None:
        rows = [
            (
                f"{float(peak.mz):.4f}",
                f"{float(peak.intensity):.1f}",
                str(peak.label),
                f"{float(getattr(peak, 'confidence', 0.0)):.2f}",
            )
            for peak in peaks
        ]
        self._fill_table(self.mass_table, rows, atom_column=None)

    def copy_current_table_tsv(self) -> str:
        """Copia la tabla espectral activa como TSV y devuelve el texto."""
        text = self.current_table_tsv()
        if text:
            QApplication.clipboard().setText(text)
        return text

    def current_table_tsv(self) -> str:
        """Devuelve la tabla espectral activa como TSV."""
        table = self.tabs.currentWidget()
        if not isinstance(table, QTableWidget):
            return ""
        headers = [table.horizontalHeaderItem(col).text() for col in range(table.columnCount())]
        lines = ["\t".join(headers)]
        for row in range(table.rowCount()):
            values = []
            for col in range(table.columnCount()):
                item = table.item(row, col)
                values.append(item.text() if item is not None else "")
            lines.append("\t".join(values))
        return "\n".join(lines)

    def export_current_table_tsv(self, path: str | None = None) -> str:
        """Exporta la tabla espectral activa como TSV y devuelve la ruta."""
        text = self.current_table_tsv()
        if not text:
            return ""
        if path is None:
            path, _ = QFileDialog.getSaveFileName(
                self,
                "Exportar tabla espectral",
                "",
                "Tabla TSV (*.tsv);;Texto (*.txt)",
            )
        if not path:
            return ""
        Path(path).write_text(text + "\n", encoding="utf-8")
        return str(path)

    def _fill_table(self, table: QTableWidget, rows: list[tuple], atom_column: int | None) -> None:
        table.blockSignals(True)
        table.setRowCount(len(rows))
        for row, values in enumerate(rows):
            atom_id = None
            if atom_column is not None:
                atom_id = int(values[atom_column])
            for col, value in enumerate(values):
                item = QTableWidgetItem(str(value))
                if atom_id is not None:
                    item.setData(Qt.ItemDataRole.UserRole, atom_id)
                table.setItem(row, col, item)
        table.blockSignals(False)

    def _emit_selected_atom(self, table: QTableWidget) -> None:
        selected = table.selectedItems()
        if not selected:
            return
        atom_id = selected[0].data(Qt.ItemDataRole.UserRole)
        if atom_id is None:
            return
        self.peak_atom_selected.emit(int(atom_id))


class CompChemDock(QDockWidget):
    """Dock MVP para generación 3D y química computacional."""

    generate_requested = pyqtSignal()
    optimize_requested = pyqtSignal()
    project_requested = pyqtSignal()
    reset_requested = pyqtSignal()
    export_xyz_requested = pyqtSignal()
    export_input_requested = pyqtSignal(str)

    def __init__(self, parent=None):
        super().__init__("3D / CompChem", parent)
        self.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea)

        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(6, 6, 6, 6)

        form = QFormLayout()
        self.backend_combo = QComboBox()
        self.backend_combo.addItem("RDKit aislado", "rdkit")
        self.backend_combo.addItem("Open Babel", "openbabel")
        self.forcefield_combo = QComboBox()
        for label, value in (("UFF", "UFF"), ("MMFF94", "MMFF94"), ("MMFF94s", "MMFF94s")):
            self.forcefield_combo.addItem(label, value)
        form.addRow("Backend", self.backend_combo)
        form.addRow("Campo de fuerza", self.forcefield_combo)
        layout.addLayout(form)

        self.status_label = QLabel("Sin conformero 3D")
        self.status_label.setWordWrap(True)
        self.status_label.setStyleSheet("color: #555555; padding: 4px;")
        layout.addWidget(self.status_label)

        self.frame_table = QTableWidget(0, 4)
        self.frame_table.setHorizontalHeaderLabels(["Paso", "Energía", "Conv.", "Mensaje"])
        self.frame_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.frame_table.verticalHeader().setVisible(False)
        self.frame_table.setAlternatingRowColors(True)
        layout.addWidget(self.frame_table)

        row1 = QHBoxLayout()
        self.btn_generate = QPushButton("Generar conformero 3D")
        self.btn_optimize = QPushButton("Optimizar")
        self.btn_reset = QPushButton("Reset/regenerar")
        self.btn_generate.clicked.connect(self.generate_requested.emit)
        self.btn_optimize.clicked.connect(self.optimize_requested.emit)
        self.btn_reset.clicked.connect(self.reset_requested.emit)
        row1.addWidget(self.btn_generate)
        row1.addWidget(self.btn_optimize)
        row1.addWidget(self.btn_reset)
        layout.addLayout(row1)

        row2 = QHBoxLayout()
        self.btn_project = QPushButton("Proyectar a 2D")
        self.btn_export_xyz = QPushButton("Exportar XYZ")
        self.btn_export_input = QPushButton("Exportar input")
        export_menu = QMenu(self.btn_export_input)
        for label, value in (("ORCA", "orca"), ("Gaussian", "gaussian"), ("NWChem", "nwchem")):
            action = export_menu.addAction(label)
            action.triggered.connect(lambda _checked=False, value=value: self.export_input_requested.emit(value))
        self.btn_export_input.setMenu(export_menu)
        self.btn_project.clicked.connect(self.project_requested.emit)
        self.btn_export_xyz.clicked.connect(self.export_xyz_requested.emit)
        row2.addWidget(self.btn_project)
        row2.addWidget(self.btn_export_xyz)
        row2.addWidget(self.btn_export_input)
        layout.addLayout(row2)

        self.setWidget(container)
        self.set_busy(False)
        self.set_has_coordinates(False)

    def selected_backend(self) -> str:
        return str(self.backend_combo.currentData() or "rdkit")

    def selected_forcefield(self) -> str:
        return str(self.forcefield_combo.currentData() or "UFF")

    def set_busy(self, busy: bool) -> None:
        for widget in (
            self.backend_combo,
            self.forcefield_combo,
            self.btn_generate,
            self.btn_optimize,
            self.btn_reset,
            self.btn_project,
            self.btn_export_xyz,
            self.btn_export_input,
        ):
            widget.setEnabled(not busy)

    def set_has_coordinates(self, has_coordinates: bool) -> None:
        self.btn_project.setEnabled(bool(has_coordinates))
        self.btn_export_xyz.setEnabled(bool(has_coordinates))
        self.btn_export_input.setEnabled(bool(has_coordinates))

    def set_status(self, text: str) -> None:
        self.status_label.setText(str(text or ""))

    def clear_frames(self) -> None:
        self.frame_table.setRowCount(0)

    def add_frame(self, frame) -> None:
        row = self.frame_table.rowCount()
        self.frame_table.insertRow(row)
        energy = getattr(frame, "energy", None)
        values = [
            str(getattr(frame, "step", row + 1)),
            "N/D" if energy is None else f"{float(energy):.6f}",
            "sí" if bool(getattr(frame, "converged", False)) else "no",
            str(getattr(frame, "message", "") or ""),
        ]
        for col, value in enumerate(values):
            self.frame_table.setItem(row, col, QTableWidgetItem(value))

    def set_result_summary(self, result, *, backend: str, cache_state: str = "miss") -> None:
        energy = getattr(result, "energy", None)
        converged = bool(getattr(result, "converged", False))
        method = str(getattr(result, "method", "") or "")
        message = str(getattr(result, "message", "") or "")
        parts = [
            f"Backend: {backend}",
            f"Método/campo: {method or self.selected_forcefield()}",
            f"Convergencia: {'sí' if converged else 'no'}",
            f"Cache: {cache_state}",
        ]
        if energy is not None:
            parts.append(f"Energía: {float(energy):.6f}")
        if message:
            parts.append(f"Mensaje: {message}")
        self.set_status(" | ".join(parts))


class ValidationDock(QDockWidget):
    """Dock de diagnóstico interactivo de valencias."""

    issue_selected = pyqtSignal(int)
    correction_requested = pyqtSignal(int, str)
    refresh_requested = pyqtSignal()
    next_requested = pyqtSignal()
    previous_requested = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__("Validación", parent)
        self.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea)

        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)

        self.info_label = QLabel("Sin validación ejecutada")
        self.info_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.info_label.setStyleSheet("color: #666666; font-style: italic; padding: 10px;")
        layout.addWidget(self.info_label)

        self.issue_table = QTableWidget(0, 6)
        self.issue_table.setHorizontalHeaderLabels(
            ["Afectado", "Sev.", "Código", "Valencia", "Mensaje", "Sugerencia"]
        )
        self.issue_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.issue_table.verticalHeader().setVisible(False)
        self.issue_table.setAlternatingRowColors(True)
        self.issue_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.issue_table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.issue_table.itemSelectionChanged.connect(self._emit_selected_issue)
        layout.addWidget(self.issue_table)

        controls = QHBoxLayout()
        controls.setContentsMargins(6, 4, 6, 6)
        self.btn_refresh = QPushButton("Validar")
        self.btn_previous = QPushButton("Anterior")
        self.btn_next = QPushButton("Siguiente")
        self.btn_copy_report = QPushButton("Copiar reporte")
        self.btn_export_report = QPushButton("Exportar reporte...")
        self.action_combo = QComboBox()
        self.btn_apply_correction = QPushButton("Aplicar corrección")
        self.btn_refresh.clicked.connect(self.refresh_requested.emit)
        self.btn_previous.clicked.connect(self.previous_requested.emit)
        self.btn_next.clicked.connect(self.next_requested.emit)
        self.btn_copy_report.clicked.connect(self.copy_report_tsv)
        self.btn_export_report.clicked.connect(self.export_report)
        self.btn_apply_correction.clicked.connect(self._emit_correction_requested)
        controls.addWidget(self.btn_refresh)
        controls.addWidget(self.btn_previous)
        controls.addWidget(self.btn_next)
        controls.addWidget(self.btn_copy_report)
        controls.addWidget(self.btn_export_report)
        controls.addWidget(self.action_combo)
        controls.addWidget(self.btn_apply_correction)
        layout.addLayout(controls)

        self.setWidget(container)
        self._issues: dict[int, object] = {}
        self.set_issues({})

    def set_issues(self, issues: dict[int, object]) -> None:
        """Carga diagnósticos de validación en la tabla."""
        self._issues = dict(issues or {})
        rows = []
        for atom_id, issue in sorted(self._issues.items()):
            element = str(getattr(issue, "element", "") or "")
            affected = f"Átomo {int(atom_id)}" + (f" ({element})" if element else "")
            valence_calc = (
                f"{float(getattr(issue, 'observed_valence', 0.0) or 0.0):.2f} = "
                f"enlaces {float(getattr(issue, 'bond_order_sum', 0.0) or 0.0):.2f} + "
                f"H asignados {int(getattr(issue, 'assigned_h', 0) or 0)} + "
                f"H implícitos {int(getattr(issue, 'implicit_h', 0) or 0)}"
            )
            rows.append(
                (
                    int(atom_id),
                    affected,
                    str(getattr(issue, "severity", "error")),
                    str(getattr(issue, "code", "")),
                    valence_calc,
                    str(getattr(issue, "message", "")),
                    str(getattr(issue, "hint", "") or getattr(issue, "suggestion", "") or ""),
                )
            )

        self.issue_table.blockSignals(True)
        self.issue_table.setRowCount(len(rows))
        for row, values in enumerate(rows):
            atom_id = int(values[0])
            for col, value in enumerate(values[1:]):
                item = QTableWidgetItem(str(value))
                item.setData(Qt.ItemDataRole.UserRole, atom_id)
                self.issue_table.setItem(row, col, item)
        self.issue_table.blockSignals(False)

        if rows:
            self.info_label.setText(f"{len(rows)} error(es) de valencia")
            self.info_label.setVisible(True)
            self.issue_table.setVisible(True)
            self.btn_copy_report.setVisible(True)
            self.btn_export_report.setVisible(True)
        else:
            self.info_label.setText("Sin errores de valencia")
            self.info_label.setVisible(True)
            self.issue_table.setVisible(False)
            self.btn_copy_report.setVisible(False)
            self.btn_export_report.setVisible(False)
        self._update_action_controls()

    def select_atom(self, atom_id: int) -> None:
        """Selecciona visualmente el diagnóstico asociado a un átomo."""
        target = int(atom_id)
        self.issue_table.blockSignals(True)
        try:
            self.issue_table.clearSelection()
            for row in range(self.issue_table.rowCount()):
                item = self.issue_table.item(row, 0)
                if item is not None and int(item.data(Qt.ItemDataRole.UserRole)) == target:
                    self.issue_table.selectRow(row)
                    self.issue_table.scrollToItem(item)
                    break
        finally:
            self.issue_table.blockSignals(False)
        self._update_action_controls(target)

    def _emit_selected_issue(self) -> None:
        selected = self.issue_table.selectedItems()
        if not selected:
            self._update_action_controls()
            return
        atom_id = selected[0].data(Qt.ItemDataRole.UserRole)
        if atom_id is None:
            self._update_action_controls()
            return
        self._update_action_controls(int(atom_id))
        self.issue_selected.emit(int(atom_id))

    def _selected_atom_id(self) -> int | None:
        selected = self.issue_table.selectedItems()
        if not selected:
            return None
        atom_id = selected[0].data(Qt.ItemDataRole.UserRole)
        return None if atom_id is None else int(atom_id)

    def _issue_actions(self, atom_id: int) -> tuple[dict[str, str], ...]:
        issue = self._issues.get(int(atom_id))
        if issue is None:
            return ()
        suggested = getattr(issue, "suggested_actions", None)
        if callable(suggested):
            return tuple(suggested())
        actions = getattr(issue, "correction_actions", ()) or ()
        return tuple(action.as_dict() for action in actions if hasattr(action, "as_dict"))

    def _update_action_controls(self, atom_id: int | None = None) -> None:
        if atom_id is None:
            atom_id = self._selected_atom_id()
        self.action_combo.blockSignals(True)
        self.action_combo.clear()
        if atom_id is not None:
            for action in self._issue_actions(int(atom_id)):
                self.action_combo.addItem(str(action.get("label", action.get("id", ""))), str(action.get("id", "")))
        self.action_combo.blockSignals(False)
        enabled = self.action_combo.count() > 0 and atom_id is not None
        self.action_combo.setVisible(enabled)
        self.btn_apply_correction.setVisible(enabled)
        self.btn_apply_correction.setEnabled(enabled)

    def _emit_correction_requested(self) -> None:
        atom_id = self._selected_atom_id()
        action_id = self.action_combo.currentData()
        if atom_id is None or not action_id:
            return
        self.correction_requested.emit(int(atom_id), str(action_id))

    def report_rows(self) -> list[dict[str, object]]:
        """Devuelve el reporte de validación en formato tabular estructurado."""
        rows: list[dict[str, object]] = []
        for atom_id, issue in sorted(self._issues.items()):
            actions = self._issue_actions(int(atom_id))
            rows.append(
                {
                    "target_type": str(getattr(issue, "target_type", "atom") or "atom"),
                    "atom_id": int(atom_id),
                    "element": str(getattr(issue, "element", "") or ""),
                    "severity": str(getattr(issue, "severity", "error") or "error"),
                    "code": str(getattr(issue, "code", "") or ""),
                    "observed_valence": float(getattr(issue, "observed_valence", 0.0) or 0.0),
                    "bond_order_sum": float(getattr(issue, "bond_order_sum", 0.0) or 0.0),
                    "assigned_h": int(getattr(issue, "assigned_h", 0) or 0),
                    "implicit_h": int(getattr(issue, "implicit_h", 0) or 0),
                    "message": str(getattr(issue, "message", "") or ""),
                    "suggestion": str(getattr(issue, "hint", "") or getattr(issue, "suggestion", "") or ""),
                    "actions": ", ".join(str(action.get("label", "")) for action in actions),
                }
            )
        return rows

    def report_tsv(self) -> str:
        """Devuelve el reporte actual como TSV."""
        rows = self.report_rows()
        headers = [
            "target_type",
            "atom_id",
            "element",
            "severity",
            "code",
            "observed_valence",
            "bond_order_sum",
            "assigned_h",
            "implicit_h",
            "message",
            "suggestion",
            "actions",
        ]
        lines = ["\t".join(headers)]
        for row in rows:
            lines.append("\t".join(str(row.get(header, "")) for header in headers))
        return "\n".join(lines)

    def copy_report_tsv(self) -> str:
        """Copia el reporte de validación como TSV."""
        text = self.report_tsv()
        QApplication.clipboard().setText(text)
        return text

    def export_report(self, path: str | None = None) -> str:
        """Exporta el reporte como TSV o JSON según extensión."""
        if path is None:
            path, _ = QFileDialog.getSaveFileName(
                self,
                "Exportar reporte de validación",
                "",
                "TSV (*.tsv);;JSON (*.json)",
            )
        if not path:
            return ""
        target = Path(path)
        if target.suffix.lower() == ".json":
            target.write_text(
                json.dumps(self.report_rows(), ensure_ascii=False, indent=2) + "\n",
                encoding="utf-8",
            )
        else:
            target.write_text(self.report_tsv() + "\n", encoding="utf-8")
        return str(target)


class AppearanceDock(QDockWidget):
    """Dock para preferencias de apariencia."""

    appearance_changed = pyqtSignal(dict)

    def __init__(self, current_style: DrawingStyle, parent=None):
        """Inicializa el dock de apariencia con el estilo actual."""
        super().__init__("Apariencia", parent)
        self.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea)

        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(8, 8, 8, 8)

        form = QFormLayout()
        self.bond_cap_combo = QComboBox()
        self.bond_cap_combo.addItem("Redondeados", "round")
        self.bond_cap_combo.addItem("Rectos", "flat")
        current_cap = (
            "round"
            if current_style.cap_style == Qt.PenCapStyle.RoundCap
            else "flat"
        )
        index = self.bond_cap_combo.findData(current_cap)
        if index >= 0:
            self.bond_cap_combo.setCurrentIndex(index)
        self.bond_cap_combo.currentIndexChanged.connect(self._emit_change)
        form.addRow("Bordes de enlace", self.bond_cap_combo)

        layout.addLayout(form)
        layout.addStretch()
        self.setWidget(container)

    def _emit_change(self) -> None:
        """Emite cambios de apariencia al resto de la aplicación."""
        self.appearance_changed.emit({"bond_caps": self.bond_cap_combo.currentData()})

    def set_bond_caps(self, value: str) -> None:
        """Actualiza el selector de extremos de enlace."""
        index = self.bond_cap_combo.findData(value)
        if index >= 0:
            self.bond_cap_combo.blockSignals(True)
            self.bond_cap_combo.setCurrentIndex(index)
            self.bond_cap_combo.blockSignals(False)


# Backwards-compatible alias
TemplatesDock = PlantillasDock
