"""
Docks de Chemuson.

Paneles laterales para plantillas, inspección de propiedades y apariencia.
"""
from PyQt6.QtWidgets import (
    QDockWidget,
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
)
from PyQt6.QtCore import Qt, pyqtSignal
from gui.style import DrawingStyle


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
        if num_atoms == 0 and num_bonds == 0 and num_text == 0:
            self.info_label.setText("Nada seleccionado")
            self.info_label.setVisible(True)
            self.prop_table.setVisible(False)
            return
            
        self.info_label.setVisible(False)
        self.prop_table.setVisible(True)
        self.prop_table.setRowCount(0)
        
        data = []
        if num_atoms == 1 and num_bonds == 0 and num_text == 0:
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
