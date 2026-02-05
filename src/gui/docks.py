"""
Chemuson Dock Widgets
Interactive panels for templates and property inspection.
"""
from PyQt6.QtWidgets import (
    QDockWidget,
    QWidget,
    QVBoxLayout,
    QLabel,
    QFormLayout,
    QComboBox,
    QTreeWidget,
    QTreeWidgetItem,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
)
from PyQt6.QtCore import Qt, pyqtSignal
from gui.style import DrawingStyle


class PlantillasDock(QDockWidget):
    """
    Dock widget displaying a list of chemical templates.
    """
    template_selected = pyqtSignal(dict)

    def __init__(self, parent=None):
        super().__init__("Plantillas", parent)
        self.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea)
        
        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)
        
        self.tree = QTreeWidget()
        self.tree.setHeaderHidden(True)
        self._populate_tree()
        self.tree.itemActivated.connect(self._emit_template)
        layout.addWidget(self.tree)
        
        self.setWidget(container)

    def _populate_tree(self) -> None:
        groups = {
            "Grupos funcionales": [
                {"name": "Alcohol", "type": "alcohol"},
                {"name": "Aldehído", "type": "aldehyde"},
                {"name": "Cetona", "type": "ketone"},
            ],
            "Aminoácidos": [
                {"name": "Glicina", "type": "glycine"},
                {"name": "Alanina", "type": "alanine"},
            ],
            "Protecciones": [
                {"name": "Boc", "type": "boc"},
                {"name": "Fmoc", "type": "fmoc"},
            ],
        }
        for group_name, items in groups.items():
            group_item = QTreeWidgetItem([group_name])
            group_item.setFlags(group_item.flags() & ~Qt.ItemFlag.ItemIsSelectable)
            for item in items:
                child = QTreeWidgetItem([item["name"]])
                child.setData(0, Qt.ItemDataRole.UserRole, item)
                group_item.addChild(child)
            self.tree.addTopLevelItem(group_item)
            group_item.setExpanded(True)

    def _emit_template(self, item: QTreeWidgetItem) -> None:
        payload = item.data(0, Qt.ItemDataRole.UserRole)
        if isinstance(payload, dict):
            self.template_selected.emit(payload)


class InspectorDock(QDockWidget):
    """
    Dock widget displaying properties of the selected object.
    """
    property_changed = pyqtSignal(str, object)

    def __init__(self, parent=None):
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
        self.update_selection(1, 0, {
            "element": atom.element,
            "id": atom.id,
            "charge": atom.charge,
            "x": atom.x,
            "y": atom.y,
        })

    def set_bond(self, bond) -> None:
        self.update_selection(0, 1, 0, {
            "order": bond.order,
            "style": bond.style,
            "aromatic": bond.is_aromatic,
        })

    def update_selection(self, num_atoms: int, num_bonds: int, num_text: int, details: dict):
        """Update the inspector with selection info."""
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
    """Dock widget for appearance preferences."""

    appearance_changed = pyqtSignal(dict)

    def __init__(self, current_style: DrawingStyle, parent=None):
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
        self.appearance_changed.emit({"bond_caps": self.bond_cap_combo.currentData()})

    def set_bond_caps(self, value: str) -> None:
        index = self.bond_cap_combo.findData(value)
        if index >= 0:
            self.bond_cap_combo.blockSignals(True)
            self.bond_cap_combo.setCurrentIndex(index)
            self.bond_cap_combo.blockSignals(False)


# Backwards-compatible alias
TemplatesDock = PlantillasDock
