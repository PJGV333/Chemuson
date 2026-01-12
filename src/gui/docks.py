"""
Chemuson Dock Widgets
Interactive panels for templates and property inspection.
"""
from PyQt6.QtWidgets import (
    QDockWidget,
    QWidget,
    QVBoxLayout,
    QLabel,
    QListWidget,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
)
from PyQt6.QtCore import Qt

class TemplatesDock(QDockWidget):
    """
    Dock widget displaying a list of chemical templates.
    """
    def __init__(self, parent=None):
        super().__init__("Plantillas", parent)
        self.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea)
        
        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)
        
        self.list_widget = QListWidget()
        self.list_widget.addItems([
            "Benceno",
            "Ciclohexano (Silla)",
            "Ciclohexano (Bote)",
            "Ciclopentadieno",
            "Nutraceúticos básicos",
            "Aminoácidos",
        ])
        layout.addWidget(self.list_widget)
        
        self.setWidget(container)


class InspectorDock(QDockWidget):
    """
    Dock widget displaying properties of the selected object.
    """
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

    def update_selection(self, num_atoms: int, num_bonds: int, details: dict = None):
        """Update the inspector with selection info."""
        if num_atoms == 0 and num_bonds == 0:
            self.info_label.setText("Nada seleccionado")
            self.info_label.setVisible(True)
            self.prop_table.setVisible(False)
            return
            
        self.info_label.setVisible(False)
        self.prop_table.setVisible(True)
        self.prop_table.setRowCount(0)
        
        data = []
        if num_atoms == 1 and num_bonds == 0:
            data = [
                ("Tipo", "Átomo"),
                ("Elemento", details.get("element", "?")),
                ("ID", str(details.get("id", "?"))),
                ("Carga", str(details.get("charge", 0))),
                ("X", f"{details.get('x', 0):.1f}"),
                ("Y", f"{details.get('y', 0):.1f}"),
            ]
        elif num_atoms == 0 and num_bonds == 1:
            data = [
                ("Tipo", "Enlace"),
                ("Orden", str(details.get("order", 1))),
                ("Estilo", str(details.get("style", "Plain"))),
                ("Aromático", "Sí" if details.get("aromatic") else "No"),
            ]
        else:
            data = [
                ("Selección", "Múltiple"),
                ("Átomos", str(num_atoms)),
                ("Enlaces", str(num_bonds)),
            ]
            
        self.prop_table.setRowCount(len(data))
        for i, (key, val) in enumerate(data):
            self.prop_table.setItem(i, 0, QTableWidgetItem(key))
            self.prop_table.setItem(i, 1, QTableWidgetItem(val))
