from PyQt6.QtWidgets import QToolBar
from PyQt6.QtGui import QAction, QActionGroup
from PyQt6.QtCore import pyqtSignal

class ChemusonToolbar(QToolBar):
    """
    Barra de herramientas para seleccionar instrumentos de dibujo.
    """
    # Señal que indica qué herramienta ha sido seleccionada
    tool_changed = pyqtSignal(str)

    def __init__(self) -> None:
        super().__init__("Herramientas de Dibujo")
        self.setMovable(False)
        self.setFloatable(False)

        # Grupo para que los botones sean mutuamente excluyentes
        self.action_group = QActionGroup(self)
        self.action_group.setExclusive(True)

        # Herramientas de átomos
        self._add_tool("C", "Átomo de Carbono", "atom_C")
        self._add_tool("N", "Átomo de Nitrógeno", "atom_N")
        self._add_tool("O", "Átomo de Oxígeno", "atom_O")
        self._add_tool("S", "Átomo de Azufre", "atom_S")

        self.addSeparator()

        # Herramientas de enlaces (por ahora simples identificadores)
        self._add_tool("Simple", "Enlace Sencillo", "bond_single")
        self._add_tool("Doble", "Enlace Doble", "bond_double")

        # Seleccionar Carbono por defecto
        carbon_action = self.findChild(QAction, "atom_C")
        if carbon_action:
            carbon_action.setChecked(True)

    def _add_tool(self, text: str, tooltip: str, internal_id: str) -> QAction:
        """Helper para añadir una acción a la toolbar."""
        action = QAction(text, self)
        action.setObjectName(internal_id)
        action.setToolTip(tooltip)
        action.setCheckable(True)
        
        self.action_group.addAction(action)
        self.addAction(action)
        
        # Conectar señal al ser activada
        action.triggered.connect(lambda: self.tool_changed.emit(internal_id))
        
        return action
