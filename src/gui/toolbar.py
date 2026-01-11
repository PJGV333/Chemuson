"""
Chemuson Toolbar
Vertical toolbar with vector icons for drawing tools.
"""
from PyQt6.QtWidgets import QToolBar, QWidget
from PyQt6.QtGui import QAction, QActionGroup
from PyQt6.QtCore import pyqtSignal, Qt, QSize

from gui.icons import (
    draw_atom_icon,
    draw_bond_icon,
    draw_generic_icon,
    draw_ring_icon,
    ATOM_COLORS,
)


class ChemusonToolbar(QToolBar):
    """
    Vertical toolbar for selecting drawing tools.
    Organized into groups: Tools, Bonds, Atoms.
    """
    # Signal emitted when a tool is selected
    tool_changed = pyqtSignal(str)

    def __init__(self, parent=None) -> None:
        super().__init__("Herramientas de Dibujo", parent)
        self.setOrientation(Qt.Orientation.Vertical)
        self.setMovable(False)
        self.setFloatable(False)
        self.setIconSize(QSize(32, 32))
        
        # Style the toolbar
        self.setStyleSheet("""
            QToolBar {
                background-color: #F5F5F5;
                border-right: 1px solid #CCCCCC;
                spacing: 2px;
                padding: 4px;
            }
            QToolButton {
                border: 1px solid transparent;
                border-radius: 4px;
                padding: 4px;
                margin: 1px;
            }
            QToolButton:hover {
                background-color: #E0E0E0;
                border: 1px solid #BBBBBB;
            }
            QToolButton:checked {
                background-color: #C8DFFF;
                border: 1px solid #6699CC;
            }
        """)

        # Exclusive action group
        self.action_group = QActionGroup(self)
        self.action_group.setExclusive(True)

        # === TOOLS GROUP ===
        self._add_section_label("Herramientas")
        self._add_tool_with_icon(
            draw_generic_icon('pointer'),
            "Seleccionar",
            "tool_select"
        )
        self._add_tool_with_icon(
            draw_generic_icon('eraser'),
            "Borrar",
            "tool_erase"
        )

        self.addSeparator()

        # === BONDS GROUP ===
        self._add_section_label("Enlaces")
        self._add_tool_with_icon(
            draw_bond_icon('single'),
            "Enlace Sencillo",
            "bond_single"
        )
        self._add_tool_with_icon(
            draw_bond_icon('double'),
            "Enlace Doble",
            "bond_double"
        )
        self._add_tool_with_icon(
            draw_ring_icon(),
            "Anillo de Benceno",
            "ring_benzene"
        )

        self.addSeparator()

        # === ATOMS GROUP ===
        self._add_section_label("Átomos")
        
        # Common atoms in order of usage
        atoms = ['C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br']
        for atom in atoms:
            self._add_tool_with_icon(
                draw_atom_icon(atom),
                f"Átomo de {self._get_atom_name(atom)} ({atom})",
                f"atom_{atom}"
            )

        # Select Carbon by default
        carbon_action = self.findChild(QAction, "atom_C")
        if carbon_action:
            carbon_action.setChecked(True)

    def _add_section_label(self, text: str) -> None:
        """Add a section label (non-interactive)."""
        # Using a simple separator approach for now
        # PyQt6 doesn't have addWidget directly for labels in toolbar easily
        pass  # Section labels handled by separators for simplicity

    def _add_tool_with_icon(self, icon, tooltip: str, internal_id: str) -> QAction:
        """Add a tool with an icon to the toolbar."""
        action = QAction(icon, "", self)
        action.setObjectName(internal_id)
        action.setToolTip(tooltip)
        action.setCheckable(True)
        
        self.action_group.addAction(action)
        self.addAction(action)
        
        # Connect signal when activated
        action.triggered.connect(lambda checked, id=internal_id: self.tool_changed.emit(id))
        
        return action

    def _get_atom_name(self, symbol: str) -> str:
        """Get the Spanish name of an element."""
        names = {
            'C': 'Carbono',
            'N': 'Nitrógeno',
            'O': 'Oxígeno',
            'S': 'Azufre',
            'P': 'Fósforo',
            'F': 'Flúor',
            'Cl': 'Cloro',
            'Br': 'Bromo',
            'H': 'Hidrógeno',
        }
        return names.get(symbol, symbol)
