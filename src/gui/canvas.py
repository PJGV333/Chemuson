from PyQt6.QtWidgets import QWidget
from PyQt6.QtGui import QPainter, QPen, QColor, QMouseEvent
from PyQt6.QtCore import Qt, QPoint
from core.molecule import Molecule

class ChemusonCanvas(QWidget):
    """
    Lienzo interactivo para dibujar moléculas.
    """
    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.molecule = Molecule()
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setBackgroundRole(QWidget.backgroundRole(self))
        self.setAutoFillBackground(True)
        
        # Estado básico para interacción
        self.last_pos = QPoint()
        self.selected_atom_idx: int | None = None
        self.hover_atom_idx: int | None = None
        self.current_tool: str = "atom_C" # Valor inicial por defecto
        self.setMouseTracking(True)
        
    def set_current_tool(self, tool_id: str) -> None:
        """Cambia la herramienta activa (C, N, O, enlace, etc)."""
        self.current_tool = tool_id
        # Limpiar selección al cambiar de herramienta para evitar enlaces accidentales
        self.selected_atom_idx = None
        self.update()

    def paintEvent(self, event) -> None:
        """Renderiza la molécula en el lienzo."""
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        
        # Fondo blanco
        painter.fillRect(self.rect(), Qt.GlobalColor.white)
        
        # Dibujar enlaces
        pen = QPen(Qt.GlobalColor.black, 2)
        painter.setPen(pen)
        for bond in self.molecule.bonds:
            a1 = self.molecule.atoms[bond.atom1_idx]
            a2 = self.molecule.atoms[bond.atom2_idx]
            painter.drawLine(int(a1.x), int(a1.y), int(a2.x), int(a2.y))
            
        # Dibujar átomos
        for i, atom in enumerate(self.molecule.atoms):
            # Resaltar si está seleccionado o hover
            is_selected = (i == self.selected_atom_idx)
            is_hover = (i == self.hover_atom_idx)
            
            if is_selected:
                painter.setBrush(QColor(200, 230, 255))
            elif is_hover:
                painter.setBrush(QColor(240, 240, 240))
            else:
                painter.setBrush(Qt.GlobalColor.white)
                
            painter.setPen(QPen(Qt.GlobalColor.black, 1))
            painter.drawEllipse(QPoint(int(atom.x), int(atom.y)), 10, 10)
            painter.drawText(QPoint(int(atom.x) - 4, int(atom.y) + 4), atom.symbol)

    def mousePressEvent(self, event: QMouseEvent) -> None:
        """Maneja el clic inicial."""
        pos = event.position()
        clicked_atom = self._get_atom_at(pos.x(), pos.y())

        if event.button() == Qt.MouseButton.LeftButton:
            if self.current_tool.startswith("atom_"):
                symbol = self.current_tool.split("_")[1]
                if clicked_atom is not None:
                    # Si ya hay un átomo, por ahora no hacemos nada o cambiamos su tipo
                    self.molecule.atoms[clicked_atom].symbol = symbol
                else:
                    # Añadir un átomo del tipo seleccionado
                    self.molecule.add_atom(symbol, pos.x(), pos.y())
                self.selected_atom_idx = None
            
            elif self.current_tool.startswith("bond_"):
                bond_type = 1 if "single" in self.current_tool else 2
                if clicked_atom is not None:
                    if self.selected_atom_idx is None:
                        self.selected_atom_idx = clicked_atom
                    elif self.selected_atom_idx != clicked_atom:
                        self.molecule.add_bond(self.selected_atom_idx, clicked_atom, order=bond_type)
                        self.selected_atom_idx = None
                    else:
                        self.selected_atom_idx = None
                else:
                    self.selected_atom_idx = None
            
            self.update()

    def mouseMoveEvent(self, event: QMouseEvent) -> None:
        """Maneja el movimiento del mouse para hover."""
        pos = event.position()
        self.hover_atom_idx = self._get_atom_at(pos.x(), pos.y())
        self.update()

    def _get_atom_at(self, x: float, y: float) -> int | None:
        """Retorna el índice del átomo en la posición (x, y) o None."""
        for i, atom in enumerate(self.molecule.atoms):
            dx = atom.x - x
            dy = atom.y - y
            if (dx*dx + dy*dy) < 15*15: # Radio de colisión un poco más grande que el dibujo
                return i
        return None
