"""
Chemuson Canvas
Page-based canvas using QGraphicsView/QGraphicsScene for document-style editing.
"""
from PyQt6.QtWidgets import (
    QGraphicsView, QGraphicsScene, QGraphicsRectItem,
    QGraphicsEllipseItem, QGraphicsLineItem, QGraphicsTextItem,
    QGraphicsDropShadowEffect, QGraphicsItem
)
from PyQt6.QtGui import (
    QPainter, QPen, QColor, QBrush, QFont, QWheelEvent
)
from PyQt6.QtCore import Qt, QRectF, QPointF

from core.molecule import Molecule


# Paper dimensions (A4, approximately)
PAPER_WIDTH = 800
PAPER_HEIGHT = 1000
PAPER_MARGIN = 40  # Margins where atoms shouldn't be placed


class AtomItem(QGraphicsEllipseItem):
    """Graphics item representing an atom."""
    
    def __init__(self, symbol: str, x: float, y: float, atom_idx: int):
        radius = 12
        super().__init__(-radius, -radius, radius * 2, radius * 2)
        self.setPos(x, y)
        self.atom_idx = atom_idx
        self.symbol = symbol
        
        # Visual style
        self.setBrush(QBrush(Qt.GlobalColor.white))
        self.setPen(QPen(QColor('#333333'), 1.5))
        
        # Make selectable
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        
        # Add text label
        self.label = QGraphicsTextItem(symbol, self)
        font = QFont('Arial', 10, QFont.Weight.Bold)
        self.label.setFont(font)
        # Center the label
        rect = self.label.boundingRect()
        self.label.setPos(-rect.width() / 2, -rect.height() / 2)
        self.label.setDefaultTextColor(QColor('#333333'))
    
    def set_selected_style(self, selected: bool):
        """Change visual style based on selection."""
        if selected:
            self.setBrush(QBrush(QColor('#C8DFFF')))
            self.setPen(QPen(QColor('#4477AA'), 2))
        else:
            self.setBrush(QBrush(Qt.GlobalColor.white))
            self.setPen(QPen(QColor('#333333'), 1.5))
    
    def set_hover_style(self, hover: bool):
        """Change visual style based on hover."""
        if hover:
            self.setBrush(QBrush(QColor('#F0F0F0')))
        else:
            self.setBrush(QBrush(Qt.GlobalColor.white))


class BondItem(QGraphicsLineItem):
    """Graphics item representing a bond."""
    
    def __init__(self, x1: float, y1: float, x2: float, y2: float, order: int = 1):
        super().__init__(x1, y1, x2, y2)
        self.order = order
        self.setPen(QPen(QColor('#333333'), 2 if order == 1 else 3))
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)


class ChemusonCanvas(QGraphicsView):
    """
    Page-based canvas for drawing molecules.
    Uses QGraphicsView/QGraphicsScene with a centered paper sheet.
    """
    
    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        
        # Create scene
        self.scene = QGraphicsScene()
        self.setScene(self.scene)
        
        # Molecule data
        self.molecule = Molecule()
        self.atom_items: list[AtomItem] = []
        self.bond_items: list[BondItem] = []
        
        # Interaction state
        self.current_tool: str = "atom_C"
        self.selected_atom_idx: int | None = None
        self.hover_atom_item: AtomItem | None = None
        
        # Setup the view
        self._setup_view()
        self._create_paper()
        
        # Enable mouse tracking for hover effects
        self.setMouseTracking(True)
        self.viewport().setMouseTracking(True)
    
    def _setup_view(self):
        """Configure the graphics view."""
        # Gray background for the "desktop"
        self.setBackgroundBrush(QBrush(QColor('#E0E0E0')))
        
        # Enable antialiasing
        self.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.setRenderHint(QPainter.RenderHint.TextAntialiasing)
        
        # Scroll behavior
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        
        # Set scene rect (allow some space around paper)
        margin = 100
        self.scene.setSceneRect(
            -margin, -margin,
            PAPER_WIDTH + 2 * margin,
            PAPER_HEIGHT + 2 * margin
        )
        
        # Drag mode for panning (when select tool is active)
        self.setDragMode(QGraphicsView.DragMode.NoDrag)
        
        # Zoom settings
        self._zoom_factor = 1.0
        self._min_zoom = 0.25
        self._max_zoom = 4.0
    
    def _create_paper(self):
        """Create the white paper rectangle with shadow."""
        # Paper rectangle
        self.paper = QGraphicsRectItem(0, 0, PAPER_WIDTH, PAPER_HEIGHT)
        self.paper.setBrush(QBrush(Qt.GlobalColor.white))
        self.paper.setPen(QPen(QColor('#CCCCCC'), 1))
        self.paper.setZValue(-10)  # Behind everything
        
        # Shadow effect
        shadow = QGraphicsDropShadowEffect()
        shadow.setBlurRadius(20)
        shadow.setColor(QColor(0, 0, 0, 80))
        shadow.setOffset(5, 5)
        self.paper.setGraphicsEffect(shadow)
        
        self.scene.addItem(self.paper)
        
        # Center view on paper
        self.centerOn(PAPER_WIDTH / 2, PAPER_HEIGHT / 2)
    
    def set_current_tool(self, tool_id: str) -> None:
        """Change the active tool."""
        self.current_tool = tool_id
        self.selected_atom_idx = None
        
        # Clear selection when changing tools
        self.scene.clearSelection()
        
        # Update cursor based on tool
        if tool_id == 'tool_select':
            self.setCursor(Qt.CursorShape.ArrowCursor)
            self.setDragMode(QGraphicsView.DragMode.RubberBandDrag)
        elif tool_id == 'tool_erase':
            self.setCursor(Qt.CursorShape.CrossCursor)
            self.setDragMode(QGraphicsView.DragMode.NoDrag)
        else:
            self.setCursor(Qt.CursorShape.CrossCursor)
            self.setDragMode(QGraphicsView.DragMode.NoDrag)
    
    def mousePressEvent(self, event) -> None:
        """Handle mouse click on the canvas."""
        scene_pos = self.mapToScene(event.pos())
        
        # Check if click is within paper bounds
        if not self._is_on_paper(scene_pos.x(), scene_pos.y()):
            super().mousePressEvent(event)
            return
        
        if event.button() == Qt.MouseButton.LeftButton:
            clicked_atom = self._get_atom_at(scene_pos.x(), scene_pos.y())
            
            if self.current_tool.startswith("atom_"):
                symbol = self.current_tool.split("_")[1]
                if clicked_atom is not None:
                    # Change existing atom type
                    self.molecule.atoms[clicked_atom].symbol = symbol
                    self.atom_items[clicked_atom].symbol = symbol
                    self.atom_items[clicked_atom].label.setPlainText(symbol)
                else:
                    # Add new atom
                    idx = self.molecule.add_atom(symbol, scene_pos.x(), scene_pos.y())
                    self._add_atom_item(symbol, scene_pos.x(), scene_pos.y(), idx)
                self.selected_atom_idx = None
            
            elif self.current_tool.startswith("bond_"):
                bond_type = 1 if "single" in self.current_tool else 2
                if clicked_atom is not None:
                    if self.selected_atom_idx is None:
                        self.selected_atom_idx = clicked_atom
                        self.atom_items[clicked_atom].set_selected_style(True)
                    elif self.selected_atom_idx != clicked_atom:
                        # Create bond
                        self.molecule.add_bond(self.selected_atom_idx, clicked_atom, order=bond_type)
                        self._add_bond_item(self.selected_atom_idx, clicked_atom, bond_type)
                        self.atom_items[self.selected_atom_idx].set_selected_style(False)
                        self.selected_atom_idx = None
                    else:
                        # Clicked same atom, deselect
                        self.atom_items[self.selected_atom_idx].set_selected_style(False)
                        self.selected_atom_idx = None
                else:
                    # Clicked empty space while bonding
                    if self.selected_atom_idx is not None:
                        self.atom_items[self.selected_atom_idx].set_selected_style(False)
                    self.selected_atom_idx = None
            
            elif self.current_tool == "tool_erase":
                if clicked_atom is not None:
                    self._remove_atom(clicked_atom)
            
            elif self.current_tool == "tool_select":
                super().mousePressEvent(event)
        else:
            super().mousePressEvent(event)
    
    def mouseMoveEvent(self, event) -> None:
        """Handle mouse movement for hover effects."""
        scene_pos = self.mapToScene(event.pos())
        
        # Update hover state
        hover_atom = self._get_atom_at(scene_pos.x(), scene_pos.y())
        
        # Clear previous hover
        if self.hover_atom_item is not None:
            if self.hover_atom_item.atom_idx != self.selected_atom_idx:
                self.hover_atom_item.set_hover_style(False)
            self.hover_atom_item = None
        
        # Set new hover
        if hover_atom is not None and hover_atom < len(self.atom_items):
            item = self.atom_items[hover_atom]
            if item.atom_idx != self.selected_atom_idx:
                item.set_hover_style(True)
            self.hover_atom_item = item
        
        super().mouseMoveEvent(event)
    
    def wheelEvent(self, event: QWheelEvent) -> None:
        """Handle mouse wheel for zooming."""
        # Zoom in/out
        if event.angleDelta().y() > 0:
            self.zoom_in()
        else:
            self.zoom_out()
    
    def zoom_in(self):
        """Zoom in the view."""
        if self._zoom_factor < self._max_zoom:
            self._zoom_factor *= 1.2
            self.scale(1.2, 1.2)
    
    def zoom_out(self):
        """Zoom out the view."""
        if self._zoom_factor > self._min_zoom:
            self._zoom_factor /= 1.2
            self.scale(1 / 1.2, 1 / 1.2)
    
    def clear_canvas(self):
        """Clear all atoms and bonds from the canvas."""
        # Remove all items except paper
        for item in self.atom_items:
            self.scene.removeItem(item)
        for item in self.bond_items:
            self.scene.removeItem(item)
        
        self.atom_items.clear()
        self.bond_items.clear()
        self.molecule.clear()
        self.selected_atom_idx = None
    
    def _add_atom_item(self, symbol: str, x: float, y: float, idx: int):
        """Add an atom graphics item."""
        item = AtomItem(symbol, x, y, idx)
        self.scene.addItem(item)
        self.atom_items.append(item)
    
    def _add_bond_item(self, idx1: int, idx2: int, order: int):
        """Add a bond graphics item between two atoms."""
        a1 = self.molecule.atoms[idx1]
        a2 = self.molecule.atoms[idx2]
        item = BondItem(a1.x, a1.y, a2.x, a2.y, order)
        item.setZValue(-5)  # Behind atoms but above paper
        self.scene.addItem(item)
        self.bond_items.append(item)
    
    def _remove_atom(self, idx: int):
        """Remove an atom and its connected bonds."""
        if idx >= len(self.atom_items):
            return
        
        # Remove atom item
        item = self.atom_items[idx]
        self.scene.removeItem(item)
        
        # This is a simplified removal - a full implementation would
        # need to properly update indices and remove connected bonds
        # For now, we just hide it
        self.molecule.atoms[idx].symbol = ""
    
    def _get_atom_at(self, x: float, y: float) -> int | None:
        """Find atom index at given position."""
        for i, atom in enumerate(self.molecule.atoms):
            if atom.symbol == "":  # Skip deleted atoms
                continue
            dx = atom.x - x
            dy = atom.y - y
            if (dx * dx + dy * dy) < 20 * 20:  # Collision radius
                return i
        return None
    
    def _is_on_paper(self, x: float, y: float) -> bool:
        """Check if coordinates are within paper bounds."""
        return (PAPER_MARGIN <= x <= PAPER_WIDTH - PAPER_MARGIN and
                PAPER_MARGIN <= y <= PAPER_HEIGHT - PAPER_MARGIN)
