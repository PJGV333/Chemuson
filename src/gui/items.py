"""
Chemuson Scene Items
QGraphicsItem subclasses for atoms and bonds with professional rendering.
"""
from __future__ import annotations

import math
from PyQt6.QtWidgets import (
    QGraphicsEllipseItem,
    QGraphicsPathItem,
    QGraphicsTextItem,
    QGraphicsItem,
)
from PyQt6.QtGui import QColor, QFont, QPainterPath, QPen, QBrush
from PyQt6.QtCore import Qt

from core.model import Atom, Bond, BondStyle


# Element colors for heteroatoms (CPK coloring scheme)
ELEMENT_COLORS = {
    'C': '#333333',   # Carbon - dark gray
    'N': '#3050F8',   # Nitrogen - blue
    'O': '#FF0D0D',   # Oxygen - red
    'S': '#FFFF30',   # Sulfur - yellow
    'P': '#FF8000',   # Phosphorus - orange
    'F': '#90E050',   # Fluorine - light green
    'Cl': '#1FF01F',  # Chlorine - green
    'Br': '#A62929',  # Bromine - dark red
    'I': '#940094',   # Iodine - purple
    'H': '#FFFFFF',   # Hydrogen - white
}


class AtomItem(QGraphicsEllipseItem):
    """Graphics item representing an atom with optional implicit hiding."""
    
    def __init__(
        self,
        atom: Atom,
        radius: float = 12.0,
        show_carbon: bool = True,
        show_hydrogen: bool = True,
    ) -> None:
        super().__init__(-radius, -radius, radius * 2, radius * 2)
        self.atom_id = atom.id
        self.element = atom.element
        self._radius = radius
        self._show_carbon = show_carbon
        self._show_hydrogen = show_hydrogen
        self._is_selected = False
        self._is_hover = False
        
        self.setPos(atom.x, atom.y)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        
        # Create text label
        self.label = QGraphicsTextItem(atom.element, self)
        font = QFont("Arial", 10, QFont.Weight.Bold)
        self.label.setFont(font)
        self._center_label()
        
        # Set color based on element
        self._set_element_color()
        
        # Apply visibility
        self._update_visibility()
    
    def _center_label(self) -> None:
        """Center the label text within the atom circle."""
        rect = self.label.boundingRect()
        self.label.setPos(-rect.width() / 2, -rect.height() / 2)
    
    def _set_element_color(self) -> None:
        """Set the label color based on element type."""
        color = ELEMENT_COLORS.get(self.element, '#333333')
        self.label.setDefaultTextColor(QColor(color))
    
    def _update_visibility(self) -> None:
        """Hide circle/label for implicit C or H based on settings."""
        hide_element = False
        
        # Check if this is an implicit carbon or hydrogen
        if self.element == "C" and not self._show_carbon:
            hide_element = True
        elif self.element == "H" and not self._show_hydrogen:
            hide_element = True
        
        if hide_element:
            # Make circle invisible but keep selectable (minimal hit area)
            self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
            self.setPen(QPen(Qt.PenStyle.NoPen))
            self.label.setVisible(False)
        else:
            # Normal visible atom with proper styling
            self._apply_normal_style()
            self.label.setVisible(True)
    
    def _apply_normal_style(self) -> None:
        """Apply normal (non-hidden) styling based on selection state."""
        if self._is_selected:
            self.setBrush(QBrush(QColor("#C8DFFF")))
            self.setPen(QPen(QColor("#4477AA"), 2))
        elif self._is_hover:
            self.setBrush(QBrush(QColor("#F0F0F0")))
            self.setPen(QPen(QColor("#333333"), 1.5))
        else:
            self.setBrush(QBrush(Qt.GlobalColor.white))
            self.setPen(QPen(QColor("#333333"), 1.5))
    
    def set_visibility_flags(self, show_carbon: bool, show_hydrogen: bool) -> None:
        """Update visibility flags and refresh display."""
        self._show_carbon = show_carbon
        self._show_hydrogen = show_hydrogen
        self._update_visibility()
    
    def set_selected(self, selected: bool) -> None:
        """Set selection state."""
        self._is_selected = selected
        self._update_visibility()
    
    def set_hover(self, hover: bool) -> None:
        """Set hover state."""
        self._is_hover = hover
        # Only update if currently visible
        is_hidden = (
            (self.element == "C" and not self._show_carbon) or
            (self.element == "H" and not self._show_hydrogen)
        )
        if not is_hidden:
            self._apply_normal_style()
    
    def set_element(self, element: str) -> None:
        """Change the element of this atom."""
        self.element = element
        self.label.setPlainText(element)
        self._center_label()
        self._set_element_color()
        self._update_visibility()


class BondItem(QGraphicsPathItem):
    """Graphics item representing a chemical bond."""
    
    def __init__(self, bond: Bond, atom1: Atom, atom2: Atom, render_aromatic_as_circle: bool = True) -> None:
        super().__init__()
        self.bond_id = bond.id
        self.a1_id = bond.a1_id
        self.a2_id = bond.a2_id
        self.order = bond.order
        self.style = bond.style
        self.stereo = bond.stereo
        self.is_aromatic = bond.is_aromatic
        self.render_aromatic_as_circle = render_aromatic_as_circle
        self.setZValue(-5)
        self.setPen(QPen(QColor("#333333"), 2))
        self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.update_positions(atom1, atom2)

    def set_bond(self, bond: Bond, atom1: Atom, atom2: Atom) -> None:
        """Update bond properties and redraw."""
        self.order = bond.order
        self.style = bond.style
        self.stereo = bond.stereo
        self.is_aromatic = bond.is_aromatic
        self.update_positions(atom1, atom2)

    def set_render_aromatic_as_circle(self, enabled: bool) -> None:
        """Toggle aromatic rendering mode."""
        self.render_aromatic_as_circle = enabled
        # We need atom positions to update, but we don't store atoms directly.
        # However, update_positions is usually called by canvas.
        # Here we can't easily trigger redraw without atoms.
        # We will rely on canvas calling update_positions or us storing current positions?
        # Better: canvas iterates and calls update_positions_from_cache or similar.
        # Actually, let's just flag it; the canvas refresh method will handle triggering update.
        pass 

    def update_positions(self, atom1: Atom, atom2: Atom) -> None:
        """Redraw the bond path based on atom positions and bond type."""
        x1, y1 = atom1.x, atom1.y
        x2, y2 = atom2.x, atom2.y
        path = QPainterPath()
        dx = x2 - x1
        dy = y2 - y1
        length = math.hypot(dx, dy) or 1.0
        nx = -dy / length  # Normal vector perpendicular to bond
        ny = dx / length

        # Aromatic bonds: if circle mode, draw as single line
        if self.is_aromatic and self.render_aromatic_as_circle:
            path.moveTo(x1, y1)
            path.lineTo(x2, y2)
            self.setPen(QPen(QColor("#333333"), 2))
            self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
            self.setPath(path)
            return

        if self.style == BondStyle.PLAIN:
            if self.order == 1:
                path.moveTo(x1, y1)
                path.lineTo(x2, y2)
            else:
                offset = 3.0
                if self.order == 2:
                    path.moveTo(x1 + nx * offset, y1 + ny * offset)
                    path.lineTo(x2 + nx * offset, y2 + ny * offset)
                    path.moveTo(x1 - nx * offset, y1 - ny * offset)
                    path.lineTo(x2 - nx * offset, y2 - ny * offset)
                else:  # Triple bond
                    path.moveTo(x1, y1)
                    path.lineTo(x2, y2)
                    path.moveTo(x1 + nx * offset, y1 + ny * offset)
                    path.lineTo(x2 + nx * offset, y2 + ny * offset)
                    path.moveTo(x1 - nx * offset, y1 - ny * offset)
                    path.lineTo(x2 - nx * offset, y2 - ny * offset)
            self.setPen(QPen(QColor("#333333"), 2))
            self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
            
        elif self.style == BondStyle.WEDGE:
            width = 10.0
            p1x = x2 + nx * width
            p1y = y2 + ny * width
            p2x = x2 - nx * width
            p2y = y2 - ny * width
            path.moveTo(x1, y1)
            path.lineTo(p1x, p1y)
            path.lineTo(p2x, p2y)
            path.closeSubpath()
            self.setPen(QPen(QColor("#333333"), 1))
            self.setBrush(QBrush(QColor("#333333")))
            
        elif self.style == BondStyle.HASHED:
            steps = 7
            max_width = 10.0
            for i in range(1, steps + 1):
                t = i / (steps + 1)
                px = x1 + dx * t
                py = y1 + dy * t
                width = max_width * t
                path.moveTo(px + nx * width, py + ny * width)
                path.lineTo(px - nx * width, py - ny * width)
            self.setPen(QPen(QColor("#333333"), 1.5))
            self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
            
        elif self.style == BondStyle.WAVY:
            segments = max(6, int(length / 6))
            amplitude = 3.0
            for i in range(segments + 1):
                t = i / segments
                px = x1 + dx * t
                py = y1 + dy * t
                offset = math.sin(t * math.pi * 4) * amplitude
                wx = px + nx * offset
                wy = py + ny * offset
                if i == 0:
                    path.moveTo(wx, wy)
                else:
                    path.lineTo(wx, wy)
            self.setPen(QPen(QColor("#333333"), 1.5))
            self.setBrush(QBrush(Qt.BrushStyle.NoBrush))

        self.setPath(path)


class AromaticCircleItem(QGraphicsEllipseItem):
    """A circle drawn inside aromatic rings to indicate aromaticity."""
    
    def __init__(self, center_x: float, center_y: float, radius: float) -> None:
        super().__init__(
            center_x - radius, center_y - radius,
            radius * 2, radius * 2
        )
        self.setPen(QPen(QColor("#333333"), 1.5))
        self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
        self.setZValue(-10)  # Behind bonds and atoms
