"""
Chemuson Scene Items
QGraphicsItem subclasses for atoms and bonds with professional rendering.
"""
from __future__ import annotations

import math
from typing import Optional
from PyQt6.QtWidgets import (
    QGraphicsEllipseItem,
    QGraphicsPathItem,
    QGraphicsTextItem,
    QGraphicsItem,
)
from PyQt6.QtGui import QColor, QFont, QPainterPath, QPen, QBrush
from PyQt6.QtCore import Qt, QRectF, QPointF

from core.model import Atom, Bond, BondStyle
from gui.style import DrawingStyle, CHEMDOODLE_LIKE


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
        style: DrawingStyle = CHEMDOODLE_LIKE,
    ) -> None:
        super().__init__(-radius, -radius, radius * 2, radius * 2)
        self.atom_id = atom.id
        self.element = atom.element
        self._radius = radius
        self._show_carbon = show_carbon
        self._show_hydrogen = show_hydrogen
        self._is_explicit = atom.is_explicit
        self._is_selected = False
        self._is_hover = False
        self._valence_error = False
        self._element_color = QColor(ELEMENT_COLORS.get(atom.element, "#333333"))
        self._charge = atom.charge
        self._style = style
        
        self.setPos(atom.x, atom.y)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        
        # Create text label
        self.label = QGraphicsTextItem(atom.element, self)
        font = QFont("Arial", 10, QFont.Weight.Bold)
        self.label.setFont(font)

        # Charge label
        self.charge_label = QGraphicsTextItem("", self)
        charge_font = QFont("Arial", 8, QFont.Weight.Bold)
        self.charge_label.setFont(charge_font)
        self.charge_label.setDefaultTextColor(QColor("#333333"))

        self._center_label()
        self._update_charge_label()
        
        # Set color based on element
        self._set_element_color()
        
        # Apply visibility
        self._update_visibility()
    
    def _center_label(self) -> None:
        """Center the label text within the atom circle."""
        rect = self.label.boundingRect()
        self.label.setPos(-rect.width() / 2, -rect.height() / 2)
        if hasattr(self, "charge_label"):
            self._position_charge_label()

    def _position_charge_label(self) -> None:
        rect = self.label.boundingRect()
        x = rect.width() / 2 + 2
        y = -rect.height() / 2 - 2
        self.charge_label.setPos(x, y)
    
    def _set_element_color(self) -> None:
        """Set the label color based on element type."""
        color = ELEMENT_COLORS.get(self.element, "#333333")
        self._element_color = QColor(color)
        self._apply_label_style()

    def _apply_label_style(self) -> None:
        font = self.label.font()
        font.setUnderline(self._valence_error)
        self.label.setFont(font)
        if self._valence_error:
            self.label.setDefaultTextColor(QColor("#C0392B"))
        else:
            self.label.setDefaultTextColor(self._element_color)
    
    def _should_hide_element(self) -> bool:
        if self._is_explicit:
            return False
        if self.element == "C" and not self._show_carbon:
            return True
        if self.element == "H" and not self._show_hydrogen:
            return True
        return False

    def _should_draw_circle(self) -> bool:
        if self.element == "H":
            return False
        return not (self._is_explicit and self.element not in {"C", "H"})

    def _update_visibility(self) -> None:
        """Hide circle/label for implicit C or H based on settings."""
        if self._should_hide_element():
            # Make circle invisible but keep selectable (minimal hit area)
            self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
            self.setPen(QPen(Qt.PenStyle.NoPen))
            self.label.setVisible(False)
            return

        self.label.setVisible(True)
        self._apply_normal_style()
    
    def _apply_normal_style(self) -> None:
        """Apply normal (non-hidden) styling based on selection state."""
        if not self._should_draw_circle():
            self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
            self.setPen(QPen(Qt.PenStyle.NoPen))
            return
        if self._is_selected:
            self.setBrush(QBrush(QColor("#C8DFFF")))
            pen = QPen(QColor("#4477AA"), self._style.stroke_px)
            pen.setCapStyle(self._style.cap_style)
            pen.setJoinStyle(self._style.join_style)
            self.setPen(pen)
        elif self._is_hover:
            self.setBrush(QBrush(QColor("#F0F0F0")))
            pen = QPen(QColor("#333333"), self._style.stroke_px)
            pen.setCapStyle(self._style.cap_style)
            pen.setJoinStyle(self._style.join_style)
            self.setPen(pen)
        else:
            self.setBrush(QBrush(QColor(self._style.atom_fill_color)))
            pen = QPen(QColor(self._style.atom_stroke_color), self._style.stroke_px)
            pen.setCapStyle(self._style.cap_style)
            pen.setJoinStyle(self._style.join_style)
            self.setPen(pen)
    
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
        if not self._should_hide_element():
            self._apply_normal_style()
    
    def set_element(self, element: str, is_explicit: Optional[bool] = None) -> None:
        """Change the element of this atom."""
        self.element = element
        if is_explicit is not None:
            self._is_explicit = is_explicit
        self.label.setPlainText(element)
        self._center_label()
        self._set_element_color()
        self._update_visibility()

    def set_valence_error(self, has_error: bool) -> None:
        """Toggle valence error underline on the label."""
        self._valence_error = has_error
        self._apply_label_style()

    def set_charge(self, charge: int) -> None:
        """Update charge label display."""
        self._charge = charge
        if charge == 0:
            self.charge_label.setVisible(False)
            return
        sign = "+" if charge > 0 else "-"
        magnitude = abs(charge)
        label = sign if magnitude == 1 else f"{sign}{magnitude}"
        self.charge_label.setPlainText(label)
        self._position_charge_label()
        self.charge_label.setVisible(True)

    def _update_charge_label(self) -> None:
        self.set_charge(self._charge)

    def set_style(self, style: DrawingStyle) -> None:
        """Update drawing style and refresh visuals."""
        self._style = style
        if not self._should_hide_element():
            self._apply_normal_style()


class BondItem(QGraphicsPathItem):
    """Graphics item representing a chemical bond."""
    
    def __init__(
        self,
        bond: Bond,
        atom1: Atom,
        atom2: Atom,
        render_aromatic_as_circle: bool = True,
        style: DrawingStyle = CHEMDOODLE_LIKE,
    ) -> None:
        super().__init__()
        self.bond_id = bond.id
        self.a1_id = bond.a1_id
        self.a2_id = bond.a2_id
        self.order = bond.order
        self.style = bond.style
        self.stereo = bond.stereo
        self.is_aromatic = bond.is_aromatic
        self.display_order = bond.display_order
        self.ring_id = bond.ring_id
        self.length_px = bond.length_px
        self.render_aromatic_as_circle = render_aromatic_as_circle
        self._style = style
        self._offset_sign = 1
        self._ring_center: QPointF | None = None
        self.setZValue(-5)
        pen = QPen(QColor(self._style.bond_color), self._style.stroke_px)
        pen.setCapStyle(self._style.cap_style)
        pen.setJoinStyle(self._style.join_style)
        self.setPen(pen)
        self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.update_positions(atom1, atom2)

    def set_bond(self, bond: Bond, atom1: Atom, atom2: Atom) -> None:
        """Update bond properties and redraw."""
        self.order = bond.order
        self.style = bond.style
        self.stereo = bond.stereo
        self.is_aromatic = bond.is_aromatic
        self.display_order = bond.display_order
        self.ring_id = bond.ring_id
        self.length_px = bond.length_px
        self.update_positions(atom1, atom2)

    def set_render_aromatic_as_circle(self, enabled: bool) -> None:
        """Toggle aromatic rendering mode."""
        self.render_aromatic_as_circle = enabled
        # Canvas refresh will call update_positions with current atom positions.

    def set_ring_context(self, ring_center: QPointF | None) -> None:
        self._ring_center = ring_center

    def set_offset_sign(self, sign: int) -> None:
        self._offset_sign = 1 if sign >= 0 else -1

    def update_positions(self, atom1: Atom, atom2: Atom) -> None:
        """Redraw the bond path based on atom positions and bond type."""
        x1, y1 = atom1.x, atom1.y
        x2, y2 = atom2.x, atom2.y
        path = QPainterPath()
        color = QColor(self._style.bond_color)
        dx = x2 - x1
        dy = y2 - y1
        length = math.hypot(dx, dy)
        if length <= 1e-6:
            self.setPath(path)
            return
        nx = -dy / length  # Normal vector perpendicular to bond
        ny = dx / length
        ux = dx / length
        uy = dy / length
        desired_length = self.length_px if self.length_px is not None else self._style.bond_length_px
        if desired_length <= 0:
            desired_length = length
        render_length = desired_length
        midx = (x1 + x2) / 2
        midy = (y1 + y2) / 2
        half = desired_length / 2
        p1x = midx - ux * half
        p1y = midy - uy * half
        p2x = midx + ux * half
        p2y = midy + uy * half

        offset_sign = self._offset_sign
        if self._ring_center is not None:
            midx = (x1 + x2) / 2
            midy = (y1 + y2) / 2
            vx = self._ring_center.x() - midx
            vy = self._ring_center.y() - midy
            offset_sign = 1 if (nx * vx + ny * vy) >= 0 else -1

        # Aromatic bonds: if circle mode, draw as single line
        if self.is_aromatic and self.render_aromatic_as_circle:
            path.moveTo(p1x, p1y)
            path.lineTo(p2x, p2y)
            pen = QPen(color, self._style.stroke_px)
            pen.setCapStyle(self._style.cap_style)
            pen.setJoinStyle(self._style.join_style)
            self.setPen(pen)
            self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
            self.setPath(path)
            return

        effective_order = self.order
        if self.is_aromatic and self.display_order is not None:
            effective_order = self.display_order

        if self.style == BondStyle.PLAIN:
            if effective_order == 1:
                path.moveTo(p1x, p1y)
                path.lineTo(p2x, p2y)
            else:
                offset = self._style.double_offset_px
                if effective_order == 2:
                    path.moveTo(p1x, p1y)
                    path.lineTo(p2x, p2y)
                    q1x = p1x + nx * offset * offset_sign + ux * self._style.inner_trim_px
                    q1y = p1y + ny * offset * offset_sign + uy * self._style.inner_trim_px
                    q2x = p2x + nx * offset * offset_sign - ux * self._style.inner_trim_px
                    q2y = p2y + ny * offset * offset_sign - uy * self._style.inner_trim_px
                    path.moveTo(q1x, q1y)
                    path.lineTo(q2x, q2y)
                else:  # Triple bond
                    path.moveTo(p1x, p1y)
                    path.lineTo(p2x, p2y)
                    path.moveTo(p1x + nx * offset, p1y + ny * offset)
                    path.lineTo(p2x + nx * offset, p2y + ny * offset)
                    path.moveTo(p1x - nx * offset, p1y - ny * offset)
                    path.lineTo(p2x - nx * offset, p2y - ny * offset)
            pen = QPen(color, self._style.stroke_px)
            pen.setCapStyle(self._style.cap_style)
            pen.setJoinStyle(self._style.join_style)
            self.setPen(pen)
            self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
            
        elif self.style == BondStyle.INTERACTION:
            path.moveTo(p1x, p1y)
            path.lineTo(p2x, p2y)
            pen = QPen(color, self._style.stroke_px, Qt.PenStyle.DotLine)
            pen.setCapStyle(self._style.cap_style)
            pen.setJoinStyle(self._style.join_style)
            self.setPen(pen)
            self.setBrush(QBrush(Qt.BrushStyle.NoBrush))

        elif self.style == BondStyle.WEDGE:
            width = self._style.wedge_width_px
            b1x = p2x + nx * (width / 2)
            b1y = p2y + ny * (width / 2)
            b2x = p2x - nx * (width / 2)
            b2y = p2y - ny * (width / 2)
            path.moveTo(p1x, p1y)
            path.lineTo(b1x, b1y)
            path.lineTo(b2x, b2y)
            path.closeSubpath()
            pen = QPen(color, 1)
            pen.setCapStyle(self._style.cap_style)
            pen.setJoinStyle(self._style.join_style)
            self.setPen(pen)
            self.setBrush(QBrush(color))
            
        elif self.style == BondStyle.HASHED:
            steps = self._style.hash_count
            for i in range(1, steps + 1):
                t = i / (steps + 1)
                px = p1x + (p2x - p1x) * t
                py = p1y + (p2y - p1y) * t
                width = self._style.hash_min_px + (self._style.hash_max_px - self._style.hash_min_px) * t
                path.moveTo(px + nx * width / 2, py + ny * width / 2)
                path.lineTo(px - nx * width / 2, py - ny * width / 2)
            pen = QPen(color, self._style.hash_stroke_px)
            pen.setCapStyle(self._style.cap_style)
            pen.setJoinStyle(self._style.join_style)
            self.setPen(pen)
            self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
            
        elif self.style == BondStyle.WAVY:
            segments = max(6, int(render_length / 6))
            amplitude = 3.0
            for i in range(segments + 1):
                t = i / segments
                px = p1x + (p2x - p1x) * t
                py = p1y + (p2y - p1y) * t
                offset = math.sin(t * math.pi * 4) * amplitude
                wx = px + nx * offset
                wy = py + ny * offset
                if i == 0:
                    path.moveTo(wx, wy)
                else:
                    path.lineTo(wx, wy)
            pen = QPen(color, self._style.stroke_px)
            pen.setCapStyle(self._style.cap_style)
            pen.setJoinStyle(self._style.join_style)
            self.setPen(pen)
            self.setBrush(QBrush(Qt.BrushStyle.NoBrush))

        self.setPath(path)

    def set_style(self, style: DrawingStyle, atom1: Atom, atom2: Atom) -> None:
        """Update drawing style and redraw bond."""
        self._style = style
        self.update_positions(atom1, atom2)


class ArrowItem(QGraphicsPathItem):
    """Graphics item representing a reaction arrow."""

    def __init__(
        self,
        start: QPointF,
        end: QPointF,
        head_at_end: bool = True,
        style: DrawingStyle = CHEMDOODLE_LIKE,
    ) -> None:
        super().__init__()
        self._style = style
        self._head_at_end = head_at_end
        pen = QPen(QColor(self._style.bond_color), self._style.stroke_px)
        pen.setCapStyle(self._style.cap_style)
        pen.setJoinStyle(self._style.join_style)
        self.setPen(pen)
        self.setBrush(QBrush(QColor(self._style.bond_color)))
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.setZValue(5)
        self.update_positions(start, end)

    def update_positions(self, start: QPointF, end: QPointF) -> None:
        dx = end.x() - start.x()
        dy = end.y() - start.y()
        length = math.hypot(dx, dy)
        if length < 1e-6:
            self.setPath(QPainterPath())
            return

        ux = dx / length
        uy = dy / length
        head_len = 12.0
        head_width = 6.0

        if self._head_at_end:
            tip = end
            tail = start
            dir_x = ux
            dir_y = uy
        else:
            tip = start
            tail = end
            dir_x = -ux
            dir_y = -uy

        base_x = tip.x() - dir_x * head_len
        base_y = tip.y() - dir_y * head_len
        nx = -dir_y
        ny = dir_x
        left = QPointF(base_x + nx * head_width, base_y + ny * head_width)
        right = QPointF(base_x - nx * head_width, base_y - ny * head_width)

        path = QPainterPath()
        path.moveTo(tail)
        path.lineTo(QPointF(base_x, base_y))
        path.moveTo(left)
        path.lineTo(tip)
        path.lineTo(right)
        path.closeSubpath()
        self.setPath(path)

    def set_style(self, style: DrawingStyle) -> None:
        self._style = style
        pen = QPen(QColor(self._style.bond_color), self._style.stroke_px)
        pen.setCapStyle(self._style.cap_style)
        pen.setJoinStyle(self._style.join_style)
        self.setPen(pen)
        self.setBrush(QBrush(QColor(self._style.bond_color)))


class BracketItem(QGraphicsPathItem):
    """Graphics item representing brackets around a region."""

    def __init__(
        self,
        rect: QRectF,
        kind: str = "[]",
        padding: float = 8.0,
        style: DrawingStyle = CHEMDOODLE_LIKE,
    ) -> None:
        super().__init__()
        self._rect = rect.adjusted(-padding, -padding, padding, padding)
        self._kind = kind
        self._style = style
        pen = QPen(QColor(self._style.bond_color), self._style.stroke_px)
        pen.setCapStyle(self._style.cap_style)
        pen.setJoinStyle(self._style.join_style)
        self.setPen(pen)
        self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
        self.setZValue(2)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self._update_path()

    def _update_path(self) -> None:
        rect = self._rect
        path = QPainterPath()
        height = rect.height()
        top = rect.top()
        bottom = rect.bottom()
        left = rect.left()
        right = rect.right()

        if self._kind == "()":
            width = max(8.0, height * 0.12)
            mid = (top + bottom) / 2
            path.moveTo(left + width, top)
            path.quadTo(left, top + height * 0.25, left, mid)
            path.quadTo(left, bottom - height * 0.25, left + width, bottom)

            path.moveTo(right - width, top)
            path.quadTo(right, top + height * 0.25, right, mid)
            path.quadTo(right, bottom - height * 0.25, right - width, bottom)
        elif self._kind == "{}":
            width = max(8.0, height * 0.12)
            mid = (top + bottom) / 2
            notch = height * 0.08
            path.moveTo(left + width, top)
            path.quadTo(left, top + height * 0.2, left + width, mid - notch)
            path.quadTo(left + width * 1.4, mid, left + width, mid + notch)
            path.quadTo(left, bottom - height * 0.2, left + width, bottom)

            path.moveTo(right - width, top)
            path.quadTo(right, top + height * 0.2, right - width, mid - notch)
            path.quadTo(right - width * 1.4, mid, right - width, mid + notch)
            path.quadTo(right, bottom - height * 0.2, right - width, bottom)
        else:  # "[]"
            arm = max(6.0, height * 0.08)
            path.moveTo(left + arm, top)
            path.lineTo(left, top)
            path.lineTo(left, bottom)
            path.lineTo(left + arm, bottom)

            path.moveTo(right - arm, top)
            path.lineTo(right, top)
            path.lineTo(right, bottom)
            path.lineTo(right - arm, bottom)

        self.setPath(path)

    def set_style(self, style: DrawingStyle) -> None:
        self._style = style
        pen = QPen(QColor(self._style.bond_color), self._style.stroke_px)
        pen.setCapStyle(self._style.cap_style)
        pen.setJoinStyle(self._style.join_style)
        self.setPen(pen)

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


class HoverAtomIndicatorItem(QGraphicsEllipseItem):
    """Amber circle overlay for hovered atoms."""

    def __init__(self, radius: float = 10.0) -> None:
        super().__init__(-radius, -radius, radius * 2, radius * 2)
        self._radius = radius
        pen = QPen(QColor("#E0A825"), 1.5)
        self.setPen(pen)
        self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
        self.setZValue(50)
        self.setVisible(False)

    def update_position(self, x: float, y: float) -> None:
        self.setPos(x, y)
        if not self.isVisible():
            self.setVisible(True)

    def hide_indicator(self) -> None:
        self.setVisible(False)


class HoverBondIndicatorItem(QGraphicsPathItem):
    """Amber parentheses overlay for hovered bonds."""

    def __init__(self, radius: float = 10.0, separation: float = 12.0) -> None:
        super().__init__()
        self._radius = radius
        self._separation = separation
        pen = QPen(QColor("#E0A825"), 1.5)
        self.setPen(pen)
        self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
        self.setZValue(50)
        self.setVisible(False)

    def update_for_bond(self, p1: QPointF, p2: QPointF) -> None:
        mid = QPointF((p1.x() + p2.x()) / 2, (p1.y() + p2.y()) / 2)
        angle = math.degrees(math.atan2(p2.y() - p1.y(), p2.x() - p1.x()))

        r = self._radius
        sep = self._separation
        path = QPainterPath()

        left_rect = QRectF(-sep - r, -r, 2 * r, 2 * r)
        right_rect = QRectF(sep - r, -r, 2 * r, 2 * r)
        path.arcMoveTo(left_rect, 60)
        path.arcTo(left_rect, 60, 240)
        path.arcMoveTo(right_rect, 120)
        path.arcTo(right_rect, 120, 240)

        self.setPath(path)
        self.setPos(mid)
        self.setRotation(angle)
        if not self.isVisible():
            self.setVisible(True)

    def hide_indicator(self) -> None:
        self.setVisible(False)


class OptimizeZoneItem(QGraphicsEllipseItem):
    """Blue translucent optimize zone around sprout anchor."""

    def __init__(self, radius: float = 28.0) -> None:
        super().__init__(-radius, -radius, radius * 2, radius * 2)
        self._radius = radius
        pen = QPen(QColor("#4A90D9"), 1.2)
        brush = QBrush(QColor(74, 144, 217, 40))
        self.setPen(pen)
        self.setBrush(brush)
        self.setZValue(30)
        self.setVisible(False)

    def update_center(self, x: float, y: float) -> None:
        self.setPos(x, y)
        if not self.isVisible():
            self.setVisible(True)

    def set_radius(self, radius: float) -> None:
        self._radius = radius
        self.setRect(-radius, -radius, radius * 2, radius * 2)

    def radius(self) -> float:
        return self._radius

    def hide_zone(self) -> None:
        self.setVisible(False)


class PreviewBondItem(QGraphicsPathItem):
    """Preview line for bond placement."""

    def __init__(self) -> None:
        super().__init__()
        pen = QPen(QColor("#4A90D9"), 1.5, Qt.PenStyle.DashLine)
        self.setPen(pen)
        self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
        self.setZValue(40)
        self.setVisible(False)

    def update_line(self, p1: QPointF, p2: QPointF) -> None:
        path = QPainterPath()
        path.moveTo(p1)
        path.lineTo(p2)
        self.setPath(path)
        if not self.isVisible():
            self.setVisible(True)

    def hide_preview(self) -> None:
        self.setVisible(False)


class PreviewRingItem(QGraphicsPathItem):
    """Preview polygon for ring placement."""

    def __init__(self) -> None:
        super().__init__()
        pen = QPen(QColor("#4A90D9"), 1.5, Qt.PenStyle.DashLine)
        self.setPen(pen)
        self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
        self.setZValue(40)
        self.setVisible(False)

    def update_polygon(self, vertices: list[QPointF]) -> None:
        if not vertices:
            self.setVisible(False)
            return
        path = QPainterPath()
        path.moveTo(vertices[0])
        for v in vertices[1:]:
            path.lineTo(v)
        path.closeSubpath()
        self.setPath(path)
        if not self.isVisible():
            self.setVisible(True)

    def hide_preview(self) -> None:
        self.setVisible(False)


class PreviewChainItem(QGraphicsPathItem):
    """Preview polyline for chain placement."""

    def __init__(self) -> None:
        super().__init__()
        pen = QPen(QColor("#4A90D9"), 1.5, Qt.PenStyle.DashLine)
        self.setPen(pen)
        self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
        self.setZValue(40)
        self.setVisible(False)

    def update_polyline(self, points: list[QPointF]) -> None:
        if len(points) < 2:
            self.setVisible(False)
            return
        path = QPainterPath()
        path.moveTo(points[0])
        for p in points[1:]:
            path.lineTo(p)
        self.setPath(path)
        if not self.isVisible():
            self.setVisible(True)

    def hide_preview(self) -> None:
        self.setVisible(False)


class PreviewChainLabelItem(QGraphicsTextItem):
    """Preview label showing chain length."""

    def __init__(self) -> None:
        super().__init__()
        self.setDefaultTextColor(QColor("#4A90D9"))
        font = QFont("Arial", 10, QFont.Weight.Bold)
        self.setFont(font)
        self.setZValue(41)
        self.setVisible(False)

    def update_label(self, text: str, pos: QPointF) -> None:
        self.setPlainText(text)
        rect = self.boundingRect()
        self.setPos(pos.x() - rect.width() / 2, pos.y() - rect.height() / 2)
        if not self.isVisible():
            self.setVisible(True)

    def hide_label(self) -> None:
        self.setVisible(False)
