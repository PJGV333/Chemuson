"""
Chemuson Icon Library
Generate vector icons programmatically using QPainter for a professional look.
"""
from PyQt6.QtWidgets import QApplication
from PyQt6.QtGui import (
    QIcon,
    QPixmap,
    QPainter,
    QColor,
    QPen,
    QFont,
    QPolygonF,
    QBrush,
    QPainterPath,
)
from PyQt6.QtCore import Qt, QSize, QPointF, QRectF
import math


# Standard icon size
ICON_SIZE = 32

# Color palette for atoms
ATOM_COLORS = {
    'C': '#333333',   # Carbon - dark gray
    'N': '#3050F8',   # Nitrogen - blue
    'O': '#FF0D0D',   # Oxygen - red
    'S': '#FFFF30',   # Sulfur - yellow
    'P': '#FF8000',   # Phosphorus - orange
    'F': '#90E050',   # Fluorine - light green
    'Cl': '#1FF01F',  # Chlorine - green
    'Br': '#A62929',  # Bromine - dark red
    'H': '#FFFFFF',   # Hydrogen - white
}


def draw_atom_icon(text: str, color: str = None) -> QIcon:
    """
    Draw an atom icon with the element symbol centered.
    
    Args:
        text: Element symbol (e.g., 'C', 'N', 'Cl')
        color: Hex color string. If None, uses ATOM_COLORS or defaults to black.
    
    Returns:
        QIcon with the element symbol
    """
    if color is None:
        color = ATOM_COLORS.get(text, '#333333')
    
    pixmap = QPixmap(ICON_SIZE, ICON_SIZE)
    pixmap.fill(Qt.GlobalColor.transparent)
    
    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    painter.setRenderHint(QPainter.RenderHint.TextAntialiasing)
    
    # Draw circular background
    painter.setBrush(QBrush(QColor(color)))
    painter.setPen(QPen(QColor('#222222'), 1.5))
    margin = 2
    painter.drawEllipse(margin, margin, ICON_SIZE - 2*margin, ICON_SIZE - 2*margin)
    
    # Determine text color based on background brightness
    bg_color = QColor(color)
    brightness = (bg_color.red() * 299 + bg_color.green() * 587 + bg_color.blue() * 114) / 1000
    text_color = '#FFFFFF' if brightness < 128 else '#000000'
    
    # Draw element symbol
    font = QFont('Arial', 14 if len(text) == 1 else 11, QFont.Weight.Bold)
    painter.setFont(font)
    painter.setPen(QColor(text_color))
    painter.drawText(QRectF(0, 0, ICON_SIZE, ICON_SIZE), Qt.AlignmentFlag.AlignCenter, text)
    
    painter.end()
    return QIcon(pixmap)


def draw_bond_icon(bond_type: str = 'single') -> QIcon:
    """
    Draw a bond icon showing a diagonal line.
    
    Args:
        bond_type: 'single' or 'double'
    
    Returns:
        QIcon with bond representation
    """
    pixmap = QPixmap(ICON_SIZE, ICON_SIZE)
    pixmap.fill(Qt.GlobalColor.transparent)
    
    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    
    margin = 6
    x1, y1 = margin, ICON_SIZE - margin
    x2, y2 = ICON_SIZE - margin, margin
    
    if bond_type == 'single':
        pen = QPen(QColor('#333333'), 3)
        painter.setPen(pen)
        painter.drawLine(x1, y1, x2, y2)

    elif bond_type == 'double':
        pen = QPen(QColor('#333333'), 2)
        painter.setPen(pen)
        # Offset lines for double bond
        offset = 3
        # Calculate perpendicular offset
        dx, dy = x2 - x1, y2 - y1
        length = math.sqrt(dx*dx + dy*dy)
        nx, ny = -dy/length * offset, dx/length * offset
        
        painter.drawLine(int(x1 + nx), int(y1 + ny), int(x2 + nx), int(y2 + ny))
        painter.drawLine(int(x1 - nx), int(y1 - ny), int(x2 - nx), int(y2 - ny))

    elif bond_type == 'triple':
        pen = QPen(QColor('#333333'), 2)
        painter.setPen(pen)
        offset = 4
        dx, dy = x2 - x1, y2 - y1
        length = math.sqrt(dx * dx + dy * dy)
        nx, ny = -dy / length * offset, dx / length * offset
        painter.drawLine(int(x1), int(y1), int(x2), int(y2))
        painter.drawLine(int(x1 + nx), int(y1 + ny), int(x2 + nx), int(y2 + ny))
        painter.drawLine(int(x1 - nx), int(y1 - ny), int(x2 - nx), int(y2 - ny))

    elif bond_type == 'wedge':
        painter.setPen(QPen(QColor('#333333'), 1))
        painter.setBrush(QBrush(QColor('#333333')))
        dx, dy = x2 - x1, y2 - y1
        length = math.sqrt(dx * dx + dy * dy)
        nx, ny = -dy / length * 4, dx / length * 4
        points = [
            QPointF(x1, y1),
            QPointF(x2 + nx, y2 + ny),
            QPointF(x2 - nx, y2 - ny),
        ]
        painter.drawPolygon(QPolygonF(points))

    elif bond_type == 'hashed':
        pen = QPen(QColor('#333333'), 1.5)
        painter.setPen(pen)
        dx, dy = x2 - x1, y2 - y1
        length = math.sqrt(dx * dx + dy * dy)
        nx, ny = -dy / length, dx / length
        steps = 6
        for i in range(1, steps + 1):
            t = i / (steps + 1)
            px = x1 + dx * t
            py = y1 + dy * t
            width = 4 * t
            painter.drawLine(
                int(px + nx * width),
                int(py + ny * width),
                int(px - nx * width),
                int(py - ny * width),
            )

    elif bond_type == 'wavy':
        pen = QPen(QColor('#333333'), 1.5)
        painter.setPen(pen)
        dx, dy = x2 - x1, y2 - y1
        length = math.sqrt(dx * dx + dy * dy)
        nx, ny = -dy / length, dx / length
        path = QPainterPath()
        segments = 8
        amplitude = 2.5
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
        painter.drawPath(path)
    
    painter.end()
    return QIcon(pixmap)


def draw_ring_icon() -> QIcon:
    """
    Draw a benzene ring (hexagon) icon.
    
    Returns:
        QIcon with hexagon shape
    """
    pixmap = QPixmap(ICON_SIZE, ICON_SIZE)
    pixmap.fill(Qt.GlobalColor.transparent)
    
    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    
    # Create hexagon points
    center = ICON_SIZE / 2
    radius = ICON_SIZE / 2 - 4
    points = []
    for i in range(6):
        angle = math.pi / 6 + i * math.pi / 3  # Start from top-right vertex
        x = center + radius * math.cos(angle)
        y = center - radius * math.sin(angle)
        points.append(QPointF(x, y))
    
    polygon = QPolygonF(points)
    
    pen = QPen(QColor('#333333'), 2)
    painter.setPen(pen)
    painter.setBrush(Qt.BrushStyle.NoBrush)
    painter.drawPolygon(polygon)
    
    # Draw inner circle (aromatic indicator)
    inner_radius = radius * 0.5
    painter.drawEllipse(QPointF(center, center), inner_radius, inner_radius)
    
    painter.end()
    return QIcon(pixmap)


def draw_generic_icon(shape: str) -> QIcon:
    """
    Draw generic tool icons (pointer, eraser, etc).
    
    Args:
        shape: 'pointer', 'eraser', 'pan', 'zoom_in', 'zoom_out', 'chain'
    
    Returns:
        QIcon with the tool shape
    """
    pixmap = QPixmap(ICON_SIZE, ICON_SIZE)
    pixmap.fill(Qt.GlobalColor.transparent)
    
    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    
    if shape == 'pointer':
        # Draw arrow cursor
        arrow_points = [
            QPointF(6, 4),    # Top
            QPointF(6, 24),   # Bottom left
            QPointF(11, 19),  # Inner corner
            QPointF(17, 27),  # Arrow tip right
            QPointF(20, 24),  # Arrow base right
            QPointF(14, 16),  # Inner corner 2
            QPointF(19, 11),  # Right wing
        ]
        polygon = QPolygonF(arrow_points)
        painter.setPen(QPen(QColor('#222222'), 1.5))
        painter.setBrush(QBrush(QColor('#4A90D9')))
        painter.drawPolygon(polygon)
    
    elif shape == 'eraser':
        # Draw eraser rectangle
        painter.setPen(QPen(QColor('#333333'), 1.5))
        painter.setBrush(QBrush(QColor('#FFB6C1')))  # Light pink
        # Draw angled eraser
        eraser_points = [
            QPointF(8, 26),
            QPointF(4, 18),
            QPointF(20, 6),
            QPointF(28, 10),
            QPointF(12, 26),
        ]
        polygon = QPolygonF(eraser_points)
        painter.drawPolygon(polygon)
        # Eraser tip
        painter.setBrush(QBrush(QColor('#FFFFFF')))
        tip_points = [
            QPointF(4, 18),
            QPointF(8, 14),
            QPointF(12, 16),
            QPointF(8, 20),
        ]
        painter.drawPolygon(QPolygonF(tip_points))
    
    elif shape == 'pan':
        # Draw hand/move icon
        painter.setPen(QPen(QColor('#333333'), 2))
        center = ICON_SIZE / 2
        arrow_len = 8
        # Four arrows
        for angle in [0, 90, 180, 270]:
            rad = math.radians(angle)
            x1 = center + 4 * math.cos(rad)
            y1 = center - 4 * math.sin(rad)
            x2 = center + arrow_len * math.cos(rad)
            y2 = center - arrow_len * math.sin(rad)
            painter.drawLine(int(x1), int(y1), int(x2), int(y2))
            # Arrow heads
            head_len = 3
            for head_angle in [angle + 30, angle - 30]:
                hrad = math.radians(head_angle + 180)
                hx = x2 + head_len * math.cos(hrad)
                hy = y2 - head_len * math.sin(hrad)
                painter.drawLine(int(x2), int(y2), int(hx), int(hy))
    
    elif shape == 'zoom_in':
        # Draw magnifying glass with plus
        painter.setPen(QPen(QColor('#333333'), 2))
        painter.setBrush(QBrush(QColor('#E0E0E0')))
        painter.drawEllipse(4, 4, 18, 18)
        # Handle
        painter.drawLine(19, 19, 27, 27)
        # Plus sign
        painter.drawLine(10, 13, 16, 13)
        painter.drawLine(13, 10, 13, 16)
    
    elif shape == 'zoom_out':
        # Draw magnifying glass with minus
        painter.setPen(QPen(QColor('#333333'), 2))
        painter.setBrush(QBrush(QColor('#E0E0E0')))
        painter.drawEllipse(4, 4, 18, 18)
        # Handle
        painter.drawLine(19, 19, 27, 27)
        # Minus sign
        painter.drawLine(10, 13, 16, 13)

    elif shape == 'chain':
        pen = QPen(QColor('#333333'), 2)
        painter.setPen(pen)
        points = [
            QPointF(6, 24),
            QPointF(12, 18),
            QPointF(18, 24),
            QPointF(24, 18),
            QPointF(28, 24),
        ]
        for i in range(len(points) - 1):
            painter.drawLine(points[i], points[i + 1])
    
    painter.end()
    return QIcon(pixmap)


# Convenience functions for common icons
def get_pointer_icon() -> QIcon:
    return draw_generic_icon('pointer')

def get_eraser_icon() -> QIcon:
    return draw_generic_icon('eraser')

def get_single_bond_icon() -> QIcon:
    return draw_bond_icon('single')

def get_double_bond_icon() -> QIcon:
    return draw_bond_icon('double')

def get_benzene_icon() -> QIcon:
    return draw_ring_icon()
