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


def _build_wavy_icon_path(
    start: QPointF,
    end: QPointF,
    cycles: int = 4,
    amplitude_ratio: float = 0.28,
) -> QPainterPath:
    path = QPainterPath()
    dx = end.x() - start.x()
    dy = end.y() - start.y()
    length = math.hypot(dx, dy)
    if length <= 1e-6:
        return path

    cycles = max(1, int(cycles))
    wavelength = length / cycles
    amplitude = wavelength * amplitude_ratio

    nx = -dy / length
    ny = dx / length
    segments = max(24, cycles * 12)
    for i in range(segments + 1):
        t = i / segments
        offset = math.sin(t * 2.0 * math.pi * cycles) * amplitude
        px = start.x() + dx * t + nx * offset
        py = start.y() + dy * t + ny * offset
        if i == 0:
            path.moveTo(px, py)
        else:
            path.lineTo(px, py)
    return path


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


def draw_glyph_icon(text: str, color: str = None) -> QIcon:
    """
    Draw a minimal glyph icon (letters, brackets, symbols).
    """
    if color is None:
        color = "#222222"

    pixmap = QPixmap(ICON_SIZE, ICON_SIZE)
    pixmap.fill(Qt.GlobalColor.transparent)

    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    painter.setRenderHint(QPainter.RenderHint.TextAntialiasing)

    font_size = 15 if len(text) == 1 else 11
    font = QFont("Arial", font_size, QFont.Weight.Bold)
    painter.setFont(font)
    painter.setPen(QColor(color))
    painter.drawText(QRectF(0, 0, ICON_SIZE, ICON_SIZE), Qt.AlignmentFlag.AlignCenter, text)

    painter.end()
    return QIcon(pixmap)


def draw_charge_icon(sign: str) -> QIcon:
    """
    Draw a circled charge icon with + or -.
    """
    pixmap = QPixmap(ICON_SIZE, ICON_SIZE)
    pixmap.fill(Qt.GlobalColor.transparent)

    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)

    center = ICON_SIZE / 2
    radius = ICON_SIZE / 2 - 6

    painter.setPen(QPen(QColor("#222222"), 2))
    painter.setBrush(Qt.BrushStyle.NoBrush)
    painter.drawEllipse(QPointF(center, center), radius, radius)

    line_len = 6
    painter.drawLine(
        QPointF(center - line_len, center),
        QPointF(center + line_len, center),
    )
    if sign == "+":
        painter.drawLine(
            QPointF(center, center - line_len),
            QPointF(center, center + line_len),
        )

    painter.end()
    return QIcon(pixmap)


def draw_electron_icon(count: int = 1, spread: float = 6.0) -> QIcon:
    """
    Draw electron dot icons (single, pair, etc).
    """
    count = max(1, int(count))
    pixmap = QPixmap(ICON_SIZE, ICON_SIZE)
    pixmap.fill(Qt.GlobalColor.transparent)

    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)

    painter.setPen(QPen(Qt.PenStyle.NoPen))
    painter.setBrush(QBrush(QColor("#222222")))

    center_x = ICON_SIZE / 2
    center_y = ICON_SIZE / 2
    radius = 2.2
    total_width = (count - 1) * spread
    start_x = center_x - total_width / 2

    for i in range(count):
        x = start_x + i * spread
        painter.drawEllipse(QPointF(x, center_y), radius, radius)

    painter.end()
    return QIcon(pixmap)


def draw_radical_charge_icon(sign: str) -> QIcon:
    """
    Draw a radical dot with a small charge sign.
    """
    pixmap = QPixmap(ICON_SIZE, ICON_SIZE)
    pixmap.fill(Qt.GlobalColor.transparent)

    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)

    # Radical dot (left-center)
    painter.setPen(QPen(Qt.PenStyle.NoPen))
    painter.setBrush(QBrush(QColor("#222222")))
    dot_pos = QPointF(ICON_SIZE * 0.38, ICON_SIZE * 0.58)
    painter.drawEllipse(dot_pos, 2.2, 2.2)

    # Charge sign (top-right)
    painter.setPen(QPen(QColor("#222222"), 2))
    sign_center = QPointF(ICON_SIZE * 0.65, ICON_SIZE * 0.40)
    line_len = 4.0
    painter.drawLine(
        QPointF(sign_center.x() - line_len, sign_center.y()),
        QPointF(sign_center.x() + line_len, sign_center.y()),
    )
    if sign == "+":
        painter.drawLine(
            QPointF(sign_center.x(), sign_center.y() - line_len),
            QPointF(sign_center.x(), sign_center.y() + line_len),
        )

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
    elif bond_type == 'aromatic':
        pen = QPen(QColor('#333333'), 2)
        painter.setPen(pen)
        # One full line plus a shorter inner line
        painter.drawLine(int(x1), int(y1), int(x2), int(y2))
        trim = 0.25
        sx = x1 + (x2 - x1) * trim
        sy = y1 + (y2 - y1) * trim
        ex = x2 - (x2 - x1) * trim
        ey = y2 - (y2 - y1) * trim
        offset = 3
        dx, dy = x2 - x1, y2 - y1
        length = math.sqrt(dx * dx + dy * dy)
        nx, ny = -dy / length * offset, dx / length * offset
        painter.drawLine(int(sx + nx), int(sy + ny), int(ex + nx), int(ey + ny))
    elif bond_type == 'interaction':
        pen = QPen(QColor('#333333'), 2, Qt.PenStyle.DotLine)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        painter.setPen(pen)
        painter.drawLine(int(x1), int(y1), int(x2), int(y2))

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
        path = _build_wavy_icon_path(QPointF(x1, y1), QPointF(x2, y2))
        painter.drawPath(path)
    
    painter.end()
    return QIcon(pixmap)


def draw_wavy_anchor_icon() -> QIcon:
    """Draw a wavy anchor stub icon."""
    pixmap = QPixmap(ICON_SIZE, ICON_SIZE)
    pixmap.fill(Qt.GlobalColor.transparent)

    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    pen = QPen(QColor('#333333'), 2)
    pen.setCapStyle(Qt.PenCapStyle.RoundCap)
    pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    painter.setPen(pen)

    x = ICON_SIZE * 0.55
    y1 = ICON_SIZE * 0.2
    y2 = ICON_SIZE * 0.8
    path = _build_wavy_icon_path(QPointF(x, y1), QPointF(x, y2), cycles=4, amplitude_ratio=0.28)
    painter.drawPath(path)

    painter.end()
    return QIcon(pixmap)


def draw_ring_icon(size: int = 6, aromatic: bool = True) -> QIcon:
    """
    Draw a ring icon with the given number of sides.
    
    Returns:
        QIcon with polygon shape
    """
    pixmap = QPixmap(ICON_SIZE, ICON_SIZE)
    pixmap.fill(Qt.GlobalColor.transparent)
    
    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    
    # Create polygon points
    center = ICON_SIZE / 2
    radius = ICON_SIZE / 2 - 5
    points = []
    sides = max(3, int(size))
    step = 2 * math.pi / sides
    start_angle = math.pi / 6  # Point-up orientation (ChemDraw-like)
    for i in range(sides):
        angle = start_angle + i * step
        x = center + radius * math.cos(angle)
        y = center - radius * math.sin(angle)
        points.append(QPointF(x, y))
    
    polygon = QPolygonF(points)
    
    pen = QPen(QColor('#333333'), 2)
    painter.setPen(pen)
    painter.setBrush(Qt.BrushStyle.NoBrush)
    painter.drawPolygon(polygon)
    
    # Draw inner circle (aromatic indicator)
    if aromatic:
        inner_radius = radius * 0.5
        painter.drawEllipse(QPointF(center, center), inner_radius, inner_radius)
    
    painter.end()
    return QIcon(pixmap)


def draw_arrow_icon(kind: str = "forward") -> QIcon:
    """
    Draw arrow icons used in the annotation palette.
    """
    pixmap = QPixmap(ICON_SIZE, ICON_SIZE)
    pixmap.fill(Qt.GlobalColor.transparent)

    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)

    open_head_kinds = {
        "forward_open",
        "retro_open",
        "both_open",
        "equilibrium_open",
        "retrosynthetic",
    }
    dashed_kinds = {
        "forward_dashed",
        "retro_dashed",
        "both_dashed",
        "equilibrium_dashed",
    }
    curved_kinds = {"curved", "curved_fishhook"}

    is_open = kind in open_head_kinds
    is_dashed = kind in dashed_kinds
    is_curved = kind in curved_kinds
    is_fishhook = kind == "curved_fishhook"

    pen = QPen(QColor("#222222"), 2)
    if is_dashed:
        pen.setStyle(Qt.PenStyle.DashLine)
    painter.setPen(pen)
    painter.setBrush(QBrush(Qt.BrushStyle.NoBrush if (is_open or is_fishhook) else QColor("#222222")))

    y = ICON_SIZE / 2
    start_x = 6
    end_x = ICON_SIZE - 6
    head_len = 6
    head_width = 4

    def draw_head(tip_x: float, tip_y: float, direction: int, head_style: str) -> None:
        base_x = tip_x - direction * head_len
        left = QPointF(base_x, tip_y - head_width)
        right = QPointF(base_x, tip_y + head_width)
        tip = QPointF(tip_x, tip_y)
        if head_style == "filled":
            painter.drawPolygon(QPolygonF([tip, left, right]))
        elif head_style == "open":
            painter.drawLine(left, tip)
            painter.drawLine(tip, right)
        elif head_style == "half":
            painter.drawLine(left, tip)

    head_style = "half" if is_fishhook else ("open" if is_open else "filled")

    if kind in {"forward", "forward_open", "forward_dashed"}:
        painter.drawLine(QPointF(start_x, y), QPointF(end_x - head_len, y))
        draw_head(end_x, y, 1, head_style)
    elif kind in {"retro", "retro_open", "retro_dashed"}:
        painter.drawLine(QPointF(start_x + head_len, y), QPointF(end_x, y))
        draw_head(start_x, y, -1, head_style)
    elif kind in {"both", "both_open", "both_dashed"}:
        painter.drawLine(QPointF(start_x + head_len, y), QPointF(end_x - head_len, y))
        draw_head(end_x, y, 1, head_style)
        draw_head(start_x, y, -1, head_style)
    elif kind in {"equilibrium", "equilibrium_dashed"}:
        offset = 4
        y_top = y - offset
        y_bottom = y + offset
        painter.drawLine(QPointF(start_x, y_top), QPointF(end_x - head_len, y_top))
        draw_head(end_x, y_top, 1, head_style)
        painter.drawLine(QPointF(start_x + head_len, y_bottom), QPointF(end_x, y_bottom))
        draw_head(start_x, y_bottom, -1, head_style)
    elif kind == "retrosynthetic":
        offset = 3
        painter.drawLine(
            QPointF(start_x, y - offset),
            QPointF(end_x - head_len, y - offset),
        )
        painter.drawLine(
            QPointF(start_x, y + offset),
            QPointF(end_x - head_len, y + offset),
        )
        draw_head(end_x, y, 1, "open")
    elif is_curved:
        control = QPointF((start_x + end_x) * 0.5, y - 6)
        path = QPainterPath()
        path.moveTo(QPointF(start_x, y))
        path.quadTo(control, QPointF(end_x - head_len, y))
        painter.drawPath(path)
        draw_head(end_x, y, 1, head_style)
    else:
        painter.drawLine(QPointF(start_x, y), QPointF(end_x, y))

    painter.end()
    return QIcon(pixmap)


def draw_generic_icon(shape: str) -> QIcon:
    """
    Draw generic tool icons (pointer, eraser, etc).
    
    Args:
        shape: 'pointer', 'eraser', 'pan', 'zoom_in', 'zoom_out', 'chain', 'lasso',
            'rotate_left', 'rotate_right', 'flip_horizontal', 'flip_vertical'
    
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
        painter.setBrush(QBrush(QColor('#F2F2F2')))
        painter.drawPolygon(polygon)
    
    elif shape == 'eraser':
        # Draw eraser rectangle
        painter.setPen(QPen(QColor('#333333'), 1.5))
        painter.setBrush(QBrush(QColor('#E0E0E0')))
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
    
    elif shape in {'rotate_left', 'rotate_right'}:
        painter.setPen(QPen(QColor('#333333'), 2))
        rect = QRectF(6, 6, 20, 20)
        clockwise = shape == 'rotate_right'
        start_angle = 45 if clockwise else 135
        span = -270 if clockwise else 270
        path = QPainterPath()
        path.arcMoveTo(rect, start_angle)
        path.arcTo(rect, start_angle, span)
        painter.drawPath(path)

        end_angle = start_angle + span
        end_rad = math.radians(end_angle)
        cx = rect.center().x()
        cy = rect.center().y()
        rx = rect.width() / 2
        ry = rect.height() / 2
        end_x = cx + rx * math.cos(end_rad)
        end_y = cy - ry * math.sin(end_rad)
        tangent_angle = end_rad + (-math.pi / 2 if clockwise else math.pi / 2)
        head_len = 5
        left = QPointF(
            end_x - head_len * math.cos(tangent_angle - 0.5),
            end_y + head_len * math.sin(tangent_angle - 0.5),
        )
        right = QPointF(
            end_x - head_len * math.cos(tangent_angle + 0.5),
            end_y + head_len * math.sin(tangent_angle + 0.5),
        )
        painter.drawLine(QPointF(end_x, end_y), left)
        painter.drawLine(QPointF(end_x, end_y), right)

    elif shape == 'flip_horizontal':
        painter.setPen(QPen(QColor('#333333'), 2))
        center = ICON_SIZE / 2
        painter.drawLine(QPointF(center, 6), QPointF(center, ICON_SIZE - 6))
        painter.drawLine(QPointF(6, center), QPointF(center - 2, center))
        painter.drawLine(QPointF(ICON_SIZE - 6, center), QPointF(center + 2, center))
        painter.drawLine(QPointF(6, center), QPointF(10, center - 3))
        painter.drawLine(QPointF(6, center), QPointF(10, center + 3))
        painter.drawLine(QPointF(ICON_SIZE - 6, center), QPointF(ICON_SIZE - 10, center - 3))
        painter.drawLine(QPointF(ICON_SIZE - 6, center), QPointF(ICON_SIZE - 10, center + 3))

    elif shape == 'flip_vertical':
        painter.setPen(QPen(QColor('#333333'), 2))
        center = ICON_SIZE / 2
        painter.drawLine(QPointF(6, center), QPointF(ICON_SIZE - 6, center))
        painter.drawLine(QPointF(center, 6), QPointF(center, center - 2))
        painter.drawLine(QPointF(center, ICON_SIZE - 6), QPointF(center, center + 2))
        painter.drawLine(QPointF(center, 6), QPointF(center - 3, 10))
        painter.drawLine(QPointF(center, 6), QPointF(center + 3, 10))
        painter.drawLine(QPointF(center, ICON_SIZE - 6), QPointF(center - 3, ICON_SIZE - 10))
        painter.drawLine(QPointF(center, ICON_SIZE - 6), QPointF(center + 3, ICON_SIZE - 10))

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

    elif shape == 'lasso':
        pen = QPen(QColor('#333333'), 1.6, Qt.PenStyle.DashLine)
        painter.setPen(pen)
        path = QPainterPath()
        path.addEllipse(QRectF(6, 8, 14, 10))
        painter.drawPath(path)
        painter.setPen(QPen(QColor('#333333'), 1.8))
        painter.drawLine(16, 18, 24, 26)
    
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
