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

from core.model import Atom, Bond


class AtomItem(QGraphicsEllipseItem):
    def __init__(self, atom: Atom, radius: float = 12.0) -> None:
        super().__init__(-radius, -radius, radius * 2, radius * 2)
        self.atom_id = atom.id
        self.element = atom.element
        self.setPos(atom.x, atom.y)

        self.setBrush(QBrush(Qt.GlobalColor.white))
        self.setPen(QPen(QColor("#333333"), 1.5))
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)

        self.label = QGraphicsTextItem(atom.element, self)
        font = QFont("Arial", 10, QFont.Weight.Bold)
        self.label.setFont(font)
        rect = self.label.boundingRect()
        self.label.setPos(-rect.width() / 2, -rect.height() / 2)
        self.label.setDefaultTextColor(QColor("#333333"))

    def set_selected(self, selected: bool) -> None:
        if selected:
            self.setBrush(QBrush(QColor("#C8DFFF")))
            self.setPen(QPen(QColor("#4477AA"), 2))
        else:
            self.setBrush(QBrush(Qt.GlobalColor.white))
            self.setPen(QPen(QColor("#333333"), 1.5))

    def set_hover(self, hover: bool) -> None:
        if hover:
            self.setBrush(QBrush(QColor("#F0F0F0")))
        else:
            self.setBrush(QBrush(Qt.GlobalColor.white))

    def set_element(self, element: str) -> None:
        self.element = element
        self.label.setPlainText(element)
        rect = self.label.boundingRect()
        self.label.setPos(-rect.width() / 2, -rect.height() / 2)


class BondItem(QGraphicsPathItem):
    def __init__(self, bond: Bond, atom1: Atom, atom2: Atom) -> None:
        super().__init__()
        self.bond_id = bond.id
        self.a1_id = bond.a1_id
        self.a2_id = bond.a2_id
        self.order = bond.order
        self.setZValue(-5)
        self.setPen(QPen(QColor("#333333"), 2))
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.update_positions(atom1, atom2)

    def update_positions(self, atom1: Atom, atom2: Atom) -> None:
        x1, y1 = atom1.x, atom1.y
        x2, y2 = atom2.x, atom2.y
        path = QPainterPath()

        if self.order == 1:
            path.moveTo(x1, y1)
            path.lineTo(x2, y2)
        else:
            dx = x2 - x1
            dy = y2 - y1
            length = math.hypot(dx, dy) or 1.0
            nx = -dy / length
            ny = dx / length
            offset = 3.0

            if self.order == 2:
                path.moveTo(x1 + nx * offset, y1 + ny * offset)
                path.lineTo(x2 + nx * offset, y2 + ny * offset)
                path.moveTo(x1 - nx * offset, y1 - ny * offset)
                path.lineTo(x2 - nx * offset, y2 - ny * offset)
            else:
                path.moveTo(x1, y1)
                path.lineTo(x2, y2)
                path.moveTo(x1 + nx * offset, y1 + ny * offset)
                path.lineTo(x2 + nx * offset, y2 + ny * offset)
                path.moveTo(x1 - nx * offset, y1 - ny * offset)
                path.lineTo(x2 - nx * offset, y2 - ny * offset)

        self.setPath(path)
