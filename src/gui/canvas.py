"""
Chemuson Canvas
Page-based canvas using QGraphicsView/QGraphicsScene for document-style editing.
"""
from __future__ import annotations

from typing import Optional

from PyQt6.QtWidgets import (
    QGraphicsView,
    QGraphicsScene,
    QGraphicsRectItem,
    QGraphicsDropShadowEffect,
)
from PyQt6.QtGui import QPainter, QPen, QColor, QBrush, QWheelEvent, QUndoStack
from PyQt6.QtCore import Qt

from core.model import ChemState, MolGraph
from gui.items import AtomItem, BondItem
from gui.commands import AddAtomCommand, AddBondCommand


# Paper dimensions (A4, approximately)
PAPER_WIDTH = 800
PAPER_HEIGHT = 1000
PAPER_MARGIN = 40  # Margins where atoms shouldn't be placed
ATOM_HIT_RADIUS = 20


class ChemusonCanvas(QGraphicsView):
    """
    Page-based canvas for drawing molecules.
    Uses QGraphicsView/QGraphicsScene with a centered paper sheet.
    """

    def __init__(self, parent=None) -> None:
        super().__init__(parent)

        self.scene = QGraphicsScene()
        self.setScene(self.scene)

        self.model = MolGraph()
        self.state = ChemState()
        self.undo_stack = QUndoStack(self)

        self.atom_items: dict[int, AtomItem] = {}
        self.bond_items: dict[int, BondItem] = {}

        self.current_tool: str = self.state.active_tool
        self.selected_atom_id: Optional[int] = None
        self.hover_atom_id: Optional[int] = None

        self._setup_view()
        self._create_paper()

        self.setMouseTracking(True)
        self.viewport().setMouseTracking(True)

    def _setup_view(self) -> None:
        self.setBackgroundBrush(QBrush(QColor("#E0E0E0")))
        self.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.setRenderHint(QPainter.RenderHint.TextAntialiasing)

        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)

        margin = 100
        self.scene.setSceneRect(
            -margin,
            -margin,
            PAPER_WIDTH + 2 * margin,
            PAPER_HEIGHT + 2 * margin,
        )

        self.setDragMode(QGraphicsView.DragMode.NoDrag)

        self._zoom_factor = 1.0
        self._min_zoom = 0.25
        self._max_zoom = 4.0

    def _create_paper(self) -> None:
        self.paper = QGraphicsRectItem(0, 0, PAPER_WIDTH, PAPER_HEIGHT)
        self.paper.setBrush(QBrush(Qt.GlobalColor.white))
        self.paper.setPen(QPen(QColor("#CCCCCC"), 1))
        self.paper.setZValue(-10)

        shadow = QGraphicsDropShadowEffect()
        shadow.setBlurRadius(20)
        shadow.setColor(QColor(0, 0, 0, 80))
        shadow.setOffset(5, 5)
        self.paper.setGraphicsEffect(shadow)

        self.scene.addItem(self.paper)
        self.centerOn(PAPER_WIDTH / 2, PAPER_HEIGHT / 2)

    def set_current_tool(self, tool_id: str) -> None:
        self.current_tool = tool_id
        self.state.active_tool = tool_id
        self._clear_selected_atom()
        if tool_id.startswith("atom_"):
            self.state.default_element = tool_id.split("_", 1)[1]

        if tool_id == "tool_select":
            self.setCursor(Qt.CursorShape.ArrowCursor)
            self.setDragMode(QGraphicsView.DragMode.RubberBandDrag)
        else:
            self.setCursor(Qt.CursorShape.CrossCursor)
            self.setDragMode(QGraphicsView.DragMode.NoDrag)

    def mousePressEvent(self, event) -> None:
        scene_pos = self.mapToScene(event.pos())

        if not self._is_on_paper(scene_pos.x(), scene_pos.y()):
            super().mousePressEvent(event)
            return

        if event.button() != Qt.MouseButton.LeftButton:
            super().mousePressEvent(event)
            return

        clicked_atom_id = self._get_atom_at(scene_pos.x(), scene_pos.y())

        if self.current_tool.startswith("atom_"):
            element = self.current_tool.split("_")[1]
            if clicked_atom_id is None:
                cmd = AddAtomCommand(self.model, self, element, scene_pos.x(), scene_pos.y())
                self.undo_stack.push(cmd)
            else:
                self._set_selected_atom(clicked_atom_id)
            return

        if self.current_tool.startswith("bond_"):
            order = 1 if "single" in self.current_tool else 2
            if clicked_atom_id is None:
                if self.selected_atom_id is not None:
                    self.undo_stack.beginMacro("Add atom and bond")
                    atom_cmd = AddAtomCommand(
                        self.model,
                        self,
                        self.state.default_element,
                        scene_pos.x(),
                        scene_pos.y(),
                    )
                    self.undo_stack.push(atom_cmd)
                    new_atom_id = atom_cmd.atom_id
                    if new_atom_id is not None:
                        bond_cmd = AddBondCommand(
                            self.model,
                            self,
                            self.selected_atom_id,
                            new_atom_id,
                            order,
                        )
                        self.undo_stack.push(bond_cmd)
                    self.undo_stack.endMacro()
                    self._clear_selected_atom()
                return

            if self.selected_atom_id is None:
                self._set_selected_atom(clicked_atom_id)
                return

            if self.selected_atom_id == clicked_atom_id:
                self._clear_selected_atom()
                return

            if self.model.find_bond_between(self.selected_atom_id, clicked_atom_id) is None:
                cmd = AddBondCommand(
                    self.model,
                    self,
                    self.selected_atom_id,
                    clicked_atom_id,
                    order,
                )
                self.undo_stack.push(cmd)
            self._clear_selected_atom()
            return

        if self.current_tool == "tool_select":
            super().mousePressEvent(event)
            return

        super().mousePressEvent(event)

    def mouseMoveEvent(self, event) -> None:
        scene_pos = self.mapToScene(event.pos())
        hover_atom_id = self._get_atom_at(scene_pos.x(), scene_pos.y())

        if self.hover_atom_id is not None and self.hover_atom_id in self.atom_items:
            if self.hover_atom_id != self.selected_atom_id:
                self.atom_items[self.hover_atom_id].set_hover(False)
            self.hover_atom_id = None

        if hover_atom_id is not None and hover_atom_id in self.atom_items:
            if hover_atom_id != self.selected_atom_id:
                self.atom_items[hover_atom_id].set_hover(True)
            self.hover_atom_id = hover_atom_id

        super().mouseMoveEvent(event)

    def wheelEvent(self, event: QWheelEvent) -> None:
        if event.angleDelta().y() > 0:
            self.zoom_in()
        else:
            self.zoom_out()

    def zoom_in(self) -> None:
        if self._zoom_factor < self._max_zoom:
            self._zoom_factor *= 1.2
            self.scale(1.2, 1.2)

    def zoom_out(self) -> None:
        if self._zoom_factor > self._min_zoom:
            self._zoom_factor /= 1.2
            self.scale(1 / 1.2, 1 / 1.2)

    def clear_canvas(self) -> None:
        self.scene.clear()
        self.model.clear()
        self.undo_stack.clear()
        self.atom_items.clear()
        self.bond_items.clear()
        self.selected_atom_id = None
        self.hover_atom_id = None
        self._create_paper()

    def add_atom_item(self, atom) -> None:
        if atom.id in self.atom_items:
            return
        item = AtomItem(atom)
        self.scene.addItem(item)
        self.atom_items[atom.id] = item

    def remove_atom_item(self, atom_id: int) -> None:
        item = self.atom_items.pop(atom_id, None)
        if item is not None:
            self.scene.removeItem(item)
        if self.selected_atom_id == atom_id:
            self.selected_atom_id = None
        if self.hover_atom_id == atom_id:
            self.hover_atom_id = None

    def add_bond_item(self, bond) -> None:
        if bond.id in self.bond_items:
            return
        atom1 = self.model.get_atom(bond.a1_id)
        atom2 = self.model.get_atom(bond.a2_id)
        item = BondItem(bond, atom1, atom2)
        self.scene.addItem(item)
        self.bond_items[bond.id] = item

    def remove_bond_item(self, bond_id: int) -> None:
        item = self.bond_items.pop(bond_id, None)
        if item is not None:
            self.scene.removeItem(item)

    def _get_atom_at(self, x: float, y: float) -> Optional[int]:
        for atom in self.model.atoms.values():
            dx = atom.x - x
            dy = atom.y - y
            if (dx * dx + dy * dy) < (ATOM_HIT_RADIUS * ATOM_HIT_RADIUS):
                return atom.id
        return None

    def _set_selected_atom(self, atom_id: int) -> None:
        if self.selected_atom_id is not None and self.selected_atom_id in self.atom_items:
            self.atom_items[self.selected_atom_id].set_selected(False)
        self.selected_atom_id = atom_id
        if atom_id in self.atom_items:
            self.atom_items[atom_id].set_selected(True)

    def _clear_selected_atom(self) -> None:
        if self.selected_atom_id is not None and self.selected_atom_id in self.atom_items:
            self.atom_items[self.selected_atom_id].set_selected(False)
        self.selected_atom_id = None

    def _is_on_paper(self, x: float, y: float) -> bool:
        return (
            PAPER_MARGIN <= x <= PAPER_WIDTH - PAPER_MARGIN
            and PAPER_MARGIN <= y <= PAPER_HEIGHT - PAPER_MARGIN
        )
