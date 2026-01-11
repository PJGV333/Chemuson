"""
Chemuson Canvas
Page-based canvas using QGraphicsView/QGraphicsScene for document-style editing.
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

from PyQt6.QtWidgets import (
    QApplication,
    QGraphicsView,
    QGraphicsScene,
    QGraphicsRectItem,
    QGraphicsDropShadowEffect,
    QGraphicsPathItem,
    QGraphicsPixmapItem,
)
from PyQt6.QtGui import (
    QPainter,
    QPainterPath,
    QPen,
    QColor,
    QBrush,
    QWheelEvent,
    QUndoStack,
    QImage,
    QPixmap,
)
from PyQt6.QtCore import Qt, QPointF, QRectF, QBuffer, QMimeData

from core.model import BondStyle, BondStereo, ChemState, MolGraph
from gui.items import AtomItem, BondItem
from gui.commands import (
    AddAtomCommand,
    AddBondCommand,
    AddRingCommand,
    ChangeAtomCommand,
    ChangeBondCommand,
    DeleteSelectionCommand,
    MoveAtomsCommand,
)
from chemio.rdkit_io import (
    molgraph_to_molfile,
    molgraph_to_smiles,
    molgraph_to_svg,
    molfile_to_molgraph,
    smiles_to_molgraph,
)


# Paper dimensions (A4, approximately)
PAPER_WIDTH = 800
PAPER_HEIGHT = 1000
PAPER_MARGIN = 40  # Margins where atoms shouldn't be placed
ATOM_HIT_RADIUS = 20
DEFAULT_BOND_LENGTH = 40.0


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
        self.bond_anchor_id: Optional[int] = None
        self.hover_atom_id: Optional[int] = None
        self._bond_last_angle: Optional[float] = None
        self._bond_zigzag_sign = 1

        self._dragging_selection = False
        self._drag_start_pos: Optional[QPointF] = None
        self._drag_start_positions: Dict[int, Tuple[float, float]] = {}
        self._drag_has_moved = False

        self._ring_dragging = False
        self._ring_anchor = None
        self._ring_preview_item: Optional[QGraphicsPathItem] = None
        self._ring_last_vertices: Optional[List[QPointF]] = None

        self._setup_view()
        self._create_paper()

        self.scene.selectionChanged.connect(self._sync_selection_from_scene)
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
        if tool_id.startswith("atom_"):
            self.state.default_element = tool_id.split("_", 1)[1]
            tool_id = "tool_atom"
        elif tool_id.startswith("bond_"):
            tool_id = "tool_bond"

        self.current_tool = tool_id
        self.state.active_tool = tool_id
        self._clear_bond_anchor()

        if tool_id == "tool_select":
            self.setCursor(Qt.CursorShape.ArrowCursor)
            self.setDragMode(QGraphicsView.DragMode.RubberBandDrag)
        else:
            self.setCursor(Qt.CursorShape.CrossCursor)
            self.setDragMode(QGraphicsView.DragMode.NoDrag)

    def set_active_bond(self, bond_spec: dict) -> None:
        self.state.active_bond_order = bond_spec.get("order", 1)
        self.state.active_bond_style = bond_spec.get("style", BondStyle.PLAIN)
        self.state.active_bond_stereo = bond_spec.get("stereo", BondStereo.NONE)
        self.state.active_bond_mode = bond_spec.get("mode", "increment")
        self.set_current_tool("tool_bond")

    def set_active_ring(self, ring_spec: dict) -> None:
        self.state.active_ring_size = ring_spec.get("size", 6)
        self.state.active_ring_aromatic = ring_spec.get("aromatic", False)
        self.set_current_tool("tool_ring")

    def set_active_element(self, element: str) -> None:
        self.state.default_element = element
        self.set_current_tool("tool_atom")

    def mousePressEvent(self, event) -> None:
        scene_pos = self.mapToScene(event.pos())

        if not self._is_on_paper(scene_pos.x(), scene_pos.y()):
            super().mousePressEvent(event)
            return

        if event.button() != Qt.MouseButton.LeftButton:
            super().mousePressEvent(event)
            return

        if self.current_tool == "tool_select":
            clicked_item = self._get_item_at(scene_pos)
            if clicked_item is None:
                if not (event.modifiers() & Qt.KeyboardModifier.ShiftModifier):
                    self.scene.clearSelection()
                super().mousePressEvent(event)
                return
            if event.modifiers() & Qt.KeyboardModifier.ShiftModifier:
                clicked_item.setSelected(not clicked_item.isSelected())
            else:
                if not clicked_item.isSelected():
                    self.scene.clearSelection()
                    clicked_item.setSelected(True)
            self._sync_selection_from_scene()
            if isinstance(clicked_item, AtomItem) and clicked_item.isSelected():
                self._begin_drag(scene_pos)
            return

        clicked_atom_id = self._get_atom_at(scene_pos.x(), scene_pos.y())
        clicked_bond_id = self._get_bond_at(scene_pos)

        if self.current_tool == "tool_atom":
            element = self.state.default_element
            if clicked_atom_id is None:
                cmd = AddAtomCommand(self.model, self, element, scene_pos.x(), scene_pos.y())
                self.undo_stack.push(cmd)
            else:
                cmd = ChangeAtomCommand(self.model, self, clicked_atom_id, element)
                self.undo_stack.push(cmd)
            return

        if self.current_tool == "tool_bond":
            if clicked_bond_id is not None:
                if self.state.active_bond_mode == "increment":
                    self._cycle_bond_order(clicked_bond_id)
                else:
                    self._apply_bond_style(clicked_bond_id)
                return

            if clicked_atom_id is None:
                if self.bond_anchor_id is None:
                    self._create_first_bond(scene_pos, event.modifiers())
                else:
                    self._extend_bond_from_anchor(scene_pos, event.modifiers())
                return

        if self.bond_anchor_id is None:
            self._set_bond_anchor(clicked_atom_id, reset_angle=True)
            return

            if self.bond_anchor_id == clicked_atom_id:
                self._clear_bond_anchor()
                return

            if self.model.find_bond_between(self.bond_anchor_id, clicked_atom_id) is None:
                self._add_bond_between(self.bond_anchor_id, clicked_atom_id)
                self._record_bond_angle_between(self.bond_anchor_id, clicked_atom_id)
            self._set_bond_anchor(clicked_atom_id, reset_angle=False)
            return

        if self.current_tool == "tool_ring":
            if clicked_bond_id is not None:
                self._start_ring_drag("bond", clicked_bond_id, scene_pos, event.modifiers())
                return
            if clicked_atom_id is not None:
                self._start_ring_drag("atom", clicked_atom_id, scene_pos, event.modifiers())
                return
            self._start_ring_drag("free", None, scene_pos, event.modifiers())
            return

        if self.current_tool == "tool_erase":
            if clicked_atom_id is not None:
                self._delete_selection({clicked_atom_id}, set())
            elif clicked_bond_id is not None:
                self._delete_selection(set(), {clicked_bond_id})
            return

        super().mousePressEvent(event)

    def mouseMoveEvent(self, event) -> None:
        scene_pos = self.mapToScene(event.pos())

        if self._ring_dragging:
            self._update_ring_preview(scene_pos, event.modifiers())
            return

        if self._dragging_selection and self._drag_start_pos is not None:
            delta = scene_pos - self._drag_start_pos
            if delta.manhattanLength() > 0:
                self._drag_has_moved = True
                for atom_id, (x, y) in self._drag_start_positions.items():
                    nx = x + delta.x()
                    ny = y + delta.y()
                    self.model.update_atom_position(atom_id, nx, ny)
                    self.update_atom_item(atom_id, nx, ny)
                self.update_bond_items_for_atoms(set(self._drag_start_positions.keys()))
            return

        hover_atom_id = self._get_atom_at(scene_pos.x(), scene_pos.y())
        if self.hover_atom_id is not None and self.hover_atom_id in self.atom_items:
            if self.hover_atom_id not in self.state.selected_atoms and self.hover_atom_id != self.bond_anchor_id:
                self.atom_items[self.hover_atom_id].set_hover(False)
            self.hover_atom_id = None

        if hover_atom_id is not None and hover_atom_id in self.atom_items:
            if hover_atom_id not in self.state.selected_atoms and hover_atom_id != self.bond_anchor_id:
                self.atom_items[hover_atom_id].set_hover(True)
            self.hover_atom_id = hover_atom_id

        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event) -> None:
        if self._ring_dragging:
            self._finalize_ring_drag(event.modifiers())
            return

        if self._dragging_selection:
            if self._drag_has_moved:
                after = {
                    atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
                    for atom_id in self._drag_start_positions
                }
                cmd = MoveAtomsCommand(
                    self.model,
                    self,
                    self._drag_start_positions,
                    after,
                    skip_first_redo=True,
                )
                self.undo_stack.push(cmd)
            self._dragging_selection = False
            self._drag_start_pos = None
            self._drag_start_positions = {}
            self._drag_has_moved = False
            return

        super().mouseReleaseEvent(event)

    def wheelEvent(self, event: QWheelEvent) -> None:
        if event.angleDelta().y() > 0:
            self.zoom_in()
        else:
            self.zoom_out()

    def keyPressEvent(self, event) -> None:
        if event.key() == Qt.Key.Key_Delete:
            self.delete_selection()
            return
        super().keyPressEvent(event)

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
        self.state.selected_atoms.clear()
        self.state.selected_bonds.clear()
        self.bond_anchor_id = None
        self.hover_atom_id = None
        self._create_paper()

    def delete_selection(self) -> None:
        if not self.state.selected_atoms and not self.state.selected_bonds:
            return
        self._delete_selection(set(self.state.selected_atoms), set(self.state.selected_bonds))

    def _delete_selection(self, atom_ids: set[int], bond_ids: set[int]) -> None:
        if not atom_ids and not bond_ids:
            return
        cmd = DeleteSelectionCommand(self.model, self, atom_ids, bond_ids)
        self.undo_stack.push(cmd)
        self.scene.clearSelection()

    def copy_to_clipboard(self) -> None:
        if not self.model.atoms:
            return
        mime = QMimeData()
        try:
            molfile = molgraph_to_molfile(self.model)
            smiles = molgraph_to_smiles(self.model)
            svg = molgraph_to_svg(self.model)
            mime.setData("chemical/x-mdl-molfile", molfile.encode("utf-8"))
            mime.setData("image/svg+xml", svg.encode("utf-8"))
            mime.setText(smiles)
        except Exception:
            pass

        image = self._render_scene_image()
        if image is not None:
            buffer = QBuffer()
            buffer.open(QBuffer.OpenModeFlag.WriteOnly)
            image.save(buffer, "PNG")
            mime.setData("image/png", buffer.data())

        QApplication.clipboard().setMimeData(mime)

    def paste_from_clipboard(self) -> None:
        clipboard = QApplication.clipboard()
        mime = clipboard.mimeData()
        if mime is None:
            return

        if mime.hasFormat("chemical/x-mdl-molfile"):
            molfile = bytes(mime.data("chemical/x-mdl-molfile")).decode("utf-8", errors="ignore")
            try:
                graph = molfile_to_molgraph(molfile)
                self._insert_molgraph(graph)
                return
            except Exception:
                pass

        if mime.hasText():
            smiles = mime.text().strip()
            if smiles:
                try:
                    graph = smiles_to_molgraph(smiles)
                    self._insert_molgraph(graph)
                    return
                except Exception:
                    pass

        if mime.hasFormat("image/png") or mime.hasImage() or mime.hasFormat("image/svg+xml"):
            self._insert_image_from_clipboard(mime)

    def _render_scene_image(self) -> Optional[QImage]:
        items = list(self.atom_items.values()) + list(self.bond_items.values())
        if not items:
            return None
        rect = items[0].sceneBoundingRect()
        for item in items[1:]:
            rect = rect.united(item.sceneBoundingRect())
        rect = rect.adjusted(-10, -10, 10, 10)
        width = max(1, int(rect.width()))
        height = max(1, int(rect.height()))
        image = QImage(width, height, QImage.Format.Format_ARGB32)
        image.fill(Qt.GlobalColor.transparent)
        painter = QPainter(image)
        self.scene.render(painter, QRectF(0, 0, width, height), rect)
        painter.end()
        return image

    def _insert_molgraph(self, graph: MolGraph) -> None:
        if not graph.atoms:
            return
        xs = [atom.x for atom in graph.atoms.values()]
        ys = [atom.y for atom in graph.atoms.values()]
        min_x, max_x = min(xs), max(xs)
        min_y, max_y = min(ys), max(ys)
        center_x = (min_x + max_x) / 2
        center_y = (min_y + max_y) / 2
        target_x = PAPER_WIDTH / 2
        target_y = PAPER_HEIGHT / 2
        dx = target_x - center_x
        dy = target_y - center_y

        self.undo_stack.beginMacro("Paste molecule")
        id_map: Dict[int, int] = {}
        for atom in graph.atoms.values():
            cmd = AddAtomCommand(
                self.model,
                self,
                atom.element,
                atom.x + dx,
                atom.y + dy,
            )
            self.undo_stack.push(cmd)
            if cmd.atom_id is not None:
                id_map[atom.id] = cmd.atom_id
        for bond in graph.bonds.values():
            a1 = id_map.get(bond.a1_id)
            a2 = id_map.get(bond.a2_id)
            if a1 is None or a2 is None:
                continue
            cmd = AddBondCommand(
                self.model,
                self,
                a1,
                a2,
                bond.order,
                bond.style,
                bond.stereo,
            )
            self.undo_stack.push(cmd)
        self.undo_stack.endMacro()

    def _insert_image_from_clipboard(self, mime: QMimeData) -> None:
        pixmap = None
        if mime.hasImage():
            image = mime.imageData()
            if isinstance(image, QImage):
                pixmap = QPixmap.fromImage(image)
            elif isinstance(image, QPixmap):
                pixmap = image
        if pixmap is None and mime.hasFormat("image/png"):
            pixmap = QPixmap()
            pixmap.loadFromData(mime.data("image/png"), "PNG")
        if pixmap is None and mime.hasFormat("image/svg+xml"):
            data = bytes(mime.data("image/svg+xml"))
            pixmap = QPixmap()
            pixmap.loadFromData(data)
        if pixmap is None or pixmap.isNull():
            return
        item = QGraphicsPixmapItem(pixmap)
        item.setZValue(20)
        item.setPos(
            PAPER_WIDTH / 2 - pixmap.width() / 2,
            PAPER_HEIGHT / 2 - pixmap.height() / 2,
        )
        self.scene.addItem(item)

    def add_atom_item(self, atom) -> None:
        if atom.id in self.atom_items:
            return
        item = AtomItem(atom)
        self.scene.addItem(item)
        self.atom_items[atom.id] = item

    def update_atom_item(self, atom_id: int, x: float, y: float) -> None:
        item = self.atom_items.get(atom_id)
        if item is not None:
            item.setPos(x, y)

    def update_atom_item_element(self, atom_id: int, element: str) -> None:
        item = self.atom_items.get(atom_id)
        if item is not None:
            item.set_element(element)

    def remove_atom_item(self, atom_id: int) -> None:
        item = self.atom_items.pop(atom_id, None)
        if item is not None:
            self.scene.removeItem(item)
        if self.bond_anchor_id == atom_id:
            self.bond_anchor_id = None
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

    def update_bond_item(self, bond_id: int) -> None:
        bond = self.model.get_bond(bond_id)
        atom1 = self.model.get_atom(bond.a1_id)
        atom2 = self.model.get_atom(bond.a2_id)
        item = self.bond_items.get(bond_id)
        if item is not None:
            item.set_bond(bond, atom1, atom2)

    def update_bond_items_for_atoms(self, atom_ids: set[int]) -> None:
        for bond in self.model.bonds.values():
            if bond.a1_id in atom_ids or bond.a2_id in atom_ids:
                self.update_bond_item(bond.id)

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

    def _get_bond_at(self, scene_pos: QPointF) -> Optional[int]:
        for item in self.scene.items(scene_pos):
            if isinstance(item, BondItem):
                return item.bond_id
        return None

    def _get_item_at(self, scene_pos: QPointF):
        for item in self.scene.items(scene_pos):
            if isinstance(item, (AtomItem, BondItem)):
                return item
        return None

    def _sync_selection_from_scene(self) -> None:
        selected_atoms = set()
        selected_bonds = set()
        for item in self.scene.selectedItems():
            if isinstance(item, AtomItem):
                selected_atoms.add(item.atom_id)
            elif isinstance(item, BondItem):
                selected_bonds.add(item.bond_id)

        self.state.selected_atoms = selected_atoms
        self.state.selected_bonds = selected_bonds

        for atom_id, item in self.atom_items.items():
            item.set_selected(atom_id in selected_atoms or atom_id == self.bond_anchor_id)

    def _set_bond_anchor(self, atom_id: int, reset_angle: bool = False) -> None:
        if self.bond_anchor_id is not None and self.bond_anchor_id in self.atom_items:
            if self.bond_anchor_id not in self.state.selected_atoms:
                self.atom_items[self.bond_anchor_id].set_selected(False)
        self.bond_anchor_id = atom_id
        if atom_id in self.atom_items:
            self.atom_items[atom_id].set_selected(True)
        if reset_angle:
            self._bond_last_angle = None
            self._bond_zigzag_sign = 1

    def _clear_bond_anchor(self) -> None:
        if self.bond_anchor_id is not None and self.bond_anchor_id in self.atom_items:
            if self.bond_anchor_id not in self.state.selected_atoms:
                self.atom_items[self.bond_anchor_id].set_selected(False)
        self.bond_anchor_id = None
        self._bond_last_angle = None
        self._bond_zigzag_sign = 1

    def _add_bond_between(self, a1_id: int, a2_id: int) -> None:
        order = self.state.active_bond_order
        style = self.state.active_bond_style
        stereo = self.state.active_bond_stereo
        if style != BondStyle.PLAIN:
            order = 1
        cmd = AddBondCommand(self.model, self, a1_id, a2_id, order, style, stereo)
        self.undo_stack.push(cmd)

    def _cycle_bond_order(self, bond_id: int) -> None:
        bond = self.model.get_bond(bond_id)
        if bond.style != BondStyle.PLAIN or bond.stereo != BondStereo.NONE:
            cmd = ChangeBondCommand(
                self.model,
                self,
                bond_id,
                new_order=1,
                new_style=BondStyle.PLAIN,
                new_stereo=BondStereo.NONE,
            )
            self.undo_stack.push(cmd)
            return
        new_order = 1 if bond.order >= 3 else bond.order + 1
        cmd = ChangeBondCommand(self.model, self, bond_id, new_order=new_order)
        self.undo_stack.push(cmd)

    def _apply_bond_style(self, bond_id: int) -> None:
        order = self.state.active_bond_order
        style = self.state.active_bond_style
        stereo = self.state.active_bond_stereo
        if style != BondStyle.PLAIN:
            order = 1
        cmd = ChangeBondCommand(
            self.model,
            self,
            bond_id,
            new_order=order,
            new_style=style,
            new_stereo=stereo,
        )
        self.undo_stack.push(cmd)

    def _create_first_bond(self, scene_pos: QPointF, modifiers: Qt.KeyboardModifiers) -> None:
        self.undo_stack.beginMacro("Add bond")
        atom_cmd = AddAtomCommand(
            self.model,
            self,
            self.state.default_element,
            scene_pos.x(),
            scene_pos.y(),
        )
        self.undo_stack.push(atom_cmd)
        anchor_id = atom_cmd.atom_id
        if anchor_id is not None:
            self._set_bond_anchor(anchor_id, reset_angle=True)
            angle = self._default_bond_angle(anchor_id)
            self._add_bond_with_new_atom(anchor_id, angle)
        self.undo_stack.endMacro()

    def _extend_bond_from_anchor(self, scene_pos: QPointF, modifiers: Qt.KeyboardModifiers) -> None:
        anchor_id = self.bond_anchor_id
        if anchor_id is None:
            return
        anchor = self.model.get_atom(anchor_id)
        dx = scene_pos.x() - anchor.x
        dy = scene_pos.y() - anchor.y
        dist = math.hypot(dx, dy)
        step = self._bond_environment_step(anchor_id)
        angle_from_click = math.atan2(dy, dx) if dist > 1e-6 else 0.0
        existing_angles = self._get_anchor_bond_angles(anchor_id)
        near_anchor = dist < DEFAULT_BOND_LENGTH * 0.6
        aligned_with_existing = False
        if existing_angles:
            nearest = min(self._angle_distance(angle_from_click, a) for a in existing_angles)
            aligned_with_existing = nearest < math.radians(15)

        if near_anchor or aligned_with_existing:
            angle = self._default_bond_angle(anchor_id)
        else:
            angle = angle_from_click
            if not (modifiers & Qt.KeyboardModifier.AltModifier):
                angle = self._snap_angle_to_environment(angle, anchor_id, step)
            self._reset_zigzag()

        self.undo_stack.beginMacro("Add bond")
        self._add_bond_with_new_atom(anchor_id, angle)
        self.undo_stack.endMacro()

    def _add_bond_with_new_atom(self, anchor_id: int, angle: float) -> None:
        anchor = self.model.get_atom(anchor_id)
        new_x = anchor.x + DEFAULT_BOND_LENGTH * math.cos(angle)
        new_y = anchor.y + DEFAULT_BOND_LENGTH * math.sin(angle)
        atom_cmd = AddAtomCommand(
            self.model,
            self,
            self.state.default_element,
            new_x,
            new_y,
        )
        self.undo_stack.push(atom_cmd)
        new_atom_id = atom_cmd.atom_id
        if new_atom_id is None:
            return
        self._add_bond_between(anchor_id, new_atom_id)
        self._record_bond_angle(angle)
        self._set_bond_anchor(new_atom_id, reset_angle=False)

    def _record_bond_angle_between(self, a1_id: int, a2_id: int) -> None:
        a1 = self.model.get_atom(a1_id)
        a2 = self.model.get_atom(a2_id)
        angle = math.atan2(a2.y - a1.y, a2.x - a1.x)
        self._record_bond_angle(angle)

    def _record_bond_angle(self, angle: float) -> None:
        self._bond_last_angle = self._normalize_angle(angle)

    def _default_bond_angle(self, anchor_id: int) -> float:
        step = self._bond_environment_step(anchor_id)
        if self._bond_last_angle is not None:
            angle = self._bond_last_angle + step * self._bond_zigzag_sign
            if step < math.pi:
                self._bond_zigzag_sign *= -1
            return self._normalize_angle(angle)

        existing = self._get_anchor_bond_angles(anchor_id)
        if not existing:
            return 0.0
        candidates = []
        for base in existing:
            candidates.append(base + step)
            candidates.append(base - step)
        return self._best_separated_angle(candidates, existing)

    def _reset_zigzag(self) -> None:
        self._bond_zigzag_sign = 1

    def _snap_angle_to_environment(self, angle: float, anchor_id: int, step: float) -> float:
        if step <= 0:
            return angle
        existing = self._get_anchor_bond_angles(anchor_id)
        if not existing:
            return self._snap_angle(angle, step)
        candidates = []
        for base in existing:
            snapped = base + round((angle - base) / step) * step
            candidates.append(snapped)
        best = min(candidates, key=lambda a: self._angle_distance(a, angle))
        return self._normalize_angle(best)

    def _snap_angle(self, angle: float, step: float) -> float:
        if step <= 0:
            return angle
        return self._normalize_angle(round(angle / step) * step)

    def _bond_environment_step(self, anchor_id: int) -> float:
        max_order = max(1, self.state.active_bond_order)
        for bond in self.model.bonds.values():
            if bond.a1_id == anchor_id or bond.a2_id == anchor_id:
                order = 2 if bond.is_aromatic else bond.order
                max_order = max(max_order, order)
        if max_order >= 3:
            return math.pi
        if max_order == 2:
            return 2 * math.pi / 3
        return math.radians(60)

    def _get_anchor_bond_angles(self, anchor_id: int) -> List[float]:
        angles = []
        anchor = self.model.get_atom(anchor_id)
        for bond in self.model.bonds.values():
            if bond.a1_id == anchor_id:
                other = self.model.get_atom(bond.a2_id)
            elif bond.a2_id == anchor_id:
                other = self.model.get_atom(bond.a1_id)
            else:
                continue
            angle = math.atan2(other.y - anchor.y, other.x - anchor.x)
            angles.append(self._normalize_angle(angle))
        return angles

    def _best_separated_angle(self, candidates: List[float], existing: List[float]) -> float:
        best = None
        best_sep = -1.0
        for candidate in candidates:
            sep = min(self._angle_distance(candidate, a) for a in existing)
            if sep > best_sep:
                best_sep = sep
                best = candidate
        return self._normalize_angle(best) if best is not None else 0.0

    def _normalize_angle(self, angle: float) -> float:
        return (angle + math.pi * 2) % (math.pi * 2)

    def _angle_distance(self, a: float, b: float) -> float:
        diff = (a - b + math.pi) % (2 * math.pi) - math.pi
        return abs(diff)

    def _begin_drag(self, scene_pos: QPointF) -> None:
        if not self.state.selected_atoms:
            return
        self._dragging_selection = True
        self._drag_start_pos = scene_pos
        self._drag_start_positions = {
            atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
            for atom_id in self.state.selected_atoms
        }
        self._drag_has_moved = False

    def _start_ring_drag(
        self,
        anchor_type: str,
        anchor_id: Optional[int],
        scene_pos: QPointF,
        modifiers: Qt.KeyboardModifiers,
    ) -> None:
        self._ring_dragging = True
        self._ring_anchor = {"type": anchor_type, "id": anchor_id, "pos": scene_pos}
        self._ring_last_vertices = None
        if self._ring_preview_item is None:
            self._ring_preview_item = QGraphicsPathItem()
            pen = QPen(QColor("#6699CC"), 1.5, Qt.PenStyle.DashLine)
            self._ring_preview_item.setPen(pen)
            self._ring_preview_item.setBrush(QBrush(Qt.BrushStyle.NoBrush))
            self._ring_preview_item.setZValue(10)
            self.scene.addItem(self._ring_preview_item)
        self._update_ring_preview(scene_pos, modifiers)

    def _update_ring_preview(self, scene_pos: QPointF, modifiers: Qt.KeyboardModifiers) -> None:
        vertices = self._compute_ring_vertices(scene_pos, modifiers)
        self._ring_last_vertices = vertices
        if self._ring_preview_item is None or not vertices:
            return
        path = QPainterPath()
        path.moveTo(vertices[0])
        for v in vertices[1:]:
            path.lineTo(v)
        path.closeSubpath()
        self._ring_preview_item.setPath(path)

    def _finalize_ring_drag(self, modifiers: Qt.KeyboardModifiers) -> None:
        vertices = self._ring_last_vertices or []
        if self._ring_preview_item is not None:
            self.scene.removeItem(self._ring_preview_item)
            self._ring_preview_item = None
        self._ring_dragging = False
        self._ring_last_vertices = None

        if not vertices or self._ring_anchor is None:
            self._ring_anchor = None
            return

        ring_size = self.state.active_ring_size
        aromatic = self.state.active_ring_aromatic
        if modifiers & Qt.KeyboardModifier.ShiftModifier:
            aromatic = not aromatic

        vertex_defs: List[Tuple[Optional[int], float, float]] = []
        if self._ring_anchor["type"] == "bond":
            bond = self.model.get_bond(self._ring_anchor["id"])
            vertex_defs.append((bond.a1_id, vertices[0].x(), vertices[0].y()))
            vertex_defs.append((bond.a2_id, vertices[1].x(), vertices[1].y()))
            for v in vertices[2:]:
                vertex_defs.append((None, v.x(), v.y()))
        elif self._ring_anchor["type"] == "free":
            for v in vertices:
                vertex_defs.append((None, v.x(), v.y()))
        else:
            vertex_defs.append((self._ring_anchor["id"], vertices[0].x(), vertices[0].y()))
            for v in vertices[1:]:
                vertex_defs.append((None, v.x(), v.y()))

        edge_defs: List[Tuple[int, int, int, BondStyle, BondStereo]] = []
        for i in range(ring_size):
            j = (i + 1) % ring_size
            order = 1
            if aromatic and i % 2 == 0:
                order = 2
            edge_defs.append((i, j, order, BondStyle.PLAIN, BondStereo.NONE))

        cmd = AddRingCommand(
            self.model,
            self,
            vertex_defs,
            edge_defs,
            element=self.state.default_element,
        )
        self.undo_stack.push(cmd)
        self._ring_anchor = None

    def _compute_ring_vertices(
        self, scene_pos: QPointF, modifiers: Qt.KeyboardModifiers
    ) -> List[QPointF]:
        if self._ring_anchor is None:
            return []

        ring_size = max(3, int(self.state.active_ring_size))
        flip = bool(modifiers & Qt.KeyboardModifier.AltModifier)

        if self._ring_anchor["type"] == "free":
            center = self._ring_anchor["pos"]
            dx = scene_pos.x() - center.x()
            dy = scene_pos.y() - center.y()
            angle0 = math.atan2(dy, dx) if (dx or dy) else 0.0
            radius = DEFAULT_BOND_LENGTH / (2 * math.sin(math.pi / ring_size))
            vertices = []
            step = 2 * math.pi / ring_size
            for i in range(ring_size):
                angle = angle0 + step * i
                vertices.append(
                    QPointF(
                        center.x() + math.cos(angle) * radius,
                        center.y() + math.sin(angle) * radius,
                    )
                )
            return vertices

        if self._ring_anchor["type"] == "bond":
            bond = self.model.get_bond(self._ring_anchor["id"])
            a1 = self.model.get_atom(bond.a1_id)
            a2 = self.model.get_atom(bond.a2_id)
            p1 = QPointF(a1.x, a1.y)
            p2 = QPointF(a2.x, a2.y)
            return self._regular_polygon_from_edge(p1, p2, ring_size, scene_pos, flip)

        anchor_atom = self.model.get_atom(self._ring_anchor["id"])
        p1 = QPointF(anchor_atom.x, anchor_atom.y)
        direction = scene_pos - p1
        length = math.hypot(direction.x(), direction.y())
        if length < 1.0:
            direction = QPointF(1.0, 0.0)
            length = 1.0
        unit = QPointF(direction.x() / length, direction.y() / length)
        bond_length = DEFAULT_BOND_LENGTH
        p2 = QPointF(p1.x() + unit.x() * bond_length, p1.y() + unit.y() * bond_length)
        return self._regular_polygon_from_edge(p1, p2, ring_size, scene_pos, flip)

    def _regular_polygon_from_edge(
        self,
        p1: QPointF,
        p2: QPointF,
        ring_size: int,
        scene_pos: QPointF,
        flip: bool,
    ) -> List[QPointF]:
        dx = p2.x() - p1.x()
        dy = p2.y() - p1.y()
        length = math.hypot(dx, dy)
        if length < 1e-6:
            return []

        radius = length / (2 * math.sin(math.pi / ring_size))
        mid = QPointF((p1.x() + p2.x()) / 2, (p1.y() + p2.y()) / 2)
        h = math.sqrt(max(radius * radius - (length / 2) ** 2, 0.0))
        nx = -dy / length
        ny = dx / length

        side = dx * (scene_pos.y() - mid.y()) - dy * (scene_pos.x() - mid.x())
        side_sign = 1 if side >= 0 else -1
        if flip:
            side_sign *= -1

        center = QPointF(mid.x() + nx * h * side_sign, mid.y() + ny * h * side_sign)
        angle0 = math.atan2(p1.y() - center.y(), p1.x() - center.x())
        step = 2 * math.pi / ring_size

        cand1 = QPointF(
            center.x() + math.cos(angle0 + step) * radius,
            center.y() + math.sin(angle0 + step) * radius,
        )
        cand2 = QPointF(
            center.x() + math.cos(angle0 - step) * radius,
            center.y() + math.sin(angle0 - step) * radius,
        )
        if (cand1 - p2).manhattanLength() <= (cand2 - p2).manhattanLength():
            direction = 1
        else:
            direction = -1

        vertices = []
        for i in range(ring_size):
            angle = angle0 + direction * step * i
            vertices.append(
                QPointF(
                    center.x() + math.cos(angle) * radius,
                    center.y() + math.sin(angle) * radius,
                )
            )
        return vertices

    def _is_on_paper(self, x: float, y: float) -> bool:
        return (
            PAPER_MARGIN <= x <= PAPER_WIDTH - PAPER_MARGIN
            and PAPER_MARGIN <= y <= PAPER_HEIGHT - PAPER_MARGIN
        )
