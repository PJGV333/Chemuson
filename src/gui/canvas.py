"""
Chemuson Canvas
Page-based canvas using QGraphicsView/QGraphicsScene for document-style editing.
"""
from __future__ import annotations

import math
from typing import Dict, Iterable, List, Optional, Tuple

from PyQt6.QtWidgets import (
    QApplication,
    QGraphicsView,
    QGraphicsScene,
    QGraphicsRectItem,
    QGraphicsDropShadowEffect,
    QGraphicsPixmapItem,
)
from PyQt6.QtGui import (
    QPainter,
    QPen,
    QColor,
    QBrush,
    QWheelEvent,
    QUndoStack,
    QImage,
    QPixmap,
)
from PyQt6.QtCore import Qt, QPointF, QRectF, QBuffer, QMimeData, pyqtSignal

from core.model import BondStyle, BondStereo, ChemState, MolGraph
from gui.items import (
    AtomItem,
    BondItem,
    AromaticCircleItem,
    HoverAtomIndicatorItem,
    HoverBondIndicatorItem,
    OptimizeZoneItem,
    PreviewBondItem,
    PreviewRingItem,
    PreviewChainItem,
    PreviewChainLabelItem,
)
from gui.style import CHEMDOODLE_LIKE, DrawingStyle
from gui.geom import (
    angle_deg,
    endpoint_from_angle_len,
    choose_optimal_direction,
    geometry_for_bond,
    candidate_directions_deg,
    filter_occupied_angles_deg,
    pick_closest_direction_deg,
    angle_distance_deg,
    segments_intersect,
    segment_min_distance,
    segments_nearly_equal,
    closest_atom,
    closest_bond,
    bond_side,
)
from gui.commands import (
    AddAtomCommand,
    AddBondCommand,
    AddRingCommand,
    AddChainCommand,
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
    kekulize_display_orders,
)


# Paper dimensions (A4, approximately)
PAPER_WIDTH = 800
PAPER_HEIGHT = 1000
PAPER_MARGIN = 40  # Margins where atoms shouldn't be placed
ATOM_HIT_RADIUS = 26
DEFAULT_BOND_LENGTH = 40.0
HOVER_ATOM_RADIUS = 20.0
HOVER_BOND_DISTANCE = 14.0
HOVER_BOND_DISTANCE = 10.0
OPTIMIZE_ZONE_SCALE = 1.2
CHAIN_MAX_BONDS = 12
ANGLE_OCCUPIED_TOLERANCE_DEG = 20.0
MIN_ATOM_DIST_SCALE = 0.65
MIN_BOND_DIST_SCALE = 0.2
COLLISION_LENGTH_BOOST = 1.2
BOND_OVERLAP_TOLERANCE_PX = 5.0


class ChemusonCanvas(QGraphicsView):
    """
    Page-based canvas for drawing molecules.
    Uses QGraphicsView/QGraphicsScene with a centered paper sheet.
    """
    
    selection_changed = pyqtSignal(int, int, dict)

    def __init__(self, parent=None) -> None:
        super().__init__(parent)

        self.scene = QGraphicsScene()
        self.setScene(self.scene)

        self.model = MolGraph()
        self.state = ChemState()
        self.undo_stack = QUndoStack(self)
        self.drawing_style: DrawingStyle = CHEMDOODLE_LIKE
        self._ring_centers: dict[int, QPointF] = {}
        self._next_ring_id = 1

        self.atom_items: dict[int, AtomItem] = {}
        self.bond_items: dict[int, BondItem] = {}
        self.aromatic_circles: list[AromaticCircleItem] = []

        self.current_tool: str = self.state.active_tool
        self.bond_anchor_id: Optional[int] = None
        self.hovered_atom_id: Optional[int] = None
        self.hovered_bond_id: Optional[int] = None
        self._last_scene_pos = QPointF(0, 0)
        self._bond_last_angle: Optional[float] = None
        self._bond_zigzag_sign = 1

        self._dragging_selection = False
        self._drag_start_pos: Optional[QPointF] = None
        self._drag_start_positions: Dict[int, Tuple[float, float]] = {}
        self._drag_has_moved = False

        self._drag_mode = "none"
        self._drag_anchor: Optional[dict] = None
        self._ring_last_vertices: Optional[List[QPointF]] = None
        self._chain_last_points: Optional[List[QPointF]] = None
        self._overlays_ready = False

        self._setup_view()
        self._create_paper()
        self._create_overlays()

        self.scene.selectionChanged.connect(self._sync_selection_from_scene)
        self.undo_stack.indexChanged.connect(self._on_undo_stack_changed)
        self.setMouseTracking(True)
        self.viewport().setMouseTracking(True)

    @property
    def graph(self) -> MolGraph:
        """Alias for model for compatibility."""
        return self.model

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

    def _create_overlays(self) -> None:
        self._hover_atom_indicator = HoverAtomIndicatorItem()
        self._hover_bond_indicator = HoverBondIndicatorItem()
        self._optimize_zone = OptimizeZoneItem()
        self._preview_bond_item = PreviewBondItem()
        self._preview_ring_item = PreviewRingItem()
        self._preview_chain_item = PreviewChainItem()
        self._preview_chain_label = PreviewChainLabelItem()

        self.scene.addItem(self._hover_atom_indicator)
        self.scene.addItem(self._hover_bond_indicator)
        self.scene.addItem(self._optimize_zone)
        self.scene.addItem(self._preview_bond_item)
        self.scene.addItem(self._preview_ring_item)
        self.scene.addItem(self._preview_chain_item)
        self.scene.addItem(self._preview_chain_label)
        self._overlays_ready = True

    def set_current_tool(self, tool_id: str) -> None:
        if tool_id.startswith("atom_"):
            self.state.default_element = tool_id.split("_", 1)[1]
            tool_id = "tool_atom"
        elif tool_id.startswith("bond_"):
            tool_id = "tool_bond"

        self.current_tool = tool_id
        self.state.active_tool = tool_id
        self._clear_bond_anchor()
        self._cancel_drag()

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
        self.state.active_bond_aromatic = bond_spec.get("aromatic", False)
        if self.state.active_bond_aromatic:
            self.state.active_bond_order = 1
            self.state.active_bond_style = BondStyle.PLAIN
            self.state.active_bond_stereo = BondStereo.NONE
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
        self._last_scene_pos = scene_pos

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

        clicked_atom_id, clicked_bond_id = self._pick_hover_target(scene_pos)

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
            if clicked_bond_id is not None and clicked_atom_id is None:
                if self.state.active_bond_mode == "increment" and not self.state.active_bond_aromatic:
                    self._cycle_bond_order(clicked_bond_id)
                else:
                    self._apply_bond_style(clicked_bond_id)
                return

            if clicked_atom_id is not None:
                self._begin_place_bond(clicked_atom_id, scene_pos)
                return
            self._begin_place_bond(None, scene_pos)
            return

        if self.current_tool == "tool_ring":
            if clicked_bond_id is not None and clicked_atom_id is None:
                self._begin_place_ring("bond", clicked_bond_id, scene_pos)
                return
            if clicked_atom_id is not None:
                self._begin_place_ring("atom", clicked_atom_id, scene_pos)
                return
            self._begin_place_ring("free", None, scene_pos)
            return

        if self.current_tool == "tool_chain":
            if clicked_atom_id is not None:
                self._begin_place_chain(clicked_atom_id, scene_pos)
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
        self._last_scene_pos = scene_pos

        if self._drag_mode == "place_bond":
            self._update_bond_preview(scene_pos, event.modifiers())
            return
        if self._drag_mode == "place_ring":
            self._update_ring_preview(scene_pos, event.modifiers())
            return
        if self._drag_mode == "place_chain":
            self._update_chain_preview(scene_pos, event.modifiers())
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

        self._update_hover(scene_pos)

        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event) -> None:
        if self._drag_mode == "place_bond":
            self._finalize_bond(event.modifiers())
            return
        if self._drag_mode == "place_ring":
            self._finalize_ring(event.modifiers())
            return
        if self._drag_mode == "place_chain":
            self._finalize_chain(event.modifiers())
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
        if self._handle_hotkeys(event):
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
        self._ring_centers.clear()
        self._next_ring_id = 1
        self.state.selected_atoms.clear()
        self.state.selected_bonds.clear()
        self.bond_anchor_id = None
        self.hovered_atom_id = None
        self.hovered_bond_id = None
        self._overlays_ready = False
        self._create_paper()
        self._create_overlays()

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
                is_aromatic=bond.is_aromatic,
            )
            self.undo_stack.push(cmd)
        self.undo_stack.endMacro()
        if any(bond.is_aromatic for bond in self.model.bonds.values()):
            self._kekulize_aromatic_bonds()

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
        # Determine visibility based on state preferences
        show_c = self.state.show_implicit_carbons or atom.element != "C"
        show_h = self.state.show_implicit_hydrogens or atom.element != "H"
        item = AtomItem(
            atom,
            show_carbon=show_c,
            show_hydrogen=show_h,
            style=self.drawing_style,
        )
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
        if self.hovered_atom_id == atom_id:
            self.hovered_atom_id = None

    def add_bond_item(self, bond) -> None:
        if bond.id in self.bond_items:
            return
        atom1 = self.model.get_atom(bond.a1_id)
        atom2 = self.model.get_atom(bond.a2_id)
        item = BondItem(
            bond,
            atom1,
            atom2,
            render_aromatic_as_circle=self.state.use_aromatic_circles,
            style=self.drawing_style,
        )
        if bond.ring_id is not None:
            item.set_ring_context(self._ring_centers.get(bond.ring_id))
        item.set_offset_sign(self._bond_offset_sign(bond))
        self.scene.addItem(item)
        self.bond_items[bond.id] = item

    def update_bond_item(self, bond_id: int) -> None:
        bond = self.model.get_bond(bond_id)
        atom1 = self.model.get_atom(bond.a1_id)
        atom2 = self.model.get_atom(bond.a2_id)
        item = self.bond_items.get(bond_id)
        if item is not None:
            if bond.ring_id is not None:
                item.set_ring_context(self._ring_centers.get(bond.ring_id))
            else:
                item.set_ring_context(None)
            item.set_offset_sign(self._bond_offset_sign(bond))
            item.set_bond(bond, atom1, atom2)

    def update_bond_items_for_atoms(self, atom_ids: set[int]) -> None:
        for bond in self.model.bonds.values():
            if bond.a1_id in atom_ids or bond.a2_id in atom_ids:
                self.update_bond_item(bond.id)

    def remove_bond_item(self, bond_id: int) -> None:
        item = self.bond_items.pop(bond_id, None)
        if item is not None:
            self.scene.removeItem(item)

    def allocate_ring_id(self) -> int:
        ring_id = self._next_ring_id
        self._next_ring_id += 1
        return ring_id

    def register_ring_center(self, ring_id: int, center: tuple[float, float]) -> None:
        self._ring_centers[ring_id] = QPointF(center[0], center[1])

    def unregister_ring_center(self, ring_id: int) -> None:
        self._ring_centers.pop(ring_id, None)

    def _bond_offset_sign(self, bond) -> int:
        atom1 = self.model.get_atom(bond.a1_id)
        atom2 = self.model.get_atom(bond.a2_id)
        p1 = QPointF(atom1.x, atom1.y)
        p2 = QPointF(atom2.x, atom2.y)
        mid = QPointF((p1.x() + p2.x()) / 2, (p1.y() + p2.y()) / 2)
        dx = p2.x() - p1.x()
        dy = p2.y() - p1.y()
        length = math.hypot(dx, dy) or 1.0
        nx = -dy / length
        ny = dx / length

        def neighbor_score(atom_id: int) -> float:
            score = 0.0
            for other_bond in self.model.bonds.values():
                if other_bond.id == bond.id:
                    continue
                if other_bond.a1_id == atom_id:
                    other = self.model.get_atom(other_bond.a2_id)
                elif other_bond.a2_id == atom_id:
                    other = self.model.get_atom(other_bond.a1_id)
                else:
                    continue
                vx = other.x - mid.x()
                vy = other.y - mid.y()
                score += nx * vx + ny * vy
            return score

        score = neighbor_score(atom1.id) + neighbor_score(atom2.id)
        return -1 if score > 0 else 1

    # -------------------------------------------------------------------------
    # Visibility and Aromatic Circle Management
    # -------------------------------------------------------------------------
    def refresh_atom_visibility(self) -> None:
        """Update visibility of all atoms based on current state settings."""
        for atom_id, item in self.atom_items.items():
            atom = self.model.get_atom(atom_id)
            show_c = self.state.show_implicit_carbons or atom.element != "C"
            show_h = self.state.show_implicit_hydrogens or atom.element != "H"
            item.set_visibility_flags(show_c, show_h)

    def refresh_aromatic_circles(self) -> None:
        """Add/remove aromatic circle items based on current state."""
        # First, update all bond items to reflect the preference
        use_circles = self.state.use_aromatic_circles
        for bond_id, item in self.bond_items.items():
            if item.is_aromatic:
                item.set_render_aromatic_as_circle(use_circles)
                # Force redraw
                self.update_bond_item(bond_id)

        # Remove existing circles
        for circle in self.aromatic_circles:
            self.scene.removeItem(circle)
        self.aromatic_circles.clear()
        
        if not use_circles:
            return
        
        # Find aromatic rings and add circles
        rings = self._detect_aromatic_rings()
        for center_x, center_y, radius in rings:
            circle = AromaticCircleItem(center_x, center_y, radius)
            self.scene.addItem(circle)
            self.aromatic_circles.append(circle)

    def _detect_aromatic_rings(self) -> list:
        """
        Detect complete aromatic rings and return their centers.
        Returns list of (center_x, center_y, radius) tuples.
        """
        # Collect aromatic bonds
        aromatic_bonds = [b for b in self.model.bonds.values() if b.is_aromatic]
        if not aromatic_bonds:
            return []
        
        # Build adjacency for aromatic atoms only
        aromatic_atoms = set()
        adjacency: dict[int, list[int]] = {}
        for bond in aromatic_bonds:
            aromatic_atoms.add(bond.a1_id)
            aromatic_atoms.add(bond.a2_id)
            adjacency.setdefault(bond.a1_id, []).append(bond.a2_id)
            adjacency.setdefault(bond.a2_id, []).append(bond.a1_id)
        
        # Find rings using simple DFS for small rings (5-7 members)
        visited_rings = set()
        rings = []
        
        for start in aromatic_atoms:
            # Try to find rings starting from this atom
            ring = self._find_ring_from(start, adjacency, max_size=7)
            if ring:
                ring_key = tuple(sorted(ring))
                if ring_key not in visited_rings:
                    visited_rings.add(ring_key)
                    rings.append(ring)
        
        # Calculate centers and radii for each ring
        result = []
        for ring in rings:
            atoms = [self.model.get_atom(aid) for aid in ring]
            cx = sum(a.x for a in atoms) / len(atoms)
            cy = sum(a.y for a in atoms) / len(atoms)
            # Radius is about 60% of distance from center to atoms
            avg_dist = sum(((a.x - cx)**2 + (a.y - cy)**2)**0.5 for a in atoms) / len(atoms)
            radius = avg_dist * 0.55
            result.append((cx, cy, radius))
        
        return result

    def _find_ring_from(self, start: int, adjacency: dict, max_size: int) -> Optional[list]:
        """Find a ring starting from the given atom using DFS."""
        # Simple BFS to find shortest cycle
        from collections import deque
        
        queue = deque([(start, [start])])
        visited = {start}
        
        while queue:
            current, path = queue.popleft()
            if len(path) > max_size:
                continue
            
            for neighbor in adjacency.get(current, []):
                if neighbor == start and len(path) >= 5:
                    return path
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append((neighbor, path + [neighbor]))
        
        return None

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
                if item.atom_id in self.model.atoms:
                    selected_atoms.add(item.atom_id)
            elif isinstance(item, BondItem):
                if item.bond_id in self.model.bonds:
                    selected_bonds.add(item.bond_id)

        self.state.selected_atoms = selected_atoms
        self.state.selected_bonds = selected_bonds

        for atom_id, item in self.atom_items.items():
            item.set_selected(atom_id in selected_atoms or atom_id == self.bond_anchor_id)
            
        # Emit selection signal
        details = {}
        if len(selected_atoms) == 1 and not selected_bonds:
            atom = self.model.get_atom(next(iter(selected_atoms)))
            details = {"type": "atom", "id": atom.id, "element": atom.element, "charge": atom.charge, "x": atom.x, "y": atom.y}
        elif len(selected_bonds) == 1 and not selected_atoms:
            bond_id = next(iter(selected_bonds))
            if bond_id in self.model.bonds:
                bond = self.model.get_bond(bond_id)
                details = {"type": "bond", "id": bond.id, "order": bond.order, "style": bond.style, "aromatic": bond.is_aromatic}
            
        self.selection_changed.emit(len(selected_atoms), len(selected_bonds), details)

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

    def _create_or_update_bond(
        self,
        a1_id: int,
        a2_id: int,
        order: int,
        style: BondStyle,
        stereo: BondStereo,
        is_aromatic: bool,
    ) -> None:
        existing = self.model.find_bond_between(a1_id, a2_id)
        if existing is not None:
            if (
                self.state.active_bond_mode == "increment"
                and not is_aromatic
                and style == BondStyle.PLAIN
                and stereo == BondStereo.NONE
            ):
                new_order = 1 if existing.order >= 3 else existing.order + 1
                new_style = existing.style
                new_stereo = existing.stereo
                new_is_aromatic = False
            else:
                new_order = order
                new_style = style
                new_stereo = stereo
                new_is_aromatic = is_aromatic
            cmd = ChangeBondCommand(
                self.model,
                self,
                existing.id,
                new_order=new_order,
                new_style=new_style,
                new_stereo=new_stereo,
                new_is_aromatic=new_is_aromatic,
            )
            self.undo_stack.push(cmd)
            if new_is_aromatic:
                self._kekulize_aromatic_bonds(seed_atoms={a1_id, a2_id})
            return

        cmd = AddBondCommand(
            self.model,
            self,
            a1_id,
            a2_id,
            order,
            style,
            stereo,
            is_aromatic=is_aromatic,
        )
        self.undo_stack.push(cmd)
        if is_aromatic:
            self._kekulize_aromatic_bonds(seed_atoms={a1_id, a2_id})

    def _add_bond_between(self, a1_id: int, a2_id: int) -> None:
        order = self.state.active_bond_order
        style = self.state.active_bond_style
        stereo = self.state.active_bond_stereo
        is_aromatic = self.state.active_bond_aromatic
        if style != BondStyle.PLAIN or is_aromatic:
            order = 1
        self._create_or_update_bond(a1_id, a2_id, order, style, stereo, is_aromatic)

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
                new_is_aromatic=False,
            )
            self.undo_stack.push(cmd)
            return
        new_order = 1 if bond.order >= 3 else bond.order + 1
        cmd = ChangeBondCommand(
            self.model,
            self,
            bond_id,
            new_order=new_order,
            new_is_aromatic=False,
        )
        self.undo_stack.push(cmd)

    def _apply_bond_style(self, bond_id: int) -> None:
        order = self.state.active_bond_order
        style = self.state.active_bond_style
        stereo = self.state.active_bond_stereo
        is_aromatic = self.state.active_bond_aromatic
        if style != BondStyle.PLAIN:
            order = 1
        if is_aromatic:
            order = 1
            style = BondStyle.PLAIN
            stereo = BondStereo.NONE
        cmd = ChangeBondCommand(
            self.model,
            self,
            bond_id,
            new_order=order,
            new_style=style,
            new_stereo=stereo,
            new_is_aromatic=is_aromatic,
        )
        self.undo_stack.push(cmd)
        if is_aromatic:
            bond = self.model.get_bond(bond_id)
            self._kekulize_aromatic_bonds(seed_atoms={bond.a1_id, bond.a2_id})

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
        if not (modifiers & Qt.KeyboardModifier.AltModifier) and len(existing_angles) <= 1:
            angle = self._default_bond_angle(anchor_id)
        else:
            near_anchor = dist < DEFAULT_BOND_LENGTH * 0.6
            aligned_with_existing = False
            if existing_angles:
                nearest = min(self._angle_distance(angle_from_click, a) for a in existing_angles)
                aligned_with_existing = nearest < math.radians(30)

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

    # -------------------------------------------------------------------------
    # Hover + Drag State
    # -------------------------------------------------------------------------
    def _pick_hover_target(self, scene_pos: QPointF) -> tuple[Optional[int], Optional[int]]:
        atoms = [(atom.id, atom.x, atom.y) for atom in self.model.atoms.values()]
        atom_id = closest_atom(scene_pos, atoms, HOVER_ATOM_RADIUS)
        if atom_id is not None:
            return atom_id, None

        bonds = []
        for bond in self.model.bonds.values():
            a1 = self.model.get_atom(bond.a1_id)
            a2 = self.model.get_atom(bond.a2_id)
            bonds.append((bond.id, QPointF(a1.x, a1.y), QPointF(a2.x, a2.y)))
        bond_id = closest_bond(scene_pos, bonds, HOVER_BOND_DISTANCE)
        return None, bond_id

    def _update_hover(self, scene_pos: QPointF) -> None:
        atom_id, bond_id = self._pick_hover_target(scene_pos)

        if self.hovered_atom_id is not None and self.hovered_atom_id in self.atom_items:
            if self.hovered_atom_id not in self.state.selected_atoms:
                self.atom_items[self.hovered_atom_id].set_hover(False)

        self.hovered_atom_id = atom_id
        self.hovered_bond_id = bond_id if atom_id is None else None

        if self.hovered_atom_id is not None and self.hovered_atom_id in self.atom_items:
            if self.hovered_atom_id not in self.state.selected_atoms:
                self.atom_items[self.hovered_atom_id].set_hover(True)
            atom = self.model.get_atom(self.hovered_atom_id)
            self._hover_atom_indicator.update_position(atom.x, atom.y)
            self._hover_bond_indicator.hide_indicator()
            return

        self._hover_atom_indicator.hide_indicator()

        if self.hovered_bond_id is not None:
            bond = self.model.get_bond(self.hovered_bond_id)
            a1 = self.model.get_atom(bond.a1_id)
            a2 = self.model.get_atom(bond.a2_id)
            self._hover_bond_indicator.update_for_bond(
                QPointF(a1.x, a1.y),
                QPointF(a2.x, a2.y),
            )
        else:
            self._hover_bond_indicator.hide_indicator()

    def _cancel_drag(self) -> None:
        self._drag_mode = "none"
        self._drag_anchor = None
        self._ring_last_vertices = None
        self._chain_last_points = None
        if self._overlays_ready:
            self._optimize_zone.hide_zone()
            self._preview_bond_item.hide_preview()
            self._preview_ring_item.hide_preview()
            self._preview_chain_item.hide_preview()
            self._preview_chain_label.hide_label()

    def _on_undo_stack_changed(self, _index: int) -> None:
        self._cancel_drag()
        self._sync_scene_with_model()
        self._kekulize_aromatic_bonds()
        self._update_hover(self._last_scene_pos)

    def _sync_scene_with_model(self) -> None:
        for item in list(self.scene.items()):
            if isinstance(item, (AtomItem, BondItem)):
                self.scene.removeItem(item)

        for atom_id in list(self.atom_items.keys()):
            self.remove_atom_item(atom_id)
        for bond_id in list(self.bond_items.keys()):
            self.remove_bond_item(bond_id)

        self._ring_centers.clear()
        ring_atoms: dict[int, set[int]] = {}
        for bond in self.model.bonds.values():
            if bond.ring_id is None:
                continue
            ring_atoms.setdefault(bond.ring_id, set()).update({bond.a1_id, bond.a2_id})
        for ring_id, atom_ids in ring_atoms.items():
            xs = [self.model.get_atom(aid).x for aid in atom_ids]
            ys = [self.model.get_atom(aid).y for aid in atom_ids]
            center = (sum(xs) / len(xs), sum(ys) / len(ys))
            self.register_ring_center(ring_id, center)

        for atom in self.model.atoms.values():
            self.add_atom_item(atom)
        for bond in self.model.bonds.values():
            self.add_bond_item(bond)

        self.refresh_atom_visibility()
        self.refresh_aromatic_circles()

    def _get_anchor_bond_angles_deg(self, anchor_id: int) -> list[float]:
        angles = []
        anchor = self.model.get_atom(anchor_id)
        p0 = QPointF(anchor.x, anchor.y)
        for bond in self.model.bonds.values():
            if bond.a1_id == anchor_id:
                other = self.model.get_atom(bond.a2_id)
            elif bond.a2_id == anchor_id:
                other = self.model.get_atom(bond.a1_id)
            else:
                continue
            p1 = QPointF(other.x, other.y)
            angles.append(angle_deg(p0, p1))
        return angles

    def _bond_geometry(self, anchor_id: int, bond_order: int, is_aromatic: bool) -> str:
        has_triple = bond_order >= 3
        has_double = bond_order == 2 or is_aromatic
        for bond in self.model.bonds.values():
            if bond.a1_id != anchor_id and bond.a2_id != anchor_id:
                continue
            if bond.order >= 3:
                has_triple = True
            if bond.order == 2 or bond.is_aromatic:
                has_double = True
        if has_triple:
            return "sp"
        if has_double:
            return "sp2"
        return "sp3"

    def _incoming_angle_deg(self, anchor_id: int) -> Optional[float]:
        angles = self._get_anchor_bond_angles_deg(anchor_id)
        if len(angles) == 1:
            return angles[0]
        return None

    def _direction_collision_metrics(
        self, anchor_id: Optional[int], p0: QPointF, p1: QPointF
    ) -> tuple[int, float, float]:
        intersections = 0
        min_atom_dist = float("inf")
        min_bond_dist = float("inf")
        atom_threshold = self.state.bond_length * MIN_ATOM_DIST_SCALE
        bond_threshold = self.state.bond_length * MIN_BOND_DIST_SCALE
        excluded_atoms: set[int] = set()
        if anchor_id is not None:
            excluded_atoms.add(anchor_id)
        target_id = closest_atom(
            p1,
            [(atom.id, atom.x, atom.y) for atom in self.model.atoms.values()],
            ATOM_HIT_RADIUS,
        )
        if target_id is not None:
            excluded_atoms.add(target_id)

        for atom in self.model.atoms.values():
            if atom.id in excluded_atoms:
                continue
            dist = math.hypot(p1.x() - atom.x, p1.y() - atom.y)
            min_atom_dist = min(min_atom_dist, dist)
            if dist < atom_threshold:
                intersections += 1

        for bond in self.model.bonds.values():
            if bond.a1_id in excluded_atoms or bond.a2_id in excluded_atoms:
                continue
            a1 = self.model.get_atom(bond.a1_id)
            a2 = self.model.get_atom(bond.a2_id)
            pA = QPointF(a1.x, a1.y)
            pB = QPointF(a2.x, a2.y)
            if segments_intersect(p0, p1, pA, pB):
                intersections += 1
            bond_dist = segment_min_distance(p0, p1, pA, pB)
            min_bond_dist = min(min_bond_dist, bond_dist)
            if bond_dist < bond_threshold:
                intersections += 1

        return intersections, min_atom_dist, min_bond_dist

    def _direction_score(
        self,
        anchor_id: Optional[int],
        p0: QPointF,
        p1: QPointF,
        angle_deg_value: float,
        mouse_angle_deg: float,
        preferred_deg: list[float],
    ) -> float:
        intersections, min_atom_dist, min_bond_dist = self._direction_collision_metrics(
            anchor_id, p0, p1
        )
        atom_threshold = self.state.bond_length * MIN_ATOM_DIST_SCALE
        bond_threshold = self.state.bond_length * MIN_BOND_DIST_SCALE

        score = angle_distance_deg(angle_deg_value, mouse_angle_deg)
        if preferred_deg:
            if min(angle_distance_deg(angle_deg_value, pref) for pref in preferred_deg) <= 15.0:
                score -= 15.0
        score += intersections * 100.0
        if min_atom_dist < atom_threshold:
            score += (atom_threshold - min_atom_dist) * 5.0
        if min_bond_dist < bond_threshold:
            score += (bond_threshold - min_bond_dist) * 5.0
        return score

    def _direction_is_valid(
        self,
        anchor_id: Optional[int],
        p0: QPointF,
        p1: QPointF,
    ) -> bool:
        intersections, min_atom_dist, min_bond_dist = self._direction_collision_metrics(
            anchor_id, p0, p1
        )
        atom_threshold = self.state.bond_length * MIN_ATOM_DIST_SCALE
        bond_threshold = self.state.bond_length * MIN_BOND_DIST_SCALE
        if intersections > 0:
            return False
        if min_atom_dist < atom_threshold:
            return False
        if min_bond_dist < bond_threshold:
            return False
        return True

    def _pick_bond_direction_deg(
        self,
        anchor: QPointF,
        anchor_id: Optional[int],
        mouse_angle_deg: float,
        bond_order: int,
        is_aromatic: bool,
        length: float,
        apply_collisions: bool = True,
        allow_length_boost: bool = True,
    ) -> tuple[float, float]:
        existing_angles = self._get_anchor_bond_angles_deg(anchor_id) if anchor_id else []
        geometry = (
            self._bond_geometry(anchor_id, bond_order, is_aromatic)
            if anchor_id is not None
            else geometry_for_bond(bond_order, is_aromatic, [])
        )
        candidates = candidate_directions_deg(geometry, existing_angles, mouse_angle_deg)
        candidates = filter_occupied_angles_deg(
            candidates, existing_angles, ANGLE_OCCUPIED_TOLERANCE_DEG
        )
        if not candidates:
            candidates = candidate_directions_deg(geometry, [], mouse_angle_deg)

        preferred: list[float] = []
        if geometry == "sp3" and anchor_id is not None:
            incoming = self._incoming_angle_deg(anchor_id)
            if incoming is not None:
                preferred = [(incoming + 60.0) % 360.0, (incoming - 60.0) % 360.0]

        if not apply_collisions:
            picked = pick_closest_direction_deg(candidates, mouse_angle_deg, preferred)
            return (picked if picked is not None else mouse_angle_deg), length

        length_candidates = (length, length * COLLISION_LENGTH_BOOST) if allow_length_boost else (length,)
        valid_choices: list[tuple[float, float, float]] = []
        for length_candidate in length_candidates:
            for angle_candidate in candidates:
                p1 = endpoint_from_angle_len(anchor, angle_candidate, length_candidate)
                if not self._direction_is_valid(anchor_id, anchor, p1):
                    continue
                score = self._direction_score(
                    anchor_id,
                    anchor,
                    p1,
                    angle_candidate,
                    mouse_angle_deg,
                    preferred,
                )
                valid_choices.append((score, angle_candidate, length_candidate))

        if valid_choices:
            _, best_angle, best_length = min(valid_choices, key=lambda item: item[0])
            return best_angle, best_length

        best_score = None
        best_angle = None
        best_length = length
        for length_candidate in length_candidates:
            for angle_candidate in candidates:
                p1 = endpoint_from_angle_len(anchor, angle_candidate, length_candidate)
                score = self._direction_score(
                    anchor_id,
                    anchor,
                    p1,
                    angle_candidate,
                    mouse_angle_deg,
                    preferred,
                )
                if best_score is None or score < best_score:
                    best_score = score
                    best_angle = angle_candidate
                    best_length = length_candidate

        return (best_angle if best_angle is not None else mouse_angle_deg), best_length

    def _select_preferred_angle(
        self, anchor_id: int, cursor_angle: float, bond_order: int, is_aromatic: bool
    ) -> float:
        anchor = self.model.get_atom(anchor_id)
        p0 = QPointF(anchor.x, anchor.y)
        angle, _ = self._pick_bond_direction_deg(
            p0,
            anchor_id,
            cursor_angle,
            bond_order,
            is_aromatic,
            self.state.bond_length,
            apply_collisions=False,
            allow_length_boost=False,
        )
        return angle

    def _find_overlapping_bond(self, p0: QPointF, p1: QPointF) -> Optional[int]:
        for bond in self.model.bonds.values():
            a1 = self.model.get_atom(bond.a1_id)
            a2 = self.model.get_atom(bond.a2_id)
            b0 = QPointF(a1.x, a1.y)
            b1 = QPointF(a2.x, a2.y)
            if segments_nearly_equal(p0, p1, b0, b1, BOND_OVERLAP_TOLERANCE_PX):
                return bond.id
        return None

    def _snap_ring_vertex(self, pos: QPointF, excluded: set[int]) -> Optional[int]:
        candidates = [
            (atom.id, atom.x, atom.y)
            for atom in self.model.atoms.values()
            if atom.id not in excluded
        ]
        return closest_atom(pos, candidates, ATOM_HIT_RADIUS)

    def _build_ring_vertex_defs(
        self, vertices: List[QPointF], anchor_type: str, anchor_id: Optional[int]
    ) -> List[Tuple[Optional[int], float, float]]:
        vertex_defs: List[Tuple[Optional[int], float, float]] = []
        used_ids: set[int] = set()

        if anchor_type == "bond" and anchor_id is not None:
            bond = self.model.get_bond(anchor_id)
            used_ids.update({bond.a1_id, bond.a2_id})
            vertex_defs.append((bond.a1_id, vertices[0].x(), vertices[0].y()))
            vertex_defs.append((bond.a2_id, vertices[1].x(), vertices[1].y()))
            remaining_vertices = vertices[2:]
        elif anchor_type == "atom" and anchor_id is not None:
            used_ids.add(anchor_id)
            vertex_defs.append((anchor_id, vertices[0].x(), vertices[0].y()))
            remaining_vertices = vertices[1:]
        else:
            remaining_vertices = vertices

        for v in remaining_vertices:
            snapped_id = self._snap_ring_vertex(v, used_ids)
            if snapped_id is not None:
                used_ids.add(snapped_id)
                vertex_defs.append((snapped_id, v.x(), v.y()))
            else:
                vertex_defs.append((None, v.x(), v.y()))

        return vertex_defs

    # -------------------------------------------------------------------------
    # Bond Placement
    # -------------------------------------------------------------------------
    def _begin_place_bond(self, anchor_id: Optional[int], scene_pos: QPointF) -> None:
        self._drag_mode = "place_bond"
        if anchor_id is None:
            self._drag_anchor = {"type": "free", "id": None, "pos": scene_pos}
            self._optimize_zone.hide_zone()
        else:
            anchor = self.model.get_atom(anchor_id)
            self._drag_anchor = {"type": "atom", "id": anchor_id, "pos": QPointF(anchor.x, anchor.y)}
            radius = self.state.bond_length * OPTIMIZE_ZONE_SCALE
            self._optimize_zone.set_radius(radius)
            self._optimize_zone.update_center(anchor.x, anchor.y)
        self._update_bond_preview(scene_pos, Qt.KeyboardModifier.NoModifier)

    def _update_bond_preview(self, scene_pos: QPointF, modifiers: Qt.KeyboardModifiers) -> None:
        if self._drag_anchor is None:
            return
        if self._drag_anchor["id"] is None:
            p0 = QPointF(self._drag_anchor["pos"].x(), self._drag_anchor["pos"].y())
        else:
            anchor = self.model.get_atom(self._drag_anchor["id"])
            p0 = QPointF(anchor.x, anchor.y)
        p1 = self._compute_bond_endpoint(p0, scene_pos, modifiers)
        self._preview_bond_item.update_line(p0, p1)

    def _finalize_bond(self, modifiers: Qt.KeyboardModifiers) -> None:
        if self._drag_anchor is None:
            self._cancel_drag()
            return
        anchor_id = self._drag_anchor["id"]
        if anchor_id is None:
            p0 = QPointF(self._drag_anchor["pos"].x(), self._drag_anchor["pos"].y())
        else:
            anchor = self.model.get_atom(anchor_id)
            p0 = QPointF(anchor.x, anchor.y)
        p1 = self._compute_bond_endpoint(p0, self._last_scene_pos, modifiers)

        atom_positions = [(atom.id, atom.x, atom.y) for atom in self.model.atoms.values()]
        target_id = closest_atom(self._last_scene_pos, atom_positions, ATOM_HIT_RADIUS)
        if target_id is None:
            target_id = closest_atom(p1, atom_positions, ATOM_HIT_RADIUS)
        if target_id == anchor_id:
            target_id = None

        if target_id is None:
            overlapping_bond = self._find_overlapping_bond(p0, p1)
            if overlapping_bond is not None:
                if self.state.active_bond_mode == "increment" and not self.state.active_bond_aromatic:
                    self._cycle_bond_order(overlapping_bond)
                else:
                    self._apply_bond_style(overlapping_bond)
                self._cancel_drag()
                return
            bond_hit = self._get_bond_at(self._last_scene_pos)
            if bond_hit is not None:
                if self.state.active_bond_mode == "increment" and not self.state.active_bond_aromatic:
                    self._cycle_bond_order(bond_hit)
                else:
                    self._apply_bond_style(bond_hit)
                self._cancel_drag()
                return

        order = self.state.active_bond_order
        style = self.state.active_bond_style
        stereo = self.state.active_bond_stereo
        is_aromatic = self.state.active_bond_aromatic
        if style != BondStyle.PLAIN or is_aromatic:
            order = 1

        if anchor_id is None:
            self.undo_stack.beginMacro("Add bond")
            anchor_cmd = AddAtomCommand(
                self.model,
                self,
                self.state.default_element,
                p0.x(),
                p0.y(),
            )
            self.undo_stack.push(anchor_cmd)
            anchor_id = anchor_cmd.atom_id
            if target_id is not None:
                self._create_or_update_bond(
                    anchor_id,
                    target_id,
                    order,
                    style,
                    stereo,
                    is_aromatic,
                )
            else:
                cmd = AddBondCommand(
                    self.model,
                    self,
                    anchor_id,
                    None,
                    order,
                    style,
                    stereo,
                    is_aromatic=is_aromatic,
                    new_atom_element=self.state.default_element,
                    new_atom_pos=(p1.x(), p1.y()),
                )
                self.undo_stack.push(cmd)
                if is_aromatic:
                    self._kekulize_aromatic_bonds()
            self.undo_stack.endMacro()
        else:
            if target_id is not None:
                self._create_or_update_bond(
                    anchor_id,
                    target_id,
                    order,
                    style,
                    stereo,
                    is_aromatic,
                )
            else:
                cmd = AddBondCommand(
                    self.model,
                    self,
                    anchor_id,
                    None,
                    order,
                    style,
                    stereo,
                    is_aromatic=is_aromatic,
                    new_atom_element=self.state.default_element,
                    new_atom_pos=(p1.x(), p1.y()),
                )
                self.undo_stack.push(cmd)
                if is_aromatic:
                    self._kekulize_aromatic_bonds()
        self._cancel_drag()

    def _compute_bond_endpoint(
        self,
        anchor: QPointF,
        scene_pos: QPointF,
        modifiers: Qt.KeyboardModifiers,
        bond_order: Optional[int] = None,
        is_aromatic: Optional[bool] = None,
    ) -> QPointF:
        dx = scene_pos.x() - anchor.x()
        dy = scene_pos.y() - anchor.y()
        dist = math.hypot(dx, dy)

        use_alt = bool(modifiers & Qt.KeyboardModifier.AltModifier)
        use_shift = bool(modifiers & Qt.KeyboardModifier.ShiftModifier)
        use_optimize = False
        if self._optimize_zone.isVisible():
            use_optimize = dist <= self._optimize_zone.radius()

        cursor_theta = angle_deg(anchor, scene_pos)
        if not self.state.fixed_angles or use_alt:
            theta = cursor_theta
            length = self.state.bond_length if self.state.fixed_lengths and not use_shift else max(5.0, dist)
            p1 = endpoint_from_angle_len(anchor, theta, length)
            snap_id = closest_atom(
                p1,
                [(atom.id, atom.x, atom.y) for atom in self.model.atoms.values()],
                ATOM_HIT_RADIUS,
            )
            if snap_id is not None and (self._drag_anchor is None or snap_id != self._drag_anchor["id"]):
                atom = self.model.get_atom(snap_id)
                return QPointF(atom.x, atom.y)
            return p1

        order = bond_order if bond_order is not None else self.state.active_bond_order
        aromatic_flag = is_aromatic if is_aromatic is not None else self.state.active_bond_aromatic
        style = (
            self.state.active_bond_style if bond_order is None and is_aromatic is None else BondStyle.PLAIN
        )
        if style != BondStyle.PLAIN or aromatic_flag:
            order = 1

        if use_optimize and self._drag_anchor["id"] is not None:
            cursor_theta = self._select_preferred_angle(
                self._drag_anchor["id"], cursor_theta, order, aromatic_flag
            )

        length = self.state.bond_length if self.state.fixed_lengths and not use_shift else max(5.0, dist)
        theta, final_length = self._pick_bond_direction_deg(
            anchor,
            self._drag_anchor["id"] if self._drag_anchor else None,
            cursor_theta,
            order,
            aromatic_flag,
            length,
            apply_collisions=True,
            allow_length_boost=self.state.fixed_lengths and not use_shift,
        )
        p1 = endpoint_from_angle_len(anchor, theta, final_length)
        snap_id = closest_atom(
            p1,
            [(atom.id, atom.x, atom.y) for atom in self.model.atoms.values()],
            ATOM_HIT_RADIUS,
        )
        if snap_id is not None and (self._drag_anchor is None or snap_id != self._drag_anchor["id"]):
            atom = self.model.get_atom(snap_id)
            return QPointF(atom.x, atom.y)
        return p1

    # -------------------------------------------------------------------------
    # Ring Placement
    # -------------------------------------------------------------------------
    def _begin_place_ring(self, anchor_type: str, anchor_id: Optional[int], scene_pos: QPointF) -> None:
        self._drag_mode = "place_ring"
        self._drag_anchor = {"type": anchor_type, "id": anchor_id, "pos": scene_pos}
        if anchor_type == "atom":
            atom = self.model.get_atom(anchor_id)
            radius = self.state.bond_length * OPTIMIZE_ZONE_SCALE
            self._optimize_zone.set_radius(radius)
            self._optimize_zone.update_center(atom.x, atom.y)
        elif anchor_type == "bond":
            bond = self.model.get_bond(anchor_id)
            a1 = self.model.get_atom(bond.a1_id)
            a2 = self.model.get_atom(bond.a2_id)
            mid = QPointF((a1.x + a2.x) / 2, (a1.y + a2.y) / 2)
            radius = self.state.bond_length * OPTIMIZE_ZONE_SCALE
            self._optimize_zone.set_radius(radius)
            self._optimize_zone.update_center(mid.x(), mid.y())
        else:
            self._optimize_zone.hide_zone()
        self._update_ring_preview(scene_pos, Qt.KeyboardModifier.NoModifier)

    def _update_ring_preview(self, scene_pos: QPointF, modifiers: Qt.KeyboardModifiers) -> None:
        vertices = self._compute_ring_vertices(scene_pos, modifiers)
        self._ring_last_vertices = vertices
        self._preview_ring_item.update_polygon(vertices)

    def _finalize_ring(self, modifiers: Qt.KeyboardModifiers) -> None:
        vertices = self._ring_last_vertices or []
        if not vertices or self._drag_anchor is None:
            self._cancel_drag()
            return

        ring_size = max(3, int(self.state.active_ring_size))
        aromatic = self.state.active_ring_aromatic

        vertex_defs = self._build_ring_vertex_defs(
            vertices, self._drag_anchor["type"], self._drag_anchor["id"]
        )

        edge_defs: List[Tuple[int, int, int, BondStyle, BondStereo, bool]] = []
        if aromatic:
            edge_defs = self._build_aromatic_edges(vertex_defs, ring_size)
        else:
            for i in range(ring_size):
                j = (i + 1) % ring_size
                edge_defs.append((i, j, 1, BondStyle.PLAIN, BondStereo.NONE, aromatic))

        cmd = AddRingCommand(
            self.model,
            self,
            vertex_defs,
            edge_defs,
            element=self.state.default_element,
        )
        self.undo_stack.push(cmd)
        if aromatic:
            self._kekulize_aromatic_bonds()
        self._cancel_drag()

    def _compute_ring_vertices(
        self, scene_pos: QPointF, modifiers: Qt.KeyboardModifiers
    ) -> List[QPointF]:
        if self._drag_anchor is None:
            return []

        ring_size = max(3, int(self.state.active_ring_size))
        anchor_type = self._drag_anchor["type"]

        if anchor_type == "free":
            center = self._drag_anchor["pos"]
            angle0 = angle_deg(center, scene_pos)
            radius = self.state.bond_length / (2 * math.sin(math.pi / ring_size))
            vertices = []
            step = 360.0 / ring_size
            for i in range(ring_size):
                theta = angle0 + step * i
                vertices.append(endpoint_from_angle_len(center, theta, radius))
            return vertices

        if anchor_type == "bond":
            bond = self.model.get_bond(self._drag_anchor["id"])
            a1 = self.model.get_atom(bond.a1_id)
            a2 = self.model.get_atom(bond.a2_id)
            p1 = QPointF(a1.x, a1.y)
            p2 = QPointF(a2.x, a2.y)
            direction = bond_side(p1, p2, scene_pos)
            return self._regular_polygon_from_edge(p1, p2, ring_size, direction)

        anchor = self.model.get_atom(self._drag_anchor["id"])
        p1 = QPointF(anchor.x, anchor.y)
        p2 = self._compute_bond_endpoint(
            p1,
            scene_pos,
            modifiers,
            bond_order=1,
            is_aromatic=self.state.active_ring_aromatic,
        )
        direction = bond_side(p1, p2, scene_pos)
        return self._regular_polygon_from_edge(p1, p2, ring_size, direction)

    def _regular_polygon_from_edge(
        self,
        p1: QPointF,
        p2: QPointF,
        ring_size: int,
        direction: int,
    ) -> List[QPointF]:
        dx = p2.x() - p1.x()
        dy = p2.y() - p1.y()
        length = math.hypot(dx, dy)
        if length < 1e-6:
            return []

        step = 360.0 / ring_size
        angle0 = angle_deg(p1, p2)
        vertices = [p1, p2]
        current = p2
        current_angle = angle0
        for _ in range(ring_size - 2):
            current_angle = (current_angle + direction * step) % 360.0
            next_point = endpoint_from_angle_len(current, current_angle, length)
            vertices.append(next_point)
            current = next_point
        return vertices

    # -------------------------------------------------------------------------
    # Chain Placement
    # -------------------------------------------------------------------------
    def _begin_place_chain(self, anchor_id: int, scene_pos: QPointF) -> None:
        self._drag_mode = "place_chain"
        anchor = self.model.get_atom(anchor_id)
        self._drag_anchor = {"type": "atom", "id": anchor_id, "pos": QPointF(anchor.x, anchor.y)}
        radius = self.state.bond_length * OPTIMIZE_ZONE_SCALE
        self._optimize_zone.set_radius(radius)
        self._optimize_zone.update_center(anchor.x, anchor.y)
        self._update_chain_preview(scene_pos, Qt.KeyboardModifier.NoModifier)

    def _update_chain_preview(self, scene_pos: QPointF, modifiers: Qt.KeyboardModifiers) -> None:
        if self._drag_anchor is None:
            return
        anchor = self.model.get_atom(self._drag_anchor["id"])
        p0 = QPointF(anchor.x, anchor.y)
        dx = scene_pos.x() - p0.x()
        dy = scene_pos.y() - p0.y()
        dist = math.hypot(dx, dy)

        if dist < 1.0:
            points = [p0]
            self._chain_last_points = points
            self._preview_chain_item.update_polyline(points)
            self._preview_chain_label.hide_label()
            return

        raw_angle = angle_deg(p0, scene_pos)
        use_alt = bool(modifiers & Qt.KeyboardModifier.AltModifier)
        use_shift = bool(modifiers & Qt.KeyboardModifier.ShiftModifier)
        use_optimize = False
        if self._optimize_zone.isVisible():
            use_optimize = dist <= self._optimize_zone.radius()

        if use_shift:
            if abs(dx) >= abs(dy):
                raw_angle = 0.0 if dx >= 0 else 180.0
            else:
                raw_angle = 90.0 if -(dy) >= 0 else 270.0

        if self.state.fixed_angles and not use_alt:
            if use_optimize:
                raw_angle = self._select_preferred_angle(self._drag_anchor["id"], raw_angle, 1, False)
            raw_angle, _ = self._pick_bond_direction_deg(
                p0,
                self._drag_anchor["id"],
                raw_angle,
                1,
                False,
                self.state.bond_length,
                apply_collisions=True,
                allow_length_boost=False,
            )

        n = max(1, min(CHAIN_MAX_BONDS, int(round(dist / self.state.bond_length))))
        points = [p0]
        geometry = self._bond_geometry(self._drag_anchor["id"], 1, False)
        zigzag = geometry == "sp3" and not use_alt and self.state.fixed_angles
        current = p0
        for i in range(1, n + 1):
            if zigzag:
                step_angle = raw_angle + (60.0 if i % 2 == 1 else -60.0)
            else:
                step_angle = raw_angle
            next_point = endpoint_from_angle_len(current, step_angle, self.state.bond_length)
            points.append(next_point)
            current = next_point

        self._chain_last_points = points
        self._preview_chain_item.update_polyline(points)
        self._preview_chain_label.update_label(str(n), points[-1] + QPointF(0, -10))

    def _finalize_chain(self, modifiers: Qt.KeyboardModifiers) -> None:
        if not self._chain_last_points or self._drag_anchor is None:
            self._cancel_drag()
            return
        anchor_id = self._drag_anchor["id"]
        new_positions = [(p.x(), p.y()) for p in self._chain_last_points[1:]]
        if not new_positions:
            self._cancel_drag()
            return
        cmd = AddChainCommand(
            self.model,
            self,
            anchor_id,
            new_positions,
            element=self.state.default_element,
        )
        self.undo_stack.push(cmd)
        self._cancel_drag()

    # -------------------------------------------------------------------------
    # Hotkeys
    # -------------------------------------------------------------------------
    def _handle_hotkeys(self, event) -> bool:
        key = event.key()
        if Qt.Key.Key_0 <= key <= Qt.Key.Key_9:
            digit = key - Qt.Key.Key_0
            if digit == 0:
                return False
            if self.hovered_atom_id is not None:
                if event.modifiers() & Qt.KeyboardModifier.ShiftModifier:
                    if digit >= 3:
                        self._create_ring_hotkey(self.hovered_atom_id, digit)
                        return True
                else:
                    self._create_chain_hotkey(self.hovered_atom_id, digit)
                    return True
            if self.hovered_bond_id is not None:
                if event.modifiers() & Qt.KeyboardModifier.ShiftModifier:
                    if digit >= 3:
                        self._create_ring_from_bond_hotkey(self.hovered_bond_id, digit)
                        return True
                else:
                    self._set_bond_order_hotkey(self.hovered_bond_id, digit)
                    return True
        return False

    def _create_chain_hotkey(self, anchor_id: int, n: int) -> None:
        anchor = self.model.get_atom(anchor_id)
        p0 = QPointF(anchor.x, anchor.y)
        angles = self._get_anchor_bond_angles_deg(anchor_id)
        theta = choose_optimal_direction(angles)
        positions = []
        for i in range(1, n + 1):
            p = endpoint_from_angle_len(p0, theta, self.state.bond_length * i)
            positions.append((p.x(), p.y()))
        cmd = AddChainCommand(self.model, self, anchor_id, positions, self.state.default_element)
        self.undo_stack.push(cmd)

    def _create_ring_hotkey(self, anchor_id: int, ring_size: int) -> None:
        anchor = self.model.get_atom(anchor_id)
        p0 = QPointF(anchor.x, anchor.y)
        angles = self._get_anchor_bond_angles_deg(anchor_id)
        theta = choose_optimal_direction(angles)
        p1 = endpoint_from_angle_len(p0, theta, self.state.bond_length)
        direction = bond_side(p0, p1, self._last_scene_pos)
        vertices = self._regular_polygon_from_edge(p0, p1, ring_size, direction)
        self._commit_ring_vertices(vertices, anchor_type="atom", anchor_id=anchor_id)

    def _create_ring_from_bond_hotkey(self, bond_id: int, ring_size: int) -> None:
        bond = self.model.get_bond(bond_id)
        a1 = self.model.get_atom(bond.a1_id)
        a2 = self.model.get_atom(bond.a2_id)
        p1 = QPointF(a1.x, a1.y)
        p2 = QPointF(a2.x, a2.y)
        direction = bond_side(p1, p2, self._last_scene_pos)
        vertices = self._regular_polygon_from_edge(p1, p2, ring_size, direction)
        self._commit_ring_vertices(vertices, anchor_type="bond", anchor_id=bond_id)

    def _commit_ring_vertices(
        self, vertices: List[QPointF], anchor_type: str, anchor_id: int
    ) -> None:
        if not vertices:
            return
        ring_size = len(vertices)
        aromatic = self.state.active_ring_aromatic
        vertex_defs = self._build_ring_vertex_defs(vertices, anchor_type, anchor_id)

        edge_defs: List[Tuple[int, int, int, BondStyle, BondStereo, bool]] = []
        if aromatic:
            edge_defs = self._build_aromatic_edges(vertex_defs, ring_size)
        else:
            for i in range(ring_size):
                j = (i + 1) % ring_size
                edge_defs.append((i, j, 1, BondStyle.PLAIN, BondStereo.NONE, aromatic))

        cmd = AddRingCommand(
            self.model,
            self,
            vertex_defs,
            edge_defs,
            element=self.state.default_element,
        )
        self.undo_stack.push(cmd)
        if aromatic:
            self._kekulize_aromatic_bonds()

    def _set_bond_order_hotkey(self, bond_id: int, order: int) -> None:
        order = max(1, min(3, order))
        cmd = ChangeBondCommand(
            self.model,
            self,
            bond_id,
            new_order=order,
            new_style=BondStyle.PLAIN,
            new_stereo=BondStereo.NONE,
            new_is_aromatic=False,
        )
        self.undo_stack.push(cmd)

    def _kekulize_aromatic_bonds(self, seed_atoms: Optional[Iterable[int]] = None) -> None:
        aromatic_bonds = [bond for bond in self.model.bonds.values() if bond.is_aromatic]
        if not aromatic_bonds:
            return

        component_atoms: Optional[set[int]] = None
        if seed_atoms is not None:
            adjacency: dict[int, list[int]] = {}
            for bond in aromatic_bonds:
                adjacency.setdefault(bond.a1_id, []).append(bond.a2_id)
                adjacency.setdefault(bond.a2_id, []).append(bond.a1_id)
            seeds = [atom_id for atom_id in seed_atoms if atom_id in adjacency]
            if not seeds:
                return
            component_atoms = set()
            stack = list(seeds)
            while stack:
                node = stack.pop()
                if node in component_atoms:
                    continue
                component_atoms.add(node)
                for neighbor in adjacency.get(node, []):
                    if neighbor not in component_atoms:
                        stack.append(neighbor)

        display_orders = kekulize_display_orders(self.model, seed_atoms=component_atoms)

        if display_orders is None:
            try:
                import networkx as nx
            except Exception:
                nx = None

            if nx is not None:
                graph = nx.Graph()
                bond_by_pair: dict[frozenset[int], int] = {}
                for bond in aromatic_bonds:
                    if component_atoms is not None and (
                        bond.a1_id not in component_atoms or bond.a2_id not in component_atoms
                    ):
                        continue
                    pair = frozenset({bond.a1_id, bond.a2_id})
                    bond_by_pair[pair] = bond.id
                    weight = 10_000_000 - (min(pair) * 1000 + max(pair))
                    graph.add_edge(bond.a1_id, bond.a2_id, weight=weight)
                matching = nx.max_weight_matching(graph, maxcardinality=True, weight="weight")
                display_orders = {bond_id: 1 for bond_id in bond_by_pair.values()}
                for u, v in matching:
                    bond_id = bond_by_pair.get(frozenset({u, v}))
                    if bond_id is not None:
                        display_orders[bond_id] = 2

        if display_orders is None:
            display_orders = {}
            used_atoms: set[int] = set()
            sorted_bonds = sorted(
                aromatic_bonds,
                key=lambda bond: (min(bond.a1_id, bond.a2_id), max(bond.a1_id, bond.a2_id)),
            )
            for bond in sorted_bonds:
                if component_atoms is not None and (
                    bond.a1_id not in component_atoms or bond.a2_id not in component_atoms
                ):
                    continue
                display_orders[bond.id] = 1
            for bond in sorted_bonds:
                if component_atoms is not None and (
                    bond.a1_id not in component_atoms or bond.a2_id not in component_atoms
                ):
                    continue
                if bond.a1_id in used_atoms or bond.a2_id in used_atoms:
                    continue
                display_orders[bond.id] = 2
                used_atoms.add(bond.a1_id)
                used_atoms.add(bond.a2_id)

        external_double_atoms: set[int] = set()
        for bond in self.model.bonds.values():
            if bond.is_aromatic:
                continue
            if bond.order >= 2:
                external_double_atoms.add(bond.a1_id)
                external_double_atoms.add(bond.a2_id)

        for bond in aromatic_bonds:
            if component_atoms is not None and (
                bond.a1_id not in component_atoms or bond.a2_id not in component_atoms
            ):
                continue
            new_order = display_orders.get(bond.id, 1)
            if bond.a1_id in external_double_atoms or bond.a2_id in external_double_atoms:
                new_order = 1
            if bond.display_order != new_order or bond.order != new_order:
                bond.display_order = new_order
                bond.order = new_order
                self.update_bond_item(bond.id)

    def _build_aromatic_edges(
        self,
        vertex_defs: List[Tuple[Optional[int], float, float]],
        ring_size: int,
    ) -> List[Tuple[int, int, int, BondStyle, BondStereo, bool]]:
        if ring_size % 2 != 0:
            return [
                (i, (i + 1) % ring_size, 1, BondStyle.PLAIN, BondStereo.NONE, True)
                for i in range(ring_size)
            ]

        constraints: dict[int, int] = {}
        for i in range(ring_size):
            j = (i + 1) % ring_size
            a_id = vertex_defs[i][0]
            b_id = vertex_defs[j][0]
            if a_id is None or b_id is None:
                continue
            existing = self.model.find_bond_between(a_id, b_id)
            if existing is not None:
                constraints[i] = 2 if existing.order >= 2 else 1

        def atom_has_double(atom_id: Optional[int], other_id: Optional[int]) -> bool:
            if atom_id is None:
                return False
            for bond in self.model.bonds.values():
                if bond.order < 2:
                    continue
                if bond.a1_id == atom_id and bond.a2_id != other_id:
                    return True
                if bond.a2_id == atom_id and bond.a1_id != other_id:
                    return True
            return False

        def score_phase(phase: int) -> int:
            score = 0
            for idx, desired in constraints.items():
                current = 2 if ((idx + phase) % 2 == 0) else 1
                if current != desired:
                    score += 1
            return score

        phase = 0 if score_phase(0) <= score_phase(1) else 1
        edges: List[Tuple[int, int, int, BondStyle, BondStereo, bool]] = []
        for i in range(ring_size):
            j = (i + 1) % ring_size
            order = 2 if ((i + phase) % 2 == 0) else 1
            a_id = vertex_defs[i][0]
            b_id = vertex_defs[j][0]
            if order == 2 and (atom_has_double(a_id, b_id) or atom_has_double(b_id, a_id)):
                order = 1
            edges.append((i, j, order, BondStyle.PLAIN, BondStereo.NONE, True))
        return edges

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

    def _is_on_paper(self, x: float, y: float) -> bool:
        return (
            PAPER_MARGIN <= x <= PAPER_WIDTH - PAPER_MARGIN
            and PAPER_MARGIN <= y <= PAPER_HEIGHT - PAPER_MARGIN
        )
