"""
Chemuson Canvas
Page-based canvas using QGraphicsView/QGraphicsScene for document-style editing.
"""
from __future__ import annotations

import math
try:
    from PyQt6 import sip
except ImportError:
    import sip
from typing import Dict, List, Optional, Tuple

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
    snap_angle_deg,
    endpoint_from_angle_len,
    choose_optimal_direction,
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
    get_aromatic_ring_geometries,
    get_visual_properties,
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
HOVER_ATOM_RADIUS = 20.0
HOVER_BOND_DISTANCE = 14.0
HOVER_BOND_DISTANCE = 10.0
OPTIMIZE_ZONE_SCALE = 1.2
CHAIN_MAX_BONDS = 12


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

    def _ensure_overlays(self) -> None:
        """Check if any overlays have been deleted by C++ and recreate them if so."""
        if not self._overlays_ready:
            return
            
        needed_recreation = False
        
        # Helper to check and recreate a single overlay item
        def check_item(attr_name, class_type):
            nonlocal needed_recreation
            item = getattr(self, attr_name, None)
            if item is None or sip.isdeleted(item):
                new_item = class_type()
                setattr(self, attr_name, new_item)
                self.scene.addItem(new_item)
                needed_recreation = True
                
        check_item("_hover_atom_indicator", HoverAtomIndicatorItem)
        check_item("_hover_bond_indicator", HoverBondIndicatorItem)
        check_item("_optimize_zone", OptimizeZoneItem)
        check_item("_preview_bond_item", PreviewBondItem)
        check_item("_preview_ring_item", PreviewRingItem)
        check_item("_preview_chain_item", PreviewChainItem)
        check_item("_preview_chain_label", PreviewChainLabelItem)
        
        if needed_recreation:
            # If we recreated anything, ensure they have proper Z values or other state
            # but usually the __init__ of the items handles defaults.
            pass

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
                self._update_ring_centers()
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
        item.set_offset_sign(self._bond_offset_sign(bond))
        self.scene.addItem(item)
        self.bond_items[bond.id] = item

    def update_bond_item(self, bond_id: int) -> None:
        bond = self.model.get_bond(bond_id)
        atom1 = self.model.get_atom(bond.a1_id)
        atom2 = self.model.get_atom(bond.a2_id)
        item = self.bond_items.get(bond_id)
        if item is not None:
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

    def _update_ring_centers(self) -> None:
        self._update_visual_properties_from_rdkit()

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
        # Use common RDKit visual update
        self._update_visual_properties_from_rdkit()

    def _update_visual_properties_from_rdkit(self) -> None:
        """
        Update visual bond orders and ring centers using RDKit logic.
        Also manages aromatic circles if enabled.
        """
        try:
            orders, centers = get_visual_properties(self.model)
            
            # 1. Update Bond Orders and Centers
            use_circles = self.state.use_aromatic_circles
            for bond_id, item in self.bond_items.items():
                # Bond Order (Kekulization)
                if use_circles:
                     # If circles are on, we usually don't force double bonds on aromatic rings,
                     # we let them be single (or aromatic specific style).
                     # But RDKit returns '2' for double in Kekule.
                     # If we want circle, we might want to RESET visual order to None (use model)
                     # or use '1' if we want strictly single lines under the circle?
                     # ChemDraw style: Circle usually implies "don't show alternating".
                     item.set_visual_order(None) 
                else:
                     # Apply Kekule order
                     item.set_visual_order(orders.get(bond_id))
                
                # Ring Center (for 'inside' normal calculation)
                # Even with circles, we might want this if we have double bonds?
                # But mostly for non-circle mode.
                center_tuple = centers.get(bond_id)
                if center_tuple is not None:
                    item.set_ring_context(QPointF(*center_tuple))
                else:
                    item.set_ring_context(None)
                    
                self.update_bond_item(bond_id)

            # 2. Update Aromatic Circles
            # Clear existing
            for circle in self.aromatic_circles:
                self.scene.removeItem(circle)
            self.aromatic_circles.clear()

            if use_circles:
                for cx, cy, radius in get_aromatic_ring_geometries(self.model):
                    circle = AromaticCircleItem(cx, cy, radius)
                    self.scene.addItem(circle)
                    self.aromatic_circles.append(circle)
                    
        except Exception:
            # Fallback safe cleanup
            pass

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

    def _add_bond_between(self, a1_id: int, a2_id: int) -> None:
        order = self.state.active_bond_order
        style = self.state.active_bond_style
        stereo = self.state.active_bond_stereo
        is_aromatic = self.state.active_bond_aromatic
        if style != BondStyle.PLAIN:
            order = 1
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
        self._ensure_overlays()
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
            self._ensure_overlays()
            try:
                self._optimize_zone.hide_zone()
                self._preview_bond_item.hide_preview()
                self._preview_ring_item.hide_preview()
                self._preview_chain_item.hide_preview()
                self._preview_chain_label.hide_label()
            except (RuntimeError, AttributeError):
                pass

    def _on_undo_stack_changed(self, _index: int) -> None:
        self._cancel_drag()
        self._update_hover(self._last_scene_pos)
        # Refresh visual properties (Kekulization/Rings) after any model change
        self._update_visual_properties_from_rdkit()

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

    def _anchor_geometry(self, anchor_id: int) -> str:
        for bond in self.model.bonds.values():
            if bond.a1_id == anchor_id or bond.a2_id == anchor_id:
                if bond.order >= 3:
                    return "sp"
                if bond.order == 2 or bond.is_aromatic:
                    return "sp2"
        return "sp3"

    def _preferred_angles(self, anchor_id: int) -> list[float]:
        geometry = self._anchor_geometry(anchor_id)
        existing = self._get_anchor_bond_angles_deg(anchor_id)
        if not existing:
            return []
        candidates: list[float] = []
        if geometry == "sp":
            for angle in existing:
                candidates.append((angle + 180.0) % 360.0)
        else:
            for angle in existing:
                candidates.append((angle + 120.0) % 360.0)
                candidates.append((angle - 120.0) % 360.0)
        return candidates

    def _select_preferred_angle(self, anchor_id: int, cursor_angle: float) -> float:
        candidates = self._preferred_angles(anchor_id)
        if not candidates:
            return cursor_angle
        def score(angle: float) -> float:
            diff = abs(((angle - cursor_angle + 180.0) % 360.0) - 180.0)
            return diff
        return min(candidates, key=score)

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

        target_id = closest_atom(
            p1,
            [(atom.id, atom.x, atom.y) for atom in self.model.atoms.values()],
            ATOM_HIT_RADIUS,
        )
        if target_id == anchor_id:
            target_id = None

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
                cmd = AddBondCommand(
                    self.model,
                    self,
                    anchor_id,
                    target_id,
                    order,
                    style,
                    stereo,
                    is_aromatic=is_aromatic,
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
            self.undo_stack.endMacro()
        else:
            if target_id is not None:
                cmd = AddBondCommand(
                    self.model,
                    self,
                    anchor_id,
                    target_id,
                    order,
                    style,
                    stereo,
                    is_aromatic=is_aromatic,
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
        self._cancel_drag()

    def _compute_bond_endpoint(
        self, anchor: QPointF, scene_pos: QPointF, modifiers: Qt.KeyboardModifiers
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
        if use_optimize and not use_alt and self._drag_anchor["id"] is not None:
            theta = self._select_preferred_angle(self._drag_anchor["id"], cursor_theta)
        else:
            theta = cursor_theta

        if self.state.fixed_angles and not use_alt:
            theta = snap_angle_deg(theta, self.state.angle_step_deg)

        if self.state.fixed_lengths and not use_shift:
            length = self.state.bond_length
        else:
            length = max(5.0, dist)

        return endpoint_from_angle_len(anchor, theta, length)

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

        vertex_defs: List[Tuple[Optional[int], float, float]] = []
        if self._drag_anchor["type"] == "bond":
            bond = self.model.get_bond(self._drag_anchor["id"])
            vertex_defs.append((bond.a1_id, vertices[0].x(), vertices[0].y()))
            vertex_defs.append((bond.a2_id, vertices[1].x(), vertices[1].y()))
            for v in vertices[2:]:
                vertex_defs.append((None, v.x(), v.y()))
        elif self._drag_anchor["type"] == "free":
            for v in vertices:
                vertex_defs.append((None, v.x(), v.y()))
        else:
            vertex_defs.append((self._drag_anchor["id"], vertices[0].x(), vertices[0].y()))
            for v in vertices[1:]:
                vertex_defs.append((None, v.x(), v.y()))

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
        p2 = self._compute_bond_endpoint(p1, scene_pos, modifiers)
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

        if use_optimize and not use_alt:
            raw_angle = self._select_preferred_angle(self._drag_anchor["id"], raw_angle)

        if use_shift:
            if abs(dx) >= abs(dy):
                raw_angle = 0.0 if dx >= 0 else 180.0
            else:
                raw_angle = 90.0 if -(dy) >= 0 else 270.0

        if self.state.fixed_angles and not use_alt:
            raw_angle = snap_angle_deg(raw_angle, self.state.angle_step_deg)

        n = max(1, min(CHAIN_MAX_BONDS, int(round(dist / self.state.bond_length))))
        points = [p0]
        geometry = self._anchor_geometry(self._drag_anchor["id"])
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
        vertex_defs: List[Tuple[Optional[int], float, float]] = []
        if anchor_type == "bond":
            bond = self.model.get_bond(anchor_id)
            vertex_defs.append((bond.a1_id, vertices[0].x(), vertices[0].y()))
            vertex_defs.append((bond.a2_id, vertices[1].x(), vertices[1].y()))
            for v in vertices[2:]:
                vertex_defs.append((None, v.x(), v.y()))
        else:
            vertex_defs.append((anchor_id, vertices[0].x(), vertices[0].y()))
            for v in vertices[1:]:
                vertex_defs.append((None, v.x(), v.y()))

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

    def _kekulize_aromatic_bonds(self) -> None:
        """Deprecated: RDKit-based visuals handle Kekulization."""
        aromatic_bonds = [bond for bond in self.model.bonds.values() if bond.is_aromatic]
        if not aromatic_bonds:
            return

        adjacency: dict[int, list[int]] = {}
        for bond in aromatic_bonds:
            adjacency.setdefault(bond.a1_id, []).append(bond.a2_id)
            adjacency.setdefault(bond.a2_id, []).append(bond.a1_id)

        color: dict[int, int] = {}
        for node in adjacency:
            if node in color:
                continue
            color[node] = 0
            stack = [node]
            while stack:
                current = stack.pop()
                for neighbor in adjacency.get(current, []):
                    if neighbor in color:
                        if color[neighbor] == color[current]:
                            return
                    else:
                        color[neighbor] = 1 - color[current]
                        stack.append(neighbor)

        u_nodes = [n for n, c in color.items() if c == 0]
        v_nodes = [n for n, c in color.items() if c == 1]
        u_index = {u: i for i, u in enumerate(u_nodes)}
        v_index = {v: i for i, v in enumerate(v_nodes)}

        degrees = {node: len(adjacency.get(node, [])) for node in adjacency}

        def edge_weight(u: int, v: int) -> int:
            du = degrees.get(u, 0)
            dv = degrees.get(v, 0)
            if du == 2 and dv == 2:
                return 3
            if du == 2 or dv == 2:
                return 2
            return 1

        class Edge:
            def __init__(self, to: int, rev: int, cap: int, cost: int) -> None:
                self.to = to
                self.rev = rev
                self.cap = cap
                self.cost = cost

        size = 2 + len(u_nodes) + len(v_nodes)
        src = 0
        sink = size - 1
        graph: list[list[Edge]] = [[] for _ in range(size)]

        def add_edge(frm: int, to: int, cap: int, cost: int) -> None:
            graph[frm].append(Edge(to, len(graph[to]), cap, cost))
            graph[to].append(Edge(frm, len(graph[frm]) - 1, 0, -cost))

        for u in u_nodes:
            add_edge(src, 1 + u_index[u], 1, 0)
        for v in v_nodes:
            add_edge(1 + len(u_nodes) + v_index[v], sink, 1, 0)

        for u in u_nodes:
            for v in adjacency.get(u, []):
                if v not in v_index:
                    continue
                w = edge_weight(u, v)
                add_edge(1 + u_index[u], 1 + len(u_nodes) + v_index[v], 1, -w)

        flow = 0
        while True:
            dist = [10 ** 9] * size
            in_queue = [False] * size
            prev_node = [-1] * size
            prev_edge = [-1] * size
            dist[src] = 0
            queue = [src]
            in_queue[src] = True
            while queue:
                node = queue.pop(0)
                in_queue[node] = False
                for i, edge in enumerate(graph[node]):
                    if edge.cap <= 0:
                        continue
                    nd = dist[node] + edge.cost
                    if nd < dist[edge.to]:
                        dist[edge.to] = nd
                        prev_node[edge.to] = node
                        prev_edge[edge.to] = i
                        if not in_queue[edge.to]:
                            in_queue[edge.to] = True
                            queue.append(edge.to)
            if dist[sink] == 10 ** 9:
                break
            v = sink
            while v != src:
                edge = graph[prev_node[v]][prev_edge[v]]
                edge.cap -= 1
                graph[v][edge.rev].cap += 1
                v = prev_node[v]
            flow += 1

        pair_u: dict[int, int | None] = {u: None for u in u_nodes}
        for u in u_nodes:
            node_u = 1 + u_index[u]
            for edge in graph[node_u]:
                if edge.to == src or edge.cap != 0:
                    continue
                if edge.to == sink:
                    continue
                if edge.to >= 1 + len(u_nodes) and edge.to < sink:
                    v = v_nodes[edge.to - 1 - len(u_nodes)]
                    pair_u[u] = v
                    break

        for bond in aromatic_bonds:
            double = False
            if color.get(bond.a1_id) == 0:
                double = pair_u.get(bond.a1_id) == bond.a2_id
            else:
                double = pair_u.get(bond.a2_id) == bond.a1_id
            order = 2 if double else 1
            if bond.order != order:
                self.model.update_bond(bond.id, order=order)
                self.update_bond_item(bond.id)

    def _build_aromatic_edges(
        self,
        vertex_defs: List[Tuple[Optional[int], float, float]],
        ring_size: int,
    ) -> List[Tuple[int, int, int, BondStyle, BondStereo, bool]]:
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

        edges: List[Tuple[int, int, int, BondStyle, BondStereo, bool]] = []
        for i in range(ring_size):
            j = (i + 1) % ring_size
            order = constraints.get(i, 1)
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
