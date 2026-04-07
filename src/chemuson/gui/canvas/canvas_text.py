from __future__ import annotations

import math
from typing import Iterable, Optional

from PyQt6.QtCore import QPointF, Qt
from PyQt6.QtGui import QBrush, QColor, QFont, QPen, QTextCharFormat, QTextCursor, QTextOption
from PyQt6.QtWidgets import QGraphicsLineItem, QGraphicsTextItem, QInputDialog

from chemuson.core.model import SIMPLE_HYDROGEN_GROUP_LABELS, BondStyle, bond_is_structural
from chemuson.gui.commands import (
    AddTextItemCommand,
    AddWavyAnchorCommand,
    ChangeAtomLabelScaleCommand,
    FormatTextItemsCommand,
)
from chemuson.gui.dialogs import AtomLabelDialog
from chemuson.gui.geom import angle_deg, angle_distance_deg, endpoint_from_angle_len, snap_angle_deg
from chemuson.gui.items import (
    ABBREVIATION_LABELS,
    ArrowItem,
    AtomItem,
    TextAnnotationItem,
    WavyAnchorItem,
)

from .canvas_chem_data import (
    CANONICAL_SIMPLE_HYDROGEN_GROUP_LABELS,
    ELECTRON_SLOT_ANGLES,
    ELEMENT_SYMBOLS,
    FUNCTIONAL_GROUP_ALIASES,
    FUNCTIONAL_GROUP_LABELS,
    HETERO_ELECTRON_ATOMS,
    IMPLICIT_H_ELEMENTS,
)
from .canvas_constants import (
    ELECTRON_ANCHOR_ROLE,
    ELECTRON_DOT_ROLE,
    ELECTRON_SCALE_ROLE,
    ELECTRON_SIDE_ROLE,
    ELECTRON_SLOT_ROLE,
    ELECTRON_SLOT_TOLERANCE_DEG,
    IMPLICIT_H_OVERLAY_ANCHOR_ROLE,
    IMPLICIT_H_OVERLAY_ANGLE_ROLE,
    LABEL_OFFSET_MIN_PX,
    LABEL_OFFSET_SCALE,
    WAVY_ANCHOR_ANGLE_ROLE,
    WAVY_ANCHOR_BOND_ROLE,
    WAVY_ANCHOR_LENGTH_ROLE,
    WAVY_ANCHOR_ROLE,
)

class CanvasTextMixin:
    def _active_text_edit_item(self) -> TextAnnotationItem | None:
        """Devuelve el item de texto actualmente en edición, si existe."""
        focus_item = self.scene.focusItem()
        if (
            isinstance(focus_item, TextAnnotationItem)
            and focus_item.scene() is self.scene
            and (
                focus_item.textInteractionFlags()
                & Qt.TextInteractionFlag.TextEditorInteraction
            )
        ):
            return focus_item
        last_item = getattr(self, "_last_text_edit_item", None)
        if (
            isinstance(last_item, TextAnnotationItem)
            and last_item.scene() is self.scene
            and (
                last_item.textInteractionFlags()
                & Qt.TextInteractionFlag.TextEditorInteraction
            )
        ):
            return last_item
        return None

    def _clear_text_cursor_selection(self, item: TextAnnotationItem) -> bool:
        """Limpia una selección interna de QTextCursor si quedó visible."""
        if item is None or item.scene() is not self.scene:
            return False
        cursor = item.textCursor()
        if not cursor.hasSelection():
            return False
        cursor.clearSelection()
        item.setTextCursor(cursor)
        return True

    def _text_item_format_snapshot(self, item: TextAnnotationItem) -> dict:
        """Captura un snapshot del documento y cursor para undo/redo."""
        cursor = item.textCursor()
        return {
            "html": item.document().toHtml(),
            "font": item.font().toString(),
            "color": item.defaultTextColor().name(QColor.NameFormat.HexArgb),
            "text_width": float(item.textWidth()),
            "editing": bool(
                item.textInteractionFlags()
                & Qt.TextInteractionFlag.TextEditorInteraction
            ),
            "selected": bool(item.isSelected()),
            "cursor_position": int(cursor.position()),
            "cursor_anchor": int(cursor.anchor()),
        }

    def _restore_text_item_format_snapshot(self, item: TextAnnotationItem, snapshot: dict) -> None:
        """Restaura documento/cursor de un texto desde un snapshot."""
        if item.scene() is not self.scene:
            return
        item.document().setHtml(snapshot.get("html", item.document().toHtml()))
        font = QFont()
        font_str = snapshot.get("font")
        if font_str:
            font.fromString(font_str)
            item.setFont(font)
        color_name = snapshot.get("color")
        if color_name:
            item.setDefaultTextColor(QColor(color_name))
        item.setTextWidth(float(snapshot.get("text_width", item.textWidth())))
        is_editing = bool(snapshot.get("editing", False))
        item.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextEditorInteraction
            if is_editing
            else Qt.TextInteractionFlag.NoTextInteraction
        )
        item.setSelected(bool(snapshot.get("selected", False)))

        text_len = len(item.toPlainText())
        anchor = max(0, min(int(snapshot.get("cursor_anchor", 0)), text_len))
        position = max(0, min(int(snapshot.get("cursor_position", anchor)), text_len))
        cursor = item.textCursor()
        cursor.setPosition(anchor)
        move_mode = (
            QTextCursor.MoveMode.KeepAnchor
            if position != anchor
            else QTextCursor.MoveMode.MoveAnchor
        )
        cursor.setPosition(position, move_mode)
        item.setTextCursor(cursor)
        if is_editing:
            item.setFocus()
        else:
            self._clear_text_cursor_selection(item)
            if self.scene.focusItem() is item:
                self.scene.clearFocus()

    def _text_items_for_formatting(self) -> list[TextAnnotationItem]:
        """Devuelve los textos afectados por una acción de formato."""
        items_to_update = [
            item for item in self.scene.selectedItems() if isinstance(item, TextAnnotationItem)
        ]
        focus_item = self.scene.focusItem()
        if isinstance(focus_item, TextAnnotationItem) and focus_item not in items_to_update:
            items_to_update.append(focus_item)
        return items_to_update

    def sync_text_selection_state(self) -> None:
        """Refresca la UI asociada al estado actual del cursor de texto."""
        self._sync_selection_from_scene()

    def _finish_active_text_editing(self, *, clear_scene_selection: bool = False) -> bool:
        """Finaliza la edición activa de texto y elimina resaltados residuales."""
        changed = False
        for selected_item in list(self.scene.selectedItems()):
            if isinstance(selected_item, TextAnnotationItem):
                changed = self._clear_text_cursor_selection(selected_item) or changed

        item = self._active_text_edit_item()
        if item is None:
            if clear_scene_selection and self.scene.selectedItems():
                self.scene.clearSelection()
                self._sync_selection_from_scene()
                return True
            if changed:
                self._sync_selection_from_scene()
            return changed

        changed = True
        item.setTextInteractionFlags(Qt.TextInteractionFlag.NoTextInteraction)
        self._clear_text_cursor_selection(item)
        if self.scene.focusItem() is item:
            self.scene.clearFocus()
        item.clearFocus()
        self.remember_text_edit_item(None)

        if item.scene() is self.scene and not item.toPlainText().strip():
            self.remove_text_item(item)
        if clear_scene_selection:
            self.scene.clearSelection()
        self._sync_selection_from_scene()
        return changed

    @staticmethod
    def _curve_factor_from_pointer(start: QPointF, end: QPointF, pointer: QPointF) -> float:
        """Calcula curvatura normalizada a partir de la distancia perpendicular."""
        dx = end.x() - start.x()
        dy = end.y() - start.y()
        length = math.hypot(dx, dy)
        if length <= 1e-6:
            return ArrowItem.CURVE_FACTOR_DEFAULT
        nx = -dy / length
        ny = dx / length
        mid_x = (start.x() + end.x()) * 0.5
        mid_y = (start.y() + end.y()) * 0.5
        offset = (pointer.x() - mid_x) * nx + (pointer.y() - mid_y) * ny
        # Higher sensitivity so users can reach strong curvature comfortably.
        factor = offset / (length * 0.65)
        return max(ArrowItem.CURVE_FACTOR_MIN, min(ArrowItem.CURVE_FACTOR_MAX, factor))

    def _constrain_annotation_endpoint(
        self,
        start: QPointF,
        end: QPointF,
        modifiers: Qt.KeyboardModifiers,
    ) -> QPointF:
        """Restringe segmentos simples a incrementos de 45° cuando se usa Shift."""
        if self.current_tool not in {"tool_arrow_line", "tool_arrow_line_dashed"}:
            return QPointF(end)
        if not (modifiers & Qt.KeyboardModifier.ShiftModifier):
            return QPointF(end)
        dx = end.x() - start.x()
        dy = end.y() - start.y()
        length = math.hypot(dx, dy)
        if length <= 1e-6:
            return QPointF(end)
        snapped_angle = snap_angle_deg(angle_deg(start, end), 45.0)
        return endpoint_from_angle_len(start, snapped_angle, length)

    def _symbol_insert_position(self, scene_pos: QPointF, atom_id: Optional[int]) -> QPointF:
        """Método auxiliar para  symbol insert position.

        Args:
            scene_pos: Descripción del parámetro.
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if atom_id is None:
            return scene_pos
        atom = self.model.get_atom(atom_id)
        center = QPointF(atom.x, atom.y)
        direction = scene_pos - center
        if direction.manhattanLength() < 2.0:
            direction = self._label_open_direction(atom_id)
        if direction.manhattanLength() < 1e-3:
            direction = QPointF(1.0, -1.0)
        length = math.hypot(direction.x(), direction.y())
        if length <= 1e-6:
            return scene_pos
        nx = direction.x() / length
        ny = direction.y() / length
        offset = max(12.0, self.state.bond_length * 0.3)
        return QPointF(center.x() + nx * offset, center.y() + ny * offset)

    def _symbol_anchor_data(
        self,
        atom_id: int,
        scene_pos: QPointF,
        offset: Optional[float] = None,
    ) -> tuple[QPointF, float]:
        """Método auxiliar para  symbol anchor data.

        Args:
            atom_id: Descripción del parámetro.
            scene_pos: Descripción del parámetro.
            offset: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atom = self.model.get_atom(atom_id)
        center = QPointF(atom.x, atom.y)
        direction = scene_pos - center
        if direction.manhattanLength() < 2.0:
            direction = self._label_open_direction(atom_id)
        if direction.manhattanLength() < 1e-3:
            direction = QPointF(1.0, -1.0)
        length = math.hypot(direction.x(), direction.y())
        if length <= 1e-6:
            return center, 0.0
        nx = direction.x() / length
        ny = direction.y() / length
        use_offset = max(12.0, self.state.bond_length * 0.3) if offset is None else offset
        pos = QPointF(center.x() + nx * use_offset, center.y() + ny * use_offset)
        angle = math.degrees(math.atan2(ny, nx))
        return pos, angle

    def _electron_dot_font_size(self, atom_id: int, scale: float) -> float:
        """Método auxiliar para  electron dot font size.

        Args:
            atom_id: Descripción del parámetro.
            scale: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = self.atom_items.get(atom_id)
        if item is not None:
            base_font = item.label.font()
        else:
            base_font = self._label_font()
        base_size = base_font.pointSizeF() or float(base_font.pointSize())
        return max(5.0, base_size * 0.55 * scale)

    def _electron_radial_offset(self, atom_id: int, scale: float) -> float:
        """Método auxiliar para  electron radial offset.

        Args:
            atom_id: Descripción del parámetro.
            scale: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = self.atom_items.get(atom_id)
        dot_size = self._electron_dot_font_size(atom_id, scale)
        if item is not None and item.label.isVisible():
            rect = item.label.boundingRect()
            label_radius = max(rect.width(), rect.height()) * 0.5
            return label_radius + dot_size * 0.45 + 2.0
        return max(14.0, self.state.bond_length * 0.25 + dot_size * 0.4)

    def _electron_pair_spread(self, atom_id: int, scale: float) -> float:
        """Método auxiliar para  electron pair spread.

        Args:
            atom_id: Descripción del parámetro.
            scale: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        dot_size = self._electron_dot_font_size(atom_id, scale)
        return max(3.0, dot_size * 0.55)

    def _wavy_anchor_length(self) -> float:
        """Método auxiliar para  wavy anchor length.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        return max(18.0, self.state.bond_length * 1.5)

    def _wavy_anchor_bond_angle(self, atom_id: int, scene_pos: QPointF) -> tuple[Optional[int], float]:
        """Método auxiliar para  wavy anchor bond angle.

        Args:
            atom_id: Descripción del parámetro.
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atom = self.model.get_atom(atom_id)
        center = QPointF(atom.x, atom.y)
        direction = scene_pos - center
        if direction.manhattanLength() < 2.0:
            direction = self._label_open_direction(atom_id)
        target_angle = None
        if direction.manhattanLength() >= 1e-3:
            target_angle = (
                math.degrees(math.atan2(direction.y(), direction.x())) + 360.0
            ) % 360.0

        bonds: list[tuple[int, float]] = []
        for bond in self.model.bonds.values():
            if bond.a1_id == atom_id:
                other = self.model.get_atom(bond.a2_id)
            elif bond.a2_id == atom_id:
                other = self.model.get_atom(bond.a1_id)
            else:
                continue
            dx = other.x - atom.x
            dy = other.y - atom.y
            angle = (math.degrees(math.atan2(dy, dx)) + 360.0) % 360.0
            bonds.append((bond.id, angle))

        if bonds:
            if target_angle is None:
                return bonds[0][0], bonds[0][1]
            best = min(bonds, key=lambda item: angle_distance_deg(item[1], target_angle))
            return best[0], best[1]

        if target_angle is None:
            return None, 0.0
        return None, target_angle

    def _position_wavy_anchor(self, item: WavyAnchorItem) -> None:
        """Método auxiliar para  position wavy anchor.

        Args:
            item: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        anchor_id = item.data(WAVY_ANCHOR_ROLE)
        base_angle = item.data(WAVY_ANCHOR_ANGLE_ROLE)
        bond_id = item.data(WAVY_ANCHOR_BOND_ROLE)
        length = item.data(WAVY_ANCHOR_LENGTH_ROLE)
        if anchor_id is None:
            return
        if anchor_id not in self.model.atoms:
            return
        atom = self.model.get_atom(anchor_id)
        use_length = float(length) if length is not None else self._wavy_anchor_length()

        angle_deg_value: Optional[float] = None
        if bond_id is not None and bond_id in self.model.bonds:
            bond = self.model.get_bond(int(bond_id))
            if bond.a1_id == anchor_id:
                other = self.model.get_atom(bond.a2_id)
            elif bond.a2_id == anchor_id:
                other = self.model.get_atom(bond.a1_id)
            else:
                other = None
            if other is not None:
                dx = other.x - atom.x
                dy = other.y - atom.y
                angle_deg_value = (math.degrees(math.atan2(dy, dx)) + 360.0) % 360.0

        if angle_deg_value is None and base_angle is not None:
            angle_deg_value = float(base_angle)
        if angle_deg_value is None:
            return

        perp_angle = math.radians(angle_deg_value + 90.0)
        nx = math.cos(perp_angle)
        ny = math.sin(perp_angle)
        center = QPointF(atom.x, atom.y)
        half = use_length * 0.5
        start = QPointF(center.x() - nx * half, center.y() - ny * half)
        end = QPointF(center.x() + nx * half, center.y() + ny * half)
        try:
            item.update_positions(start, end)
        except RuntimeError:
            self._wavy_anchors.discard(item)

    def _update_wavy_anchors_for_atom(self, atom_id: int) -> None:
        """Método auxiliar para  update wavy anchors for atom.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        for item in list(self._wavy_anchors):
            if item.data(WAVY_ANCHOR_ROLE) != atom_id:
                continue
            self._position_wavy_anchor(item)

    def _bond_angles_for_atom(self, atom_id: int) -> list[float]:
        """Método auxiliar para  bond angles for atom.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atom = self.model.get_atom(atom_id)
        angles: list[float] = []
        for bond in self.model.bonds.values():
            if bond.a1_id == atom_id:
                other = self.model.get_atom(bond.a2_id)
            elif bond.a2_id == atom_id:
                other = self.model.get_atom(bond.a1_id)
            else:
                continue
            dx = other.x - atom.x
            dy = other.y - atom.y
            angle = (math.degrees(math.atan2(dy, dx)) + 360.0) % 360.0
            angles.append(angle)
        return angles

    def _occupied_electron_slots(self, atom_id: int) -> set[int]:
        """Método auxiliar para  occupied electron slots.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        occupied: set[int] = set()
        for item in self._electron_dots:
            if item.data(ELECTRON_ANCHOR_ROLE) != atom_id:
                continue
            slot_idx = item.data(ELECTRON_SLOT_ROLE)
            if slot_idx is None:
                continue
            occupied.add(int(slot_idx))
        return occupied

    def _candidate_electron_slots(self, atom_id: int, mode: str) -> list[int]:
        """Método auxiliar para  candidate electron slots.

        Args:
            atom_id: Descripción del parámetro.
            mode: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        bond_angles = self._bond_angles_for_atom(atom_id)
        bond_count = len(bond_angles)
        if mode != "lone_pair":
            return list(range(len(ELECTRON_SLOT_ANGLES)))

        # Lone pair tool: 90-degree spacing when free, tetrahedral when bonded.
        if bond_count == 0:
            return [0, 2, 4, 6]

        ideal_angles: list[float] = []
        if bond_count == 1:
            b = bond_angles[0]
            ideal_angles = [b + 109.5, b - 109.5, b + 180.0]
        elif bond_count == 2:
            b1, b2 = bond_angles[0], bond_angles[1]
            ideal_angles = [
                b1 + 109.5,
                b1 - 109.5,
                b2 + 109.5,
                b2 - 109.5,
            ]
        else:
            ideal_angles = []

        candidates: list[int] = []
        if ideal_angles:
            for ang in ideal_angles:
                ang = (ang + 360.0) % 360.0
                best_idx = None
                best_dist = 1e9
                for idx, slot_angle in enumerate(ELECTRON_SLOT_ANGLES):
                    dist = angle_distance_deg(slot_angle, ang)
                    if dist < best_dist:
                        best_dist = dist
                        best_idx = idx
                if best_idx is not None and best_idx not in candidates:
                    candidates.append(best_idx)

        if not candidates:
            candidates = list(range(len(ELECTRON_SLOT_ANGLES)))

        # Remove slots blocked by bonds
        filtered: list[int] = []
        for idx in candidates:
            slot_angle = ELECTRON_SLOT_ANGLES[idx]
            blocked = any(
                angle_distance_deg(slot_angle, b) <= ELECTRON_SLOT_TOLERANCE_DEG
                for b in bond_angles
            )
            if not blocked:
                filtered.append(idx)
        return filtered if filtered else candidates

    def _select_electron_slot(
        self,
        atom_id: int,
        scene_pos: QPointF,
        candidate_slots: Optional[list[int]] = None,
    ) -> Optional[int]:
        """Método auxiliar para  select electron slot.

        Args:
            atom_id: Descripción del parámetro.
            scene_pos: Descripción del parámetro.
            candidate_slots: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atom = self.model.get_atom(atom_id)
        center = QPointF(atom.x, atom.y)
        direction = scene_pos - center
        if direction.manhattanLength() < 2.0:
            direction = self._label_open_direction(atom_id)
        if direction.manhattanLength() < 1e-3:
            direction = QPointF(1.0, 0.0)
        target_angle = (math.degrees(math.atan2(direction.y(), direction.x())) + 360.0) % 360.0

        bond_angles = self._bond_angles_for_atom(atom_id)
        occupied = self._occupied_electron_slots(atom_id)

        if candidate_slots is None:
            candidate_slots = list(range(len(ELECTRON_SLOT_ANGLES)))

        best_idx: Optional[int] = None
        best_dist = 1e9
        for idx in candidate_slots:
            slot_angle = ELECTRON_SLOT_ANGLES[idx]
            if idx in occupied:
                continue
            blocked = any(angle_distance_deg(slot_angle, b) <= ELECTRON_SLOT_TOLERANCE_DEG for b in bond_angles)
            if blocked:
                continue
            dist = angle_distance_deg(slot_angle, target_angle)
            if dist < best_dist:
                best_dist = dist
                best_idx = idx
        return best_idx

    def _electron_slot_position(self, atom_id: int, slot_idx: int, scale: float) -> tuple[QPointF, QPointF]:
        """Método auxiliar para  electron slot position.

        Args:
            atom_id: Descripción del parámetro.
            slot_idx: Descripción del parámetro.
            scale: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atom = self.model.get_atom(atom_id)
        angle = math.radians(ELECTRON_SLOT_ANGLES[slot_idx])
        nx = math.cos(angle)
        ny = math.sin(angle)
        px, py = -ny, nx
        radial_offset = self._electron_radial_offset(atom_id, scale)
        center = QPointF(atom.x, atom.y)
        pos = QPointF(center.x() + nx * radial_offset, center.y() + ny * radial_offset)
        tangent = QPointF(px, py)
        return pos, tangent

    def _position_electron_dot(self, item: TextAnnotationItem) -> None:
        """Método auxiliar para  position electron dot.

        Args:
            item: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if not item.data(ELECTRON_DOT_ROLE):
            return
        anchor_id = item.data(ELECTRON_ANCHOR_ROLE)
        slot_idx = item.data(ELECTRON_SLOT_ROLE)
        side = item.data(ELECTRON_SIDE_ROLE) or 0
        if anchor_id is None or slot_idx is None:
            return
        scale = float(item.data(ELECTRON_SCALE_ROLE) or 1.0)
        pos, tangent = self._electron_slot_position(int(anchor_id), int(slot_idx), scale)
        offset = self._electron_pair_spread(int(anchor_id), scale) * int(side)
        pos = QPointF(pos.x() + tangent.x() * offset, pos.y() + tangent.y() * offset)
        # Update font size to match atom label scale
        item_font = item.font()
        item_font.setPointSizeF(self._electron_dot_font_size(int(anchor_id), scale))
        item_font.setBold(False)
        item.setFont(item_font)
        rect = item.boundingRect()
        item.setPos(pos.x() - rect.width() / 2, pos.y() - rect.height() / 2)

    def _update_electron_dots_for_atom(self, atom_id: int) -> None:
        """Método auxiliar para  update electron dots for atom.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        for item in list(self._electron_dots):
            if item.data(ELECTRON_ANCHOR_ROLE) != atom_id:
                continue
            self._position_electron_dot(item)

    def _create_dot_item(
        self,
        pos: QPointF,
        atom_item: Optional[AtomItem],
        scale: float,
        anchor_atom_id: Optional[int] = None,
        slot_idx: Optional[int] = None,
        side: int = 0,
    ) -> None:
        """Método auxiliar para  create dot item.

        Args:
            pos: Descripción del parámetro.
            atom_item: Descripción del parámetro.
            scale: Descripción del parámetro.
            anchor_atom_id: Descripción del parámetro.
            slot_idx: Descripción del parámetro.
            side: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = TextAnnotationItem("•", 0.0, 0.0)
        self._apply_text_settings(item)
        item.setDefaultTextColor(QColor("#222222"))
        item.setTextInteractionFlags(Qt.TextInteractionFlag.NoTextInteraction)
        try:
            item.document().setDocumentMargin(0)
            item.document().setDefaultStyleSheet("body { background: transparent; }")
        except Exception:
            pass
        cursor = item.textCursor()
        cursor.clearSelection()
        item.setTextCursor(cursor)
        fmt = QTextCharFormat()
        fmt.setBackground(QBrush(Qt.BrushStyle.NoBrush))
        cursor.select(cursor.SelectionType.Document)
        cursor.mergeCharFormat(fmt)

        if scale and abs(scale - 1.0) > 1e-3:
            font = item.font()
            font.setPointSizeF(font.pointSizeF() * scale)
            item.setFont(font)
        if anchor_atom_id is not None:
            font = item.font()
            font.setPointSizeF(self._electron_dot_font_size(anchor_atom_id, scale))
            font.setBold(False)
            item.setFont(font)

        rect = item.boundingRect()
        item.setPos(pos.x() - rect.width() / 2, pos.y() - rect.height() / 2)
        item.setData(ELECTRON_DOT_ROLE, True)
        if anchor_atom_id is not None:
            item.setData(ELECTRON_ANCHOR_ROLE, anchor_atom_id)
        if slot_idx is not None:
            item.setData(ELECTRON_SLOT_ROLE, slot_idx)
            item.setData(ELECTRON_SIDE_ROLE, int(side))
        item.setData(ELECTRON_SCALE_ROLE, float(scale))
        self._electron_dots.add(item)

        self.undo_stack.push(AddTextItemCommand(self, item))

    def _insert_electron_dots(
        self,
        scene_pos: QPointF,
        atom_id: Optional[int],
        count: int,
        scale: float,
        spread: Optional[float] = None,
        mode: str = "default",
    ) -> None:
        """Método auxiliar para  insert electron dots.

        Args:
            scene_pos: Descripción del parámetro.
            atom_id: Descripción del parámetro.
            count: Descripción del parámetro.
            scale: Descripción del parámetro.
            spread: Descripción del parámetro.
            mode: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if atom_id is None:
            return
        atom = self.model.get_atom(atom_id)
        if atom.element not in HETERO_ELECTRON_ATOMS:
            return

        slot_candidates = self._candidate_electron_slots(atom_id, mode)
        slot_idx = self._select_electron_slot(atom_id, scene_pos, slot_candidates)
        if slot_idx is None:
            return

        pos, tangent = self._electron_slot_position(atom_id, slot_idx, scale)
        if count <= 1:
            self._create_dot_item(
                pos,
                self.atom_items.get(atom_id),
                scale,
                anchor_atom_id=atom_id,
                slot_idx=slot_idx,
                side=0,
            )
            return

        base_spread = self._electron_pair_spread(atom_id, scale)
        if mode == "lone_pair":
            base_spread *= 1.25
        use_spread = base_spread if spread is None else spread
        dx = tangent.x() * use_spread
        dy = tangent.y() * use_spread
        p1 = QPointF(pos.x() + dx, pos.y() + dy)
        p2 = QPointF(pos.x() - dx, pos.y() - dy)
        self._create_dot_item(
            p1,
            self.atom_items.get(atom_id),
            scale,
            anchor_atom_id=atom_id,
            slot_idx=slot_idx,
            side=1,
        )
        self._create_dot_item(
            p2,
            self.atom_items.get(atom_id),
            scale,
            anchor_atom_id=atom_id,
            slot_idx=slot_idx,
            side=-1,
        )

    def _insert_wavy_anchor(self, scene_pos: QPointF, atom_id: int) -> None:
        """Método auxiliar para  insert wavy anchor.

        Args:
            scene_pos: Descripción del parámetro.
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if atom_id not in self.model.atoms:
            return
        bond_id, base_angle = self._wavy_anchor_bond_angle(atom_id, scene_pos)
        length = self._wavy_anchor_length()

        atom = self.model.get_atom(atom_id)
        perp_angle = math.radians(base_angle + 90.0)
        nx = math.cos(perp_angle)
        ny = math.sin(perp_angle)
        center = QPointF(atom.x, atom.y)
        half = length * 0.5
        start = QPointF(center.x() - nx * half, center.y() - ny * half)
        end = QPointF(center.x() + nx * half, center.y() + ny * half)

        item = WavyAnchorItem(start, end, style=self.drawing_style)
        item.setData(WAVY_ANCHOR_ROLE, atom_id)
        item.setData(WAVY_ANCHOR_ANGLE_ROLE, base_angle)
        item.setData(WAVY_ANCHOR_LENGTH_ROLE, length)
        if bond_id is not None:
            item.setData(WAVY_ANCHOR_BOND_ROLE, int(bond_id))
        self._wavy_anchors.add(item)
        self.undo_stack.push(AddWavyAnchorCommand(self, item))

    def _insert_symbol_item(
        self,
        text: str,
        scene_pos: QPointF,
        atom_id: Optional[int],
        scale: float,
        anchor_to_atom: bool,
        rotate: bool = False,
    ) -> None:
        """Método auxiliar para  insert symbol item.

        Args:
            text: Descripción del parámetro.
            scene_pos: Descripción del parámetro.
            atom_id: Descripción del parámetro.
            scale: Descripción del parámetro.
            anchor_to_atom: Descripción del parámetro.
            rotate: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        pos = scene_pos
        angle = 0.0
        if anchor_to_atom and atom_id is not None:
            pos, angle = self._symbol_anchor_data(atom_id, scene_pos)
        item = TextAnnotationItem(text, 0.0, 0.0)
        self._apply_text_settings(item)
        item.setDefaultTextColor(QColor("#222222"))
        item.setTextInteractionFlags(Qt.TextInteractionFlag.NoTextInteraction)
        try:
            item.document().setDocumentMargin(0)
            item.document().setDefaultStyleSheet("body { background: transparent; }")
        except Exception:
            pass
        cursor = item.textCursor()
        cursor.clearSelection()
        item.setTextCursor(cursor)
        fmt = QTextCharFormat()
        fmt.setBackground(QBrush(Qt.BrushStyle.NoBrush))
        cursor.select(cursor.SelectionType.Document)
        cursor.mergeCharFormat(fmt)
        if scale and abs(scale - 1.0) > 1e-3:
            font = item.font()
            font.setPointSizeF(font.pointSizeF() * scale)
            item.setFont(font)
        rect = item.boundingRect()
        if anchor_to_atom and atom_id is not None and atom_id in self.atom_items:
            atom_item = self.atom_items[atom_id]
            local = atom_item.mapFromScene(QPointF(pos.x(), pos.y()))
            item.setParentItem(atom_item)
            item.setPos(local.x() - rect.width() / 2, local.y() - rect.height() / 2)
        else:
            item.setPos(pos.x() - rect.width() / 2, pos.y() - rect.height() / 2)
        if rotate:
            item.setRotation(angle)
        self.undo_stack.push(AddTextItemCommand(self, item))

    @staticmethod
    def _normalize_atom_label(text: str) -> str | None:
        """Método auxiliar para  normalize atom label.

        Args:
            text: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        cleaned = text.strip()
        if not cleaned:
            return None
        allowed = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+-()_^")
        if any(ch not in allowed for ch in cleaned):
            return None
        alias = FUNCTIONAL_GROUP_ALIASES.get(cleaned.lower())
        if alias:
            return alias
        if "_" in cleaned or "^" in cleaned:
            return cleaned
        if cleaned.isalpha() and len(cleaned) <= 2:
            if len(cleaned) == 1:
                normalized = cleaned.upper()
            else:
                normalized = cleaned[0].upper() + cleaned[1:].lower()
            if normalized in ELEMENT_SYMBOLS:
                return normalized
        return cleaned

    def resolve_atom_label_spec(self, label: str) -> dict:
        """Resuelve una etiqueta de UI a la especificación química persistente."""
        normalized = self._normalize_atom_label(label) or label.strip()
        if not normalized:
            return {
                "element": "",
                "group_h_cap": None,
                "explicit_h": None,
                "no_implicit": False,
            }
        simple_group = SIMPLE_HYDROGEN_GROUP_LABELS.get(normalized)
        if simple_group is None:
            simple_group = SIMPLE_HYDROGEN_GROUP_LABELS.get(normalized.upper())
        if simple_group is not None:
            element, group_h_cap = simple_group
            return {
                "element": element,
                "group_h_cap": int(group_h_cap),
                "explicit_h": None,
                "no_implicit": True,
            }
        return {
            "element": normalized,
            "group_h_cap": None,
            "explicit_h": None,
            "no_implicit": False,
        }

    def _editable_atom_label(self, atom) -> str:
        """Devuelve la etiqueta canónica para edición en diálogos."""
        cap = getattr(atom, "group_h_cap", None)
        if cap is not None:
            canonical = CANONICAL_SIMPLE_HYDROGEN_GROUP_LABELS.get(
                (atom.element, int(cap))
            )
            if canonical:
                return canonical
        return atom.element

    def _atom_label_items(self) -> list[str]:
        """Método auxiliar para  atom label items.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        elements = ["C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "H"]
        return elements + FUNCTIONAL_GROUP_LABELS

    def _label_anchor_candidates(self, label: str) -> list[str]:
        """Método auxiliar para  label anchor candidates.

        Args:
            label: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        cleaned = label.strip()
        if not cleaned:
            return []
        tokens = self._group_label_tokens(cleaned)
        if tokens:
            seen: set[str] = set()
            ordered: list[str] = []
            for symbol, _token in tokens:
                if symbol in seen:
                    continue
                seen.add(symbol)
                ordered.append(symbol)
            return ordered
        ordered: list[str] = []
        seen: set[str] = set()
        i = 0
        while i < len(cleaned):
            ch = cleaned[i]
            if not ch.isalpha():
                i += 1
                continue
            if ch.isupper():
                symbol = None
                if i + 1 < len(cleaned) and cleaned[i + 1].islower():
                    candidate = cleaned[i : i + 2]
                    if candidate in ELEMENT_SYMBOLS:
                        symbol = candidate
                        i += 2
                if symbol is None and ch in ELEMENT_SYMBOLS:
                    symbol = ch
                    i += 1
                if symbol is None:
                    i += 1
                    continue
                if symbol not in seen:
                    seen.add(symbol)
                    ordered.append(symbol)
                continue
            i += 1
        return ordered

    def _prompt_anchor_for_atom(self, atom_id: int) -> None:
        """Método auxiliar para  prompt anchor for atom.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atom = self.model.get_atom(atom_id)
        if atom is None:
            return
        label = atom.element
        candidates = self._label_anchor_candidates(label)
        if not candidates:
            return
        items = ["Auto"] + candidates
        current = self._group_anchor_overrides.get(atom_id) or "Auto"
        index = items.index(current) if current in items else 0
        value, ok = QInputDialog.getItem(
            self,
            "Átomo de unión",
            "Átomo de unión:",
            items,
            index,
            True,
        )
        if not ok:
            return
        value = value.strip()
        anchor = None
        if value and value != "Auto":
            anchor = self._normalize_atom_label(value)
        self.set_anchor_override(atom_id, anchor)
        self._refresh_atom_label(atom_id)

    def _cycle_anchor_override(self, atom_id: int) -> bool:
        """Método auxiliar para  cycle anchor override.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atom = self.model.get_atom(atom_id)
        if atom is None:
            return False
        if atom.element in ELEMENT_SYMBOLS:
            return False
        candidates = self._label_anchor_candidates(atom.element)
        if not candidates:
            return False
        current = self._group_anchor_overrides.get(atom_id)
        if current in candidates:
            idx = (candidates.index(current) + 1) % len(candidates)
        else:
            idx = 0
        self.set_anchor_override(atom_id, candidates[idx])
        self._refresh_atom_label(atom_id)
        return True

    def _prompt_atom_label(
        self, current: str, atom_id: Optional[int] = None, initial: str | None = None
    ) -> tuple[Optional[str], Optional[str]]:
        """Método auxiliar para  prompt atom label.

        Args:
            current: Descripción del parámetro.
            atom_id: Descripción del parámetro.
            initial: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        items = self._atom_label_items()
        seed = (initial or "").strip() or current
        if seed:
            if seed not in items:
                items = [seed] + items
            current_index = items.index(seed)
        else:
            current_index = 0
        anchor = self._group_anchor_overrides.get(atom_id) if atom_id is not None else None
        dialog = AtomLabelDialog(
            seed,
            anchor,
            items,
            ELEMENT_SYMBOLS,
            self,
        )
        if not dialog.exec():
            return None, None
        value, anchor_value = dialog.value()
        normalized = self._normalize_atom_label(value)
        if not normalized:
            return None, None
        anchor_normalized = None
        if anchor_value:
            anchor_normalized = self._normalize_atom_label(anchor_value)
        # If the user chose an anchor for a group label, reorder the label now
        # so the visible text reflects the chosen attachment atom.
        if (
            anchor_normalized
            and atom_id is not None
            and normalized not in ELEMENT_SYMBOLS
        ):
            normalized, _ = self._reflow_group_label(
                normalized, atom_id, anchor_normalized
            )
        return normalized, anchor_normalized

    def add_text_item(self, item: TextAnnotationItem) -> None:
        """Añade text item.

        Args:
            item: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self.ensure_graphics_item_opacity(item)
        if item.scene() is not self.scene:
            self.scene.addItem(item)

    def remove_text_item(self, item: TextAnnotationItem) -> None:
        """Elimina text item.

        Args:
            item: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if item.scene() is self.scene:
            self.scene.removeItem(item)
        if item in self._electron_dots:
            self._electron_dots.discard(item)

    def readd_text_item(self, item: TextAnnotationItem) -> None:
        """Método auxiliar para readd text item.

        Args:
            item: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self.ensure_graphics_item_opacity(item)
        if item.scene() is not self.scene:
            self.scene.addItem(item)
        if item.data(ELECTRON_DOT_ROLE):
            self._electron_dots.add(item)
            self._position_electron_dot(item)

    def _label_font(self) -> QFont:
        """Método auxiliar para  label font.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        size = self.state.label_font_size
        if size <= 0:
            size = 10.0
        font = QFont(self.state.label_font_family)
        font.setPointSizeF(size)
        font.setBold(self.state.label_font_bold)
        font.setItalic(self.state.label_font_italic)
        font.setUnderline(self.state.label_font_underline)
        return font

    def _atom_label_scale_factor(self, atom_id: int) -> float:
        """Devuelve la escala efectiva de etiqueta para un átomo."""
        atom = self.model.atoms.get(atom_id)
        value = getattr(atom, "label_scale", None) if atom is not None else None
        try:
            scale = float(value) if value is not None else 1.0
        except Exception:
            scale = 1.0
        return max(0.2, scale)

    def _label_font_for_atom(self, atom_id: int) -> QFont:
        """Construye la fuente efectiva para la etiqueta de un átomo."""
        font = self._label_font()
        size = font.pointSizeF()
        if size <= 0.0:
            size = float(font.pointSize()) if font.pointSize() > 0 else 10.0
        font.setPointSizeF(size * self._atom_label_scale_factor(atom_id))
        return font

    def _effective_label_font_size(self, atom_id: int) -> float:
        """Devuelve el tamaño efectivo actual de la etiqueta de un átomo."""
        font = self._label_font_for_atom(atom_id)
        size = font.pointSizeF()
        if size <= 0.0:
            size = float(font.pointSize()) if font.pointSize() > 0 else 10.0
        return max(1.0, float(size))

    def _selected_structure_label_atom_ids(self) -> list[int]:
        """Devuelve átomos estructurales afectados por cambios de tamaño de etiquetas."""
        return sorted(
            atom_id
            for atom_id in self._selected_atom_ids_for_transform()
            if atom_id in self.model.atoms
        )

    def _label_scale_from_target_size(self, size_pt: float) -> Optional[float]:
        """Convierte un tamaño objetivo en escala local respecto a la fuente global."""
        base_size = float(self.state.label_font_size)
        if base_size <= 0.0:
            base_size = 10.0
        target_size = max(1.0, float(size_pt))
        scale = target_size / base_size
        if abs(scale - 1.0) < 0.02:
            return None
        return max(0.2, scale)

    def label_font(self) -> QFont:
        """Método auxiliar para label font.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        return self._label_font()

    def current_label_size_value(self) -> float:
        """Devuelve el tamaño actual a mostrar en el diálogo de etiquetas."""
        selected_atom_ids = self._selected_structure_label_atom_ids()
        if selected_atom_ids:
            return self._effective_label_font_size(selected_atom_ids[0])
        size = float(self.state.label_font_size)
        return size if size > 0.0 else 10.0

    def apply_label_font(self, font: QFont) -> None:
        """Aplica label font.

        Args:
            font: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        size = font.pointSizeF()
        if size <= 0:
            size = font.pointSize()
        if size <= 0:
            size = 10.0
        self.state.label_font_family = font.family()
        self.state.label_font_size = float(size)
        self.state.label_font_bold = font.bold()
        self.state.label_font_italic = font.italic()
        self.state.label_font_underline = font.underline()
        self.refresh_label_fonts()

    def apply_label_font_size(self, size_pt: float) -> bool:
        """Aplica tamaño de etiqueta a selección estructural o globalmente."""
        target_size = max(6.0, float(size_pt))
        selected_atom_ids = self._selected_structure_label_atom_ids()
        if not selected_atom_ids:
            font = self.label_font()
            font.setPointSizeF(target_size)
            self.apply_label_font(font)
            return True

        changes: list[tuple[int, Optional[float], Optional[float]]] = []
        for atom_id in selected_atom_ids:
            atom = self.model.get_atom(atom_id)
            old_scale = getattr(atom, "label_scale", None)
            new_scale = self._label_scale_from_target_size(target_size)
            if self._optional_float_equal(old_scale, new_scale, tol=0.01):
                continue
            changes.append((atom_id, old_scale, new_scale))

        if not changes:
            return False

        if len(changes) > 1:
            self.undo_stack.beginMacro("Set selected label size")
        for atom_id, old_scale, new_scale in changes:
            self.undo_stack.push(
                ChangeAtomLabelScaleCommand(
                    self.model,
                    self,
                    atom_id,
                    new_scale,
                    old_scale=old_scale,
                )
            )
        if len(changes) > 1:
            self.undo_stack.endMacro()
        return True

    def adjust_label_font_size(self, delta: float) -> bool:
        """Aumenta o reduce tamaño de etiquetas en selección o globalmente."""
        delta_value = float(delta)
        selected_atom_ids = self._selected_structure_label_atom_ids()
        if not selected_atom_ids:
            font = self.label_font()
            size = font.pointSizeF()
            if size <= 0.0:
                size = float(font.pointSize()) if font.pointSize() > 0 else 10.0
            font.setPointSizeF(max(6.0, float(size) + delta_value))
            self.apply_label_font(font)
            return True

        changes: list[tuple[int, Optional[float], Optional[float]]] = []
        for atom_id in selected_atom_ids:
            atom = self.model.get_atom(atom_id)
            old_scale = getattr(atom, "label_scale", None)
            current_size = self._effective_label_font_size(atom_id)
            new_scale = self._label_scale_from_target_size(max(6.0, current_size + delta_value))
            if self._optional_float_equal(old_scale, new_scale, tol=0.01):
                continue
            changes.append((atom_id, old_scale, new_scale))

        if not changes:
            return False

        if len(changes) > 1:
            self.undo_stack.beginMacro("Adjust selected label size")
        for atom_id, old_scale, new_scale in changes:
            self.undo_stack.push(
                ChangeAtomLabelScaleCommand(
                    self.model,
                    self,
                    atom_id,
                    new_scale,
                    old_scale=old_scale,
                )
            )
        if len(changes) > 1:
            self.undo_stack.endMacro()
        return True

    def refresh_label_fonts(self, atom_ids: Optional[Iterable[int]] = None) -> None:
        """Método auxiliar para refresh label fonts.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if atom_ids is None:
            atom_ids = list(self.atom_items.keys())
        atom_ids = list(atom_ids)
        for atom_id in atom_ids:
            item = self.atom_items.get(atom_id)
            if item is None:
                continue
            item.set_label_font(self._label_font_for_atom(atom_id))
        self.refresh_atom_labels(atom_ids)
        self.recompute_numbering()

    def set_use_element_colors(self, use_element_colors: bool) -> None:
        """Actualiza el estado de use element colors.

        Args:
            use_element_colors: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self.state.use_element_colors = use_element_colors
        for item in self.atom_items.values():
            item.set_use_element_colors(use_element_colors)
        self.refresh_atom_labels()

    def refresh_atom_labels(self, atom_ids: Optional[Iterable[int]] = None) -> None:
        """Método auxiliar para refresh atom labels.

        Args:
            atom_ids: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if atom_ids is None:
            atom_ids = list(self.atom_items.keys())
        for atom_id in atom_ids:
            self._refresh_atom_label(atom_id)

    def _refresh_atom_label(self, atom_id: int) -> None:
        """Método auxiliar para  refresh atom label.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = self.atom_items.get(atom_id)
        atom = self.model.atoms.get(atom_id)
        if item is None or atom is None:
            return
        label, anchor, offset = self._build_atom_label(atom)
        item.set_display_label(label, anchor)
        item.set_label_offset(offset)
        self._update_bond_label_shrinks({atom_id})

    def _assigned_hydrogen_count(self, atom_id: int) -> int:
        """Devuelve H asignados al átomo excluyendo H como nodos separados."""
        counter = getattr(self.model, "assigned_hydrogen_count", None)
        if callable(counter):
            try:
                return int(max(0, counter(atom_id)))
            except Exception:
                pass
        atom = self.model.atoms.get(atom_id)
        return int(getattr(atom, "explicit_h", 0) or 0) if atom is not None else 0

    def _resolved_display_element(self, atom) -> str:
        """Resuelve alias legacy simples a su elemento visible real."""
        display_element = atom.element
        if display_element not in ELEMENT_SYMBOLS:
            legacy_spec = self.resolve_atom_label_spec(display_element)
            legacy_element = str(legacy_spec.get("element") or "")
            if legacy_element in ELEMENT_SYMBOLS and legacy_spec.get("group_h_cap") is not None:
                display_element = legacy_element
        return display_element

    def _inline_label_hydrogen_count(self, atom_id: int, display_element: str) -> int:
        """Cuenta los H que se dibujan inline junto al símbolo del heteroátomo."""
        atom = self.model.atoms.get(atom_id)
        if (
            atom is None
            or display_element not in ELEMENT_SYMBOLS
            or bool(getattr(atom, "is_coordination_center", False))
        ):
            return 0
        assigned_h = self._assigned_hydrogen_count(atom_id)
        implicit_h = self._implicit_hydrogen_count(atom_id, display_element)
        if (
            implicit_h > 0
            and display_element in {"O", "S"}
            and self._atom_degree(atom_id) >= 2
        ):
            implicit_h = 0
        inline_implicit_h = 0 if self.state.show_implicit_hydrogens else implicit_h
        return int(max(0, assigned_h + inline_implicit_h))

    def _build_atom_label(self, atom) -> tuple[str, Optional[str], QPointF]:
        """Método auxiliar para  build atom label.

        Args:
            atom: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        label = atom.element
        display_element = self._resolved_display_element(atom)
        anchor: Optional[str] = None
        anchor_override = self._group_anchor_overrides.get(atom.id)
        is_coordination_center = bool(getattr(atom, "is_coordination_center", False))
        if display_element in ELEMENT_SYMBOLS:
            if (
                not is_coordination_center
                and display_element == "C"
                and not (
                self.state.show_implicit_carbons or atom.is_explicit
                )
            ):
                return display_element, None, QPointF(0.0, 0.0)
            total_h = self._inline_label_hydrogen_count(atom.id, display_element)
            if (
                total_h > 0
                and display_element != "H"
                and not is_coordination_center
            ):
                h_text = "H" if total_h == 1 else f"H{total_h}"
                if self._prefer_prefix_h(atom.id, display_element, total_h):
                    label = f"{h_text}{display_element}"
                    anchor = display_element
                else:
                    label = f"{display_element}{h_text}"
            else:
                label = display_element
        else:
            label, anchor = self._reflow_group_label(label, atom.id, anchor_override)
            if anchor is None and anchor_override and anchor_override in label:
                anchor = anchor_override
        offset = self._label_offset(atom.id)
        return label, anchor, offset

    def _atom_degree(self, atom_id: int) -> int:
        """Método auxiliar para  atom degree.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        degree_getter = getattr(self.model, "atom_degree", None)
        if callable(degree_getter):
            try:
                return int(degree_getter(atom_id))
            except Exception:
                pass
        return sum(
            1
            for bond in self.model.bonds.values()
            if bond_is_structural(bond)
            and (bond.a1_id == atom_id or bond.a2_id == atom_id)
        )

    def _is_disposable_orphan_atom(self, atom_id: int) -> bool:
        """Indica si el átomo es un placeholder invisible que no debe interferir con edición."""
        checker = getattr(self.model, "is_disposable_orphan_atom", None)
        if callable(checker):
            return bool(checker(atom_id))
        atom = self.model.atoms.get(atom_id)
        return bool(
            atom is not None
            and atom.element == "C"
            and not atom.is_explicit
            and self._atom_degree(atom_id) <= 0
        )

    def _interactive_atom_candidates(self) -> list[tuple[int, float, float]]:
        """Devuelve átomos válidos para hover, snap y colocación de enlaces."""
        return [
            (atom.id, atom.x, atom.y)
            for atom in self.model.atoms.values()
            if not self._is_disposable_orphan_atom(atom.id)
        ]

    def _reflow_group_label(
        self, label: str, atom_id: int, anchor_override: Optional[str] = None
    ) -> tuple[str, Optional[str]]:
        """Método auxiliar para  reflow group label.

        Args:
            label: Descripción del parámetro.
            atom_id: Descripción del parámetro.
            anchor_override: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        cleaned = label.strip()
        if not cleaned:
            return label, None

        charge = ""
        if cleaned and cleaned[-1] in "+-":
            charge = cleaned[-1]
            cleaned = cleaned[:-1]

        tokens = self._group_label_tokens(cleaned)
        # If the user explicitly chose an anchor atom, force it to the front.
        # This keeps labels like "TBSO" consistent with an O-anchored bond.
        if anchor_override:
            anchor_first = True
        else:
            direction = self._label_open_direction(atom_id)
            anchor_first = direction.x() >= -0.2

        if anchor_override and anchor_override in cleaned:
            if not tokens or len(tokens) < 2:
                # Use bond direction to decide if anchor goes first or last
                cleaned = self._move_anchor_in_label(cleaned, anchor_override, anchor_first)
                return f"{cleaned}{charge}", anchor_override
            if not any(symbol == anchor_override for symbol, _token in tokens):
                # Anchor letter is inside an abbreviation token (e.g., TBSO -> B/S).
                cleaned = self._move_anchor_in_label(cleaned, anchor_override, anchor_first)
                return f"{cleaned}{charge}", anchor_override
        if not tokens or len(tokens) < 2:
            return label, None

        if anchor_override and any(symbol == anchor_override for symbol, _token in tokens):
            if anchor_first:
                while tokens[0][0] != anchor_override:
                    tokens = tokens[1:] + tokens[:1]
            else:
                while tokens[-1][0] != anchor_override:
                    tokens = tokens[1:] + tokens[:1]
            anchor_symbol = anchor_override
        else:
            anchor_symbol = tokens[0][0]
            if not anchor_first:
                # Move first token to end
                tokens = tokens[1:] + tokens[:1]
                anchor_symbol = tokens[-1][0]
        cleaned = "".join(token for _symbol, token in tokens)
        return f"{cleaned}{charge}", anchor_symbol

    def _move_anchor_in_label(self, label: str, anchor: str, anchor_first: bool) -> str:
        """Método auxiliar para  move anchor in label.

        Args:
            label: Descripción del parámetro.
            anchor: Descripción del parámetro.
            anchor_first: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if not label or not anchor:
            return label
        if anchor_first:
            idx = label.find(anchor)
        else:
            idx = label.rfind(anchor)
        if idx < 0:
            return label
        before = label[:idx]
        after = label[idx + len(anchor):]
        core = before + after
        return f"{anchor}{core}" if anchor_first else f"{core}{anchor}"

    def get_anchor_override(self, atom_id: int) -> Optional[str]:
        """Método auxiliar para get anchor override.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        return self._group_anchor_overrides.get(atom_id)

    def set_anchor_override(self, atom_id: int, anchor: Optional[str]) -> None:
        """Actualiza el estado de anchor override.

        Args:
            atom_id: Descripción del parámetro.
            anchor: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if anchor:
            self._group_anchor_overrides[atom_id] = anchor
        else:
            self._group_anchor_overrides.pop(atom_id, None)

    def _group_label_tokens(self, label: str) -> Optional[list[tuple[str, str]]]:
        """Método auxiliar para  group label tokens.

        Args:
            label: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        tokens: list[tuple[str, str]] = []
        i = 0
        sorted_abbrs = sorted(list(ABBREVIATION_LABELS), key=len, reverse=True)
        
        while i < len(label):
            ch = label[i]
            symbol = None
            
            # 1. Try abbreviations first
            for abbr in sorted_abbrs:
                if label.startswith(abbr, i):
                    symbol = abbr
                    i += len(abbr)
                    break
            
            # 2. Try element symbols if no abbreviation matched
            if symbol is None:
                if not ch.isupper():
                    return None
                if i + 1 < len(label) and label[i + 1].islower():
                    candidate = label[i : i + 2]
                    if candidate in ELEMENT_SYMBOLS:
                        symbol = candidate
                        i += 2
                if symbol is None and ch in ELEMENT_SYMBOLS:
                    symbol = ch
                    i += 1
            
            if symbol is None:
                return None
                
            token = symbol
            # Collect numbers and script markers
            while i < len(label):
                if label[i].isdigit():
                    token += label[i]
                    i += 1
                    continue
                if label[i] in "_^":
                    marker = label[i]
                    token += marker
                    i += 1
                    while i < len(label) and label[i].isalnum():
                        token += label[i]
                        i += 1
                    continue
                break
            tokens.append((symbol, token))
            
        return tokens if tokens else None

    def _implicit_hydrogen_count(self, atom_id: int, element: str) -> int:
        """Método auxiliar para  implicit hydrogen count.

        Args:
            atom_id: Descripción del parámetro.
            element: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if element not in IMPLICIT_H_ELEMENTS:
            return 0
        try:
            return int(max(0, self.model.implicit_h_count(atom_id)))
        except Exception:
            return 0

    def _label_open_direction(self, atom_id: int) -> QPointF:
        """Método auxiliar para  label open direction.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        def collect_vectors(skip_hydrogen: bool) -> list[tuple[float, float]]:
            """Método auxiliar para collect vectors.

            Args:
                skip_hydrogen: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
            atom = self.model.get_atom(atom_id)
            vectors: list[tuple[float, float]] = []
            for bond in self.model.bonds.values():
                if bond.a1_id == atom_id:
                    other = self.model.get_atom(bond.a2_id)
                elif bond.a2_id == atom_id:
                    other = self.model.get_atom(bond.a1_id)
                else:
                    continue
                if skip_hydrogen and other.element == "H":
                    continue
                dx = other.x - atom.x
                dy = other.y - atom.y
                length = math.hypot(dx, dy)
                if length <= 1e-6:
                    continue
                vectors.append((dx / length, dy / length))
            return vectors

        vectors = collect_vectors(skip_hydrogen=True)
        if not vectors:
            vectors = collect_vectors(skip_hydrogen=False)
        if not vectors:
            return QPointF(0.0, 0.0)
        sum_x = sum(v[0] for v in vectors)
        sum_y = sum(v[1] for v in vectors)
        if abs(sum_x) + abs(sum_y) < 1e-3:
            if len(vectors) == 1:
                sum_x, sum_y = -vectors[0][0], -vectors[0][1]
            else:
                return QPointF(0.0, 0.0)
        else:
            sum_x, sum_y = -sum_x, -sum_y
        length = math.hypot(sum_x, sum_y)
        if length <= 1e-6:
            return QPointF(0.0, 0.0)
        return QPointF(sum_x / length, sum_y / length)

    def _label_offset(self, atom_id: int) -> QPointF:
        """Método auxiliar para  label offset.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atom = self.model.get_atom(atom_id)
        if atom is None:
            return QPointF(0.0, 0.0)
        base_offset = self._base_label_offset(atom_id)
        extra_offset = self._custom_bond_length_label_offset(atom_id, base_offset=base_offset)
        return base_offset + extra_offset

    def _base_label_offset(self, atom_id: int) -> QPointF:
        """Calcula el offset base de etiqueta sin ajustes por longitud personalizada."""
        atom = self.model.get_atom(atom_id)
        if atom is None:
            return QPointF(0.0, 0.0)
        display_element = self._resolved_display_element(atom)
        total_h = self._inline_label_hydrogen_count(atom_id, display_element)
        if (
            display_element in ELEMENT_SYMBOLS
            and display_element not in {"C", "H"}
            and total_h > 0
        ):
            direction = self._label_open_direction(atom_id)
            if direction == QPointF(0.0, 0.0):
                return QPointF(0.0, 0.0)
            length = math.hypot(direction.x(), direction.y())
            if length <= 1e-6:
                return QPointF(0.0, 0.0)
            direction = QPointF(direction.x() / length, direction.y() / length)
            if self._atom_degree(atom_id) >= 2 and abs(direction.y()) < 0.35:
                biased_y = 0.35 if direction.y() >= 0.0 else -0.35
                biased_length = math.hypot(direction.x(), biased_y)
                if biased_length > 1e-6:
                    direction = QPointF(direction.x() / biased_length, biased_y / biased_length)
            size = self._label_font_for_atom(atom_id).pointSizeF()
            if size <= 0:
                size = self._label_font().pointSizeF()
            if size <= 0:
                size = 10.0
            offset = max(3.0, size * 0.38)
            return QPointF(direction.x() * offset, direction.y() * offset)
        if any(ch.isalpha() for ch in atom.element):
            return QPointF(0.0, 0.0)
        direction = self._label_open_direction(atom_id)
        if direction == QPointF(0.0, 0.0):
            return QPointF(0.0, 0.0)
        size = self._label_font().pointSizeF()
        if size <= 0:
            size = 10.0
        offset = max(LABEL_OFFSET_MIN_PX, size * LABEL_OFFSET_SCALE)
        return QPointF(direction.x() * offset, direction.y() * offset)

    def _label_layout_center(self, atom_id: int, offset: QPointF) -> Optional[QPointF]:
        """Estima el centro en escena de la etiqueta usando un offset propuesto."""
        item = self.atom_items.get(atom_id)
        if item is None or not item.label.isVisible():
            return None
        rect = item.label.mapRectToParent(item.label.boundingRect())
        current_offset = QPointF(getattr(item, "_label_offset", QPointF(0.0, 0.0)))
        center = rect.center()
        center += QPointF(offset.x() - current_offset.x(), offset.y() - current_offset.y())
        item_pos = item.pos()
        return QPointF(item_pos.x() + center.x(), item_pos.y() + center.y())

    def _label_axis_extent(self, atom_id: int, ux: float, uy: float) -> float:
        """Semiextensión proyectada de la etiqueta/carga sobre el eje del enlace."""
        item = self.atom_items.get(atom_id)
        if item is None or not item.label.isVisible():
            return 0.0
        rect = item.label.mapRectToParent(item.label.boundingRect())
        if item.charge_label.isVisible():
            rect = rect.united(item.charge_label.mapRectToParent(item.charge_label.boundingRect()))
        if rect.isNull():
            return 0.0
        return 0.5 * (abs(float(ux)) * rect.width() + abs(float(uy)) * rect.height())

    def _custom_bond_length_label_offset(
        self,
        atom_id: int,
        *,
        base_offset: Optional[QPointF] = None,
    ) -> QPointF:
        """Desplaza una etiqueta hacia el extremo dibujado sin invadir la opuesta."""
        atom = self.model.atoms.get(atom_id)
        if atom is None or not self._atom_label_visible(atom):
            return QPointF(0.0, 0.0)
        if base_offset is None:
            base_offset = self._base_label_offset(atom_id)

        offset_x = 0.0
        offset_y = 0.0
        total_weight = 0.0
        for bond in self.model.bonds.values():
            if bond.a1_id == atom_id:
                other_atom = self.model.atoms.get(bond.a2_id)
            elif bond.a2_id == atom_id:
                other_atom = self.model.atoms.get(bond.a1_id)
            else:
                continue
            if other_atom is None:
                continue
            custom_length = getattr(bond, "length_px", None)
            if custom_length is None:
                continue
            try:
                target_length = float(custom_length)
            except Exception:
                continue
            if not math.isfinite(target_length) or target_length <= 1e-6:
                continue

            dx = other_atom.x - atom.x
            dy = other_atom.y - atom.y
            actual_length = math.hypot(dx, dy)
            if actual_length <= 1e-6:
                continue

            desired_shift = (actual_length - max(1.0, target_length)) * 0.5
            if abs(desired_shift) <= 1e-4:
                continue

            applied_shift = desired_shift
            ux = dx / actual_length
            uy = dy / actual_length
            if desired_shift > 0.0 and self._atom_label_visible(other_atom):
                other_base_offset = self._base_label_offset(other_atom.id)
                center = self._label_layout_center(atom_id, base_offset)
                other_center = self._label_layout_center(other_atom.id, other_base_offset)
                if center is not None and other_center is not None:
                    available = (
                        (other_center.x() - center.x()) * ux
                        + (other_center.y() - center.y()) * uy
                    )
                    required_gap = (
                        self._label_axis_extent(atom_id, ux, uy)
                        + self._label_axis_extent(other_atom.id, ux, uy)
                        + max(2.0, float(self.drawing_style.stroke_px))
                    )
                    max_inward_shift = max(0.0, (available - required_gap) * 0.5)
                    applied_shift = min(desired_shift, max_inward_shift)

            weight = abs(applied_shift)
            offset_x += ux * applied_shift * weight
            offset_y += uy * applied_shift * weight
            total_weight += weight

        if total_weight <= 1e-6:
            return QPointF(0.0, 0.0)
        return QPointF(offset_x / total_weight, offset_y / total_weight)

    def _neighbor_angles_deg(self, atom_id: int) -> list[float]:
        """Método auxiliar para  neighbor angles deg.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atom = self.model.get_atom(atom_id)
        origin = QPointF(atom.x, atom.y)
        angles: list[float] = []
        for bond in self.model.bonds.values():
            if bond.a1_id == atom_id:
                other = self.model.get_atom(bond.a2_id)
            elif bond.a2_id == atom_id:
                other = self.model.get_atom(bond.a1_id)
            else:
                continue
            if other.element == "H":
                continue
            angles.append(angle_deg(origin, QPointF(other.x, other.y)))
        return angles

    def _atom_hybridization(self, atom_id: int) -> str:
        """Método auxiliar para  atom hybridization.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        has_double = False
        has_triple = False
        has_aromatic = False
        for bond in self.model.bonds.values():
            if bond.a1_id != atom_id and bond.a2_id != atom_id:
                continue
            if bond.is_aromatic:
                has_aromatic = True
            order = bond.display_order if bond.display_order is not None else bond.order
            if order >= 3:
                has_triple = True
            elif order == 2:
                has_double = True
        if has_triple:
            return "sp"
        if has_double or has_aromatic:
            return "sp2"
        return "sp3"

    def _angle_of_vector(self, vx: float, vy: float) -> float:
        # Invert Y to convert from screen coordinates (Y down) to math coordinates (Y up)
        # This matches the convention used in geom.angle_deg and endpoint_from_angle_len
        """Método auxiliar para  angle of vector.

        Args:
            vx: Descripción del parámetro.
            vy: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        return math.degrees(math.atan2(-vy, vx))

    def _implicit_h_angles(self, atom_id: int, count: int) -> list[float]:
        """
        Calculate angles for implicit hydrogen atoms to be drawn.
        
        This function finds the best positions for implicit H atoms by:
        1. Finding all gaps between existing bonds
        2. Distributing H atoms into the largest gaps
        """
        if count <= 0:
            return []
        
        atom = self.model.get_atom(atom_id)
        origin = QPointF(atom.x, atom.y)
        
        # Collect all existing bond angles using angle_deg from geom.py
        # This function correctly handles screen coordinates (Y down -> math Y up)
        # Skip H atoms since we're computing where to place implicit H
        neighbor_angles: list[float] = []
        for bond in self.model.bonds.values():
            if bond.a1_id == atom_id:
                other = self.model.get_atom(bond.a2_id)
            elif bond.a2_id == atom_id:
                other = self.model.get_atom(bond.a1_id)
            else:
                continue
            # Skip existing hydrogen atoms
            if other.element == "H":
                continue
            # Use angle_deg from geom.py for consistent coordinate handling
            a = angle_deg(origin, QPointF(other.x, other.y))
            neighbor_angles.append(a)
        
        # No existing bonds - use default positions
        if not neighbor_angles:
            return [0.0, 120.0, 240.0, 300.0][:count]
        
        # Find gaps and place H atoms in the largest gaps
        return self._find_best_h_positions(neighbor_angles, count)

    def _find_best_h_positions(self, occupied_angles: list[float], count: int) -> list[float]:
        """
        Find the best positions for implicit H atoms using a bisector strategy.
        
        Args:
            occupied_angles: List of angles (in degrees) where bonds already exist
            count: Number of H atoms to place
            
        Returns:
            List of angles where H atoms should be placed
        """
        if count <= 0:
            return []
        
        if not occupied_angles:
            # Default positions (triangular/tetrahedral-ish)
            if count == 1: return [0.0]
            if count == 2: return [120.0, 240.0]
            return [0.0, 120.0, 240.0]
        
        # Local function to normalize angle to [0, 360) degrees
        def norm_deg(a: float) -> float:
            """Método auxiliar para norm deg.

            Args:
                a: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
            v = a % 360.0
            return v + 360.0 if v < 0 else v
        
        # Sort angles
        sorted_angles = sorted(norm_deg(a) for a in occupied_angles)
        n = len(sorted_angles)
        
        # Find the single largest gap
        max_gap = 0.0
        best_mid = 0.0
        
        for i in range(n):
            start = sorted_angles[i]
            end = sorted_angles[(i + 1) % n]
            
            if i == n - 1:
                gap_size = (360.0 - start) + end
            else:
                gap_size = end - start
            
            if gap_size > max_gap:
                max_gap = gap_size
                # Bisector of the gap
                best_mid = norm_deg(start + gap_size / 2.0)
        
        result_angles: list[float] = []
        
        # Distribute H atoms symmetrically around the bisector of the largest gap
        if count == 1:
            # Place directly in the middle
            result_angles.append(best_mid)
        elif count == 2:
            # Place in a "V" shape (±30 degrees from bisector)
            # This looks good for sp3 chains (tetrahedral projection)
            spread = 30.0
            if max_gap < 60.0: # If gap is very tight, reduce spread
                spread = max_gap / 3.0
            result_angles.append(norm_deg(best_mid - spread))
            result_angles.append(norm_deg(best_mid + spread))
        else:
            # Count >= 3 (e.g. terminal CH3)
            # Distribute in a fan shape (0, ±60 degrees) or similar
            # Standard "Trident" look: center, +60, -60
            result_angles.append(best_mid)
            result_angles.append(norm_deg(best_mid - 60.0))
            result_angles.append(norm_deg(best_mid + 60.0))
            
            # If we needed more than 3, we'd add more, but C is usually max 4 bonds
            if count > 3:
                # Fallback purely even distribution if strange valency
                start_fan = best_mid - (count - 1) * 30.0
                return [norm_deg(start_fan + i * 60.0) for i in range(count)]

        return result_angles

    def _clear_implicit_h_overlays(self, atom_ids: Optional[Iterable[int]] = None) -> None:
        """Método auxiliar para  clear implicit h overlays.

        Args:
            atom_ids: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if atom_ids is None:
            atom_ids = list(self._implicit_h_overlays.keys())
        for atom_id in list(atom_ids):
            overlays = self._implicit_h_overlays.pop(atom_id, [])
            for line_item, text_item in overlays:
                if line_item.scene() is self.scene:
                    self.scene.removeItem(line_item)
                if text_item.scene() is self.scene:
                    self.scene.removeItem(text_item)

    def _refresh_implicit_h_overlays(self, atom_ids: Optional[Iterable[int]] = None) -> None:
        """Método auxiliar para  refresh implicit h overlays.

        Args:
            atom_ids: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if atom_ids is None:
            atom_ids = list(self.atom_items.keys())
        self._clear_implicit_h_overlays(atom_ids)
        if not self.state.show_implicit_hydrogens:
            return
        for atom_id in atom_ids:
            atom = self.model.atoms.get(atom_id)
            if atom is None or atom.element == "H":
                continue
            font = self._label_font_for_atom(atom_id)
            count = self._implicit_hydrogen_count(atom_id, atom.element)
            if count <= 0:
                continue
            bond_length = max(1.0, self._local_bond_length(atom_id))
            angles = self._implicit_h_angles(atom_id, count)
            overlays: list[tuple[QGraphicsLineItem, QGraphicsTextItem]] = []
            for angle_value in angles:
                # Calculate full position for text center
                pos = endpoint_from_angle_len(QPointF(0.0, 0.0), angle_value, bond_length)
                rad = math.radians(angle_value)
                ux = math.cos(rad)
                uy = -math.sin(rad)
                
                # Setup text item first to get its size
                text_item = QGraphicsTextItem("H")
                text_item.setFont(font)
                text_item.setDefaultTextColor(QColor("#000000"))
                text_item.setZValue(10)
                text_item.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
                text_item.setData(IMPLICIT_H_OVERLAY_ANCHOR_ROLE, atom_id)
                text_item.setData(IMPLICIT_H_OVERLAY_ANGLE_ROLE, float(angle_value))
                text_item.setParentItem(self.atom_items[atom_id])
                
                rect = text_item.boundingRect()
                text_item.setPos(pos.x() - rect.width() / 2, pos.y() - rect.height() / 2)
                
                # Calculate shortened line end to avoid overlap with text
                # Shrink by roughly half the text width plus a small padding
                shrink = rect.width() / 2.0 + 3.0
                line_end_len = max(0.0, bond_length - shrink)

                # Also shrink the line start to avoid crossing atom labels (e.g., visible carbons).
                start_shrink = self._label_shrink_for_atom(atom_id, ux, uy)
                line_start_len = max(0.0, start_shrink)
                if line_end_len < line_start_len:
                    line_end_len = line_start_len
                line_start = endpoint_from_angle_len(QPointF(0.0, 0.0), angle_value, line_start_len)
                line_end = endpoint_from_angle_len(QPointF(0.0, 0.0), angle_value, line_end_len)
                
                line_item = QGraphicsLineItem(
                    line_start.x(),
                    line_start.y(),
                    line_end.x(),
                    line_end.y(),
                )
                pen = QPen(QColor(self.drawing_style.bond_color), self.drawing_style.stroke_px)
                pen.setCapStyle(self.drawing_style.cap_style)
                pen.setJoinStyle(self.drawing_style.join_style)
                line_item.setPen(pen)
                line_item.setZValue(-5)
                line_item.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
                line_item.setData(IMPLICIT_H_OVERLAY_ANCHOR_ROLE, atom_id)
                line_item.setData(IMPLICIT_H_OVERLAY_ANGLE_ROLE, float(angle_value))
                line_item.setParentItem(self.atom_items[atom_id])

                overlays.append((line_item, text_item))
            if overlays:
                self._implicit_h_overlays[atom_id] = overlays

    def _local_bond_length(self, atom_id: int) -> float:
        """Método auxiliar para  local bond length.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atom = self.model.get_atom(atom_id)
        if atom is None:
            return float(self.state.bond_length)
        lengths: list[float] = []
        for bond in self.model.bonds.values():
            if bond.a1_id == atom_id:
                other_id = bond.a2_id
            elif bond.a2_id == atom_id:
                other_id = bond.a1_id
            else:
                continue
            other = self.model.get_atom(other_id)
            if other is None:
                continue
            length = math.hypot(other.x - atom.x, other.y - atom.y)
            if length > 1e-6:
                lengths.append(length)
        if lengths:
            return sum(lengths) / len(lengths)
        return float(self.state.bond_length)

    def _update_bond_label_shrinks(self, atom_ids: set[int]) -> None:
        """Método auxiliar para  update bond label shrinks.

        Args:
            atom_ids: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        for bond in self.model.bonds.values():
            if bond.a1_id not in atom_ids and bond.a2_id not in atom_ids:
                continue
            item = self.bond_items.get(bond.id)
            if item is None:
                continue
            self._configure_bond_rendering(bond, item)
            atom1 = self.model.get_atom(bond.a1_id)
            atom2 = self.model.get_atom(bond.a2_id)
            dx = atom2.x - atom1.x
            dy = atom2.y - atom1.y
            length = math.hypot(dx, dy)
            if length <= 1e-6:
                item.set_label_shrink(0.0, 0.0)
                item.update_positions(atom1, atom2)
                continue
            ux = dx / length
            uy = dy / length
            shrink_start_extra = 0.0
            shrink_end_extra = 0.0
            if bond.style == BondStyle.WEDGE:
                # A wedge is narrow at `a1` and wide at `a2`; only the wide end
                # needs extra label clearance. Applying the same padding to the
                # tip creates visible gaps on labeled heteroatoms.
                shrink_end_extra = max(
                    0.0,
                    self._bond_render_width(bond) - self.drawing_style.stroke_px,
                )
            shrink_start = self._label_shrink_for_atom(
                bond.a1_id,
                ux,
                uy,
                extra_pad=shrink_start_extra,
            )
            shrink_end = self._label_shrink_for_atom(
                bond.a2_id,
                -ux,
                -uy,
                extra_pad=shrink_end_extra,
            )
            item.set_label_shrink(shrink_start, shrink_end)
            item.update_positions(atom1, atom2)
        self._refresh_implicit_h_overlays(atom_ids)

    def _label_shrink_for_atom(
        self,
        atom_id: int,
        ux: float,
        uy: float,
        *,
        extra_pad: float = 0.0,
    ) -> float:
        """Método auxiliar para  label shrink for atom.

        Args:
            atom_id: Descripción del parámetro.
            ux: Descripción del parámetro.
            uy: Descripción del parámetro.
            extra_pad: Padding adicional para geometrías con ancho real.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = self.atom_items.get(atom_id)
        atom = self.model.atoms.get(atom_id)
        if (
            atom is not None
            and bool(getattr(atom, "is_coordination_center", False))
            and not bool(getattr(atom, "sphere_transparent", False))
        ):
            if item is not None and hasattr(item, "_coordination_draw_radius"):
                try:
                    return max(0.0, float(item._coordination_draw_radius()) + 1.5)
                except Exception:
                    pass
            configured = getattr(atom, "sphere_radius", None)
            return max(0.0, float(configured) if configured is not None else 16.0)
        if item is None or not item.label.isVisible():
            return 0.0
        rect = item.label.mapRectToParent(item.label.boundingRect())
        if rect.isNull() or rect.width() <= 0 or rect.height() <= 0:
            return 0.0
        # Increased padding to clear label characters (especially "C")
        # 6.0px provides a comfortable margin for standard font sizes.
        pad = 6.0
        display_element = self._resolved_display_element(atom) if atom is not None else ""
        if atom is not None and display_element in ELEMENT_SYMBOLS and display_element not in {"C", "H"}:
            pad = 3.5
            if self._inline_label_hydrogen_count(atom_id, display_element) > 0:
                pad = 5.0
        pad += max(0.0, float(extra_pad))
        rect = rect.adjusted(-pad, -pad, pad, pad)
        distance = self._ray_ellipse_distance(rect, ux, uy)
        return distance if distance is not None else 0.0

    def _apply_text_settings(self, item: TextAnnotationItem, property_name: str = "all") -> None:
        """Apply text format settings to the item, ideally merging only what changed."""
        s = self._current_text_settings
        
        # Prepare a partial format for merging
        fmt = QTextCharFormat()
        
        if property_name in ("all", "family"):
            fmt.setFontFamilies([s["family"]])
        if property_name in ("all", "size"):
            fmt.setFontPointSize(s["size"])
        if property_name in ("all", "bold"):
            fmt.setFontWeight(QFont.Weight.Bold if s["bold"] else QFont.Weight.Normal)
        if property_name in ("all", "italic"):
            fmt.setFontItalic(s["italic"])
        if property_name in ("all", "underline"):
            fmt.setFontUnderline(s["underline"])
        if property_name in ("all", "sub", "sup"):
            if s["sub"]:
                fmt.setVerticalAlignment(QTextCharFormat.VerticalAlignment.AlignSubScript)
            elif s["sup"]:
                fmt.setVerticalAlignment(QTextCharFormat.VerticalAlignment.AlignSuperScript)
            else:
                fmt.setVerticalAlignment(QTextCharFormat.VerticalAlignment.AlignNormal)
        
        if property_name in ("all", "color"):
            fmt.setForeground(QBrush(s["color"]))
        
        # Get Cursor
        cursor = item.textCursor()
        is_editing = (item.textInteractionFlags() & Qt.TextInteractionFlag.TextEditorInteraction)

        if not is_editing:
            # Not editing -> Apply to whole document
            cursor.select(cursor.SelectionType.Document)
            # For whole document, also update default properties
            if property_name in ("all", "family", "size", "bold", "italic", "underline"):
                font = item.font()
                if property_name in ("all", "family"): font.setFamily(s["family"])
                if property_name in ("all", "size"): font.setPointSize(s["size"])
                if property_name in ("all", "bold"): font.setBold(s["bold"])
                if property_name in ("all", "italic"): font.setItalic(s["italic"])
                if property_name in ("all", "underline"): font.setUnderline(s["underline"])
                item.setFont(font)
            if property_name in ("all", "color"):
                item.setDefaultTextColor(s["color"])
        
        # Apply the format
        cursor.mergeCharFormat(fmt)
        if is_editing and property_name in ("sub", "sup") and cursor.hasSelection():
            end_pos = cursor.selectionEnd()
            cursor.setPosition(end_pos)
            cursor.setCharFormat(fmt)
        item.setTextCursor(cursor)
        if not is_editing:
            self._clear_text_cursor_selection(item)

    def update_text_format(self, family: str, size: int, b: bool, i: bool, u: bool, sub: bool, sup: bool, property_name: str = "all") -> None:
        """Update current text settings and apply to selected text items."""
        before_settings = dict(self._current_text_settings)
        self._current_text_settings.update({
            "family": family, "size": size,
            "bold": b, "italic": i, "underline": u,
            "sub": sub, "sup": sup
        })
        items_to_update = self._text_items_for_formatting()
        if not items_to_update:
            self.sync_text_selection_state()
            return
        self.undo_stack.push(
            FormatTextItemsCommand(
                self,
                items_to_update,
                before_settings,
                self._current_text_settings,
                property_name,
            )
        )

    def update_text_alignment(self, alignment: Qt.AlignmentFlag) -> None:
        """Actualiza text alignment.

        Args:
            alignment: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        items_to_update = self.scene.selectedItems()
        focus_item = self.scene.focusItem()
        if isinstance(focus_item, TextAnnotationItem) and focus_item not in items_to_update:
            items_to_update.append(focus_item)

        for item in items_to_update:
            if not isinstance(item, TextAnnotationItem):
                continue
            width = float(item.textWidth())
            if alignment != Qt.AlignmentFlag.AlignLeft and (
                not math.isfinite(width) or width <= 0.0
            ):
                item.setTextWidth(max(60.0, self._text_effective_width(item)))
            cursor = item.textCursor()
            is_editing = (
                item.textInteractionFlags()
                & Qt.TextInteractionFlag.TextEditorInteraction
            )
            if not is_editing:
                cursor.select(cursor.SelectionType.Document)
            block_format = cursor.blockFormat()
            block_format.setAlignment(alignment)
            cursor.setBlockFormat(block_format)
            item.setTextCursor(cursor)
            if not is_editing:
                self._clear_text_cursor_selection(item)
            doc = item.document()
            option = doc.defaultTextOption()
            option.setWrapMode(QTextOption.WrapMode.WrapAtWordBoundaryOrAnywhere)
            option.setAlignment(alignment)
            doc.setDefaultTextOption(option)
            item.update()
        self.sync_text_selection_state()

    def update_text_color(self, color: QColor) -> None:
        """Update current text color and apply to selection."""
        self._current_text_settings["color"] = color
        
        items_to_update = self.scene.selectedItems()
        focus_item = self.scene.focusItem()
        if isinstance(focus_item, TextAnnotationItem) and focus_item not in items_to_update:
            items_to_update.append(focus_item)

        for item in items_to_update:
            if isinstance(item, TextAnnotationItem):
                self._apply_text_settings(item, "color")
