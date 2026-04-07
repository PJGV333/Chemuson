from __future__ import annotations

import math
from typing import Optional

from PyQt6.QtCore import QPointF, Qt

from chemuson.core.model import BondStereo, BondStyle
from chemuson.gui.commands import AddAtomCommand, AddBondCommand
from chemuson.gui.geom import angle_deg, closest_atom, endpoint_from_angle_len

from .canvas_constants import ATOM_HIT_RADIUS, OPTIMIZE_ZONE_SCALE


class CanvasBondDragMixin:
    def _cancel_drag(self) -> None:
        self._drag_mode = "none"
        self._drag_anchor = None
        self._flex_bond_pending = None
        self._ring_last_vertices = None
        self._chain_last_points = None
        self._drag_free_orientation = False
        self._bond_drag_start = None
        self._rotation_dragging = False
        self._scale_dragging = False
        self._rotation_center = None
        self._rotation_start_angle = None
        self._rotation_start_positions = {}
        self._rotation_start_arrow_positions = {}
        self._is_rotating_3d = False
        self._rotation_3d_start_view_pos = None
        self._rotation_3d_before_positions = {}
        self._rotation_3d_click_atom_id = None
        self._rotation_3d_has_moved = False
        self._rotation_3d_drag_start_pitch_deg = 0.0
        self._rotation_3d_drag_start_yaw_deg = 0.0
        self._drag_start_selection_bbox = None
        self._drag_affected_bond_ids = set()
        self._drag_affects_ring_centers = False
        self._release_interaction_mouse()
        self._clear_selection_drag()
        if self._overlays_ready:
            self._optimize_zone.hide_zone()
            self._preview_bond_item.hide_preview()
            self._preview_ring_item.hide_preview()
            self._preview_chain_item.hide_preview()
            self._preview_chain_label.hide_label()
            self._preview_arrow_item.hide_preview()

    def _on_undo_stack_changed(self, _index: int) -> None:
        self._cancel_drag()
        if self._force_full_scene_sync_on_undo_change:
            snapshot = self._selection_snapshot()
            self._force_full_scene_sync_on_undo_change = False
            self._sync_scene_with_model()
            self._restore_selection_snapshot(snapshot)
        self._update_hover(self._last_scene_pos)
        self._update_selection_overlay()

    def _begin_place_bond(self, anchor_id: Optional[int], scene_pos: QPointF) -> None:
        self._drag_mode = "place_bond"
        self._drag_free_orientation = False
        self._bond_drag_start = scene_pos
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
        if self._should_use_default_bond_angle(modifiers, scene_pos):
            p1 = self._compute_default_bond_endpoint(p0, self._drag_anchor["id"])
        else:
            p1 = self._compute_bond_endpoint(p0, scene_pos, modifiers)
        if self.state.active_bond_style == BondStyle.FLEX:
            curve_1, curve_2 = self._flex_curve_from_pointer(p0, p1, scene_pos, modifiers)
            self._preview_bond_item.update_flex_curve(p0, p1, curve_1, curve_2)
        else:
            self._preview_bond_item.update_line(p0, p1)

    def _flex_curve_from_pointer(
        self,
        start: QPointF,
        end: QPointF,
        pointer: QPointF,
        modifiers: Qt.KeyboardModifiers,
    ) -> tuple[float, float]:
        dx = end.x() - start.x()
        dy = end.y() - start.y()
        length = math.hypot(dx, dy)
        if length <= 1e-6:
            return 0.22, 0.22
        ux = dx / length
        uy = dy / length
        nx = -uy
        ny = ux
        mid_x = (start.x() + end.x()) * 0.5
        mid_y = (start.y() + end.y()) * 0.5
        signed = ((pointer.x() - mid_x) * nx + (pointer.y() - mid_y) * ny) / length
        magnitude = min(0.95, max(0.12, abs(signed)))
        sign = -1.0 if signed < 0.0 else 1.0
        if abs(signed) < 0.02:
            sign = 1.0
        c1 = sign * magnitude
        c2 = -c1 if modifiers & Qt.KeyboardModifier.ShiftModifier else c1
        return c1, c2

    def _begin_flex_adjust_mode(
        self,
        p0: QPointF,
        p1: QPointF,
        anchor_id: Optional[int],
        target_id: Optional[int],
        use_default: bool,
        modifiers: Qt.KeyboardModifiers,
    ) -> None:
        curve_1, curve_2 = self._flex_curve_from_pointer(p0, p1, self._last_scene_pos, modifiers)
        self._flex_bond_pending = {
            "p0": QPointF(p0),
            "p1": QPointF(p1),
            "anchor_id": anchor_id,
            "target_id": target_id,
            "use_default": bool(use_default),
            "curve_1": float(curve_1),
            "curve_2": float(curve_2),
        }
        self._drag_mode = "flex_adjust"
        self._drag_anchor = None
        self._drag_free_orientation = False
        self._bond_drag_start = None
        if self._overlays_ready:
            self._optimize_zone.hide_zone()
            self._preview_bond_item.update_flex_curve(p0, p1, curve_1, curve_2)

    def _update_flex_bond_preview(
        self,
        scene_pos: QPointF,
        modifiers: Qt.KeyboardModifiers,
    ) -> None:
        pending = self._flex_bond_pending
        if not pending:
            return
        p0 = pending["p0"]
        p1 = pending["p1"]
        curve_1, curve_2 = self._flex_curve_from_pointer(p0, p1, scene_pos, modifiers)
        pending["curve_1"] = float(curve_1)
        pending["curve_2"] = float(curve_2)
        if self._overlays_ready:
            self._preview_bond_item.update_flex_curve(p0, p1, curve_1, curve_2)

    def _commit_flex_bond(self) -> None:
        pending = self._flex_bond_pending
        if not pending:
            self._cancel_drag()
            return

        p0 = QPointF(pending["p0"])
        p1 = QPointF(pending["p1"])
        anchor_id = pending["anchor_id"]
        target_id = pending["target_id"]
        use_default = bool(pending["use_default"])
        curve_1 = float(pending["curve_1"])
        curve_2 = float(pending["curve_2"])

        order = 1
        style = BondStyle.FLEX
        stereo = BondStereo.NONE
        is_aromatic = False

        if anchor_id is None:
            self.undo_stack.beginMacro("Add flexible bond")
            anchor_cmd = AddAtomCommand(
                self.model,
                self,
                self.state.default_element,
                p0.x(),
                p0.y(),
                expected_bonds=1,
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
                    flex_curve_1=curve_1,
                    flex_curve_2=curve_2,
                )
                target_atom = self.model.get_atom(target_id)
                final_p1 = QPointF(target_atom.x, target_atom.y)
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
                    flex_curve_1=curve_1,
                    flex_curve_2=curve_2,
                )
                self.undo_stack.push(cmd)
                final_p1 = QPointF(p1)
            self.undo_stack.endMacro()
            self._update_bond_angle_state(p0, final_p1, use_default, anchor_id)
            self._cancel_drag()
            return

        if target_id is not None:
            self._create_or_update_bond(
                anchor_id,
                target_id,
                order,
                style,
                stereo,
                is_aromatic,
                flex_curve_1=curve_1,
                flex_curve_2=curve_2,
            )
            target_atom = self.model.get_atom(target_id)
            final_p1 = QPointF(target_atom.x, target_atom.y)
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
                flex_curve_1=curve_1,
                flex_curve_2=curve_2,
            )
            self.undo_stack.push(cmd)
            final_p1 = QPointF(p1)
        self._update_bond_angle_state(p0, final_p1, use_default, anchor_id)
        self._cancel_drag()

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
        use_default = self._should_use_default_bond_angle(modifiers, self._last_scene_pos)
        if use_default:
            p1 = self._compute_default_bond_endpoint(p0, anchor_id)
        else:
            p1 = self._compute_bond_endpoint(p0, self._last_scene_pos, modifiers)

        atom_positions = self._interactive_atom_candidates()
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

        final_p1 = p1
        if target_id is not None:
            target_atom = self.model.get_atom(target_id)
            final_p1 = QPointF(target_atom.x, target_atom.y)

        if style == BondStyle.FLEX:
            self._begin_flex_adjust_mode(p0, final_p1, anchor_id, target_id, use_default, modifiers)
            return

        if anchor_id is None:
            self.undo_stack.beginMacro("Add bond")
            anchor_cmd = AddAtomCommand(
                self.model,
                self,
                self.state.default_element,
                p0.x(),
                p0.y(),
                expected_bonds=1,
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
            self._update_bond_angle_state(p0, final_p1, use_default, anchor_id)
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
            self._update_bond_angle_state(p0, final_p1, use_default, anchor_id)
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
        free_drag = self._drag_mode == "place_bond" and self._drag_free_orientation
        if not self.state.fixed_angles or use_alt or free_drag:
            theta = cursor_theta
            length = self.state.bond_length if self.state.fixed_lengths and not use_shift else max(5.0, dist)
            p1 = endpoint_from_angle_len(anchor, theta, length)
            snap_id = closest_atom(p1, self._interactive_atom_candidates(), ATOM_HIT_RADIUS)
            if snap_id is not None and (self._drag_anchor is None or snap_id != self._drag_anchor["id"]):
                atom = self.model.get_atom(snap_id)
                return QPointF(atom.x, atom.y)
            return p1

        order = bond_order if bond_order is not None else self.state.active_bond_order
        aromatic_flag = is_aromatic if is_aromatic is not None else self.state.active_bond_aromatic
        style = self.state.active_bond_style if bond_order is None and is_aromatic is None else BondStyle.PLAIN
        if style != BondStyle.PLAIN or aromatic_flag:
            order = 1

        if use_optimize and self._drag_anchor["id"] is not None:
            cursor_theta = self._select_preferred_angle(
                self._drag_anchor["id"],
                cursor_theta,
                order,
                aromatic_flag,
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
            allow_length_boost=False,
        )
        p1 = endpoint_from_angle_len(anchor, theta, final_length)
        snap_id = closest_atom(p1, self._interactive_atom_candidates(), ATOM_HIT_RADIUS)
        if snap_id is not None and (self._drag_anchor is None or snap_id != self._drag_anchor["id"]):
            atom = self.model.get_atom(snap_id)
            return QPointF(atom.x, atom.y)
        return p1
