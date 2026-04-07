from __future__ import annotations

import math
from typing import Optional

from PyQt6.QtCore import QPointF, Qt

from chemuson.core.model import BondStereo, BondStyle
from chemuson.gui.commands import AddAtomCommand, AddBondCommand, ChangeBondCommand


class CanvasBondModelOpsMixin:
    def _infer_coordination_donor_atom(
        self,
        a1_id: int,
        a2_id: int,
        preferred: Optional[int] = None,
    ) -> int:
        if preferred in {a1_id, a2_id}:
            return int(preferred)
        atom1 = self.model.get_atom(a1_id)
        atom2 = self.model.get_atom(a2_id)
        a1_is_center = bool(getattr(atom1, "is_coordination_center", False))
        a2_is_center = bool(getattr(atom2, "is_coordination_center", False))
        if a1_is_center and not a2_is_center:
            return a2_id
        if a2_is_center and not a1_is_center:
            return a1_id
        return a1_id

    @staticmethod
    def _parse_bond_style_payload(bond_data: dict) -> BondStyle:
        raw_style = bond_data.get("style")
        if raw_style is None:
            raw_style = bond_data.get("type", BondStyle.PLAIN.value)
        try:
            return BondStyle(raw_style)
        except Exception:
            return BondStyle.PLAIN

    @staticmethod
    def _parse_bond_stereo_payload(bond_data: dict) -> BondStereo:
        raw_stereo = bond_data.get("stereo", BondStereo.NONE.value)
        try:
            return BondStereo(raw_stereo)
        except Exception:
            return BondStereo.NONE

    def _create_or_update_bond(
        self,
        a1_id: int,
        a2_id: int,
        order: int,
        style: BondStyle,
        stereo: BondStereo,
        is_aromatic: bool,
        flex_curve_1: Optional[float] = None,
        flex_curve_2: Optional[float] = None,
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
                if existing.is_aromatic and not is_aromatic and style != BondStyle.PLAIN:
                    new_is_aromatic = True
            new_donor_atom_id = None
            if new_style == BondStyle.COORDINATION:
                new_donor_atom_id = self._infer_coordination_donor_atom(
                    existing.a1_id,
                    existing.a2_id,
                    preferred=getattr(existing, "donor_atom_id", None),
                )
            new_flex_curve_1 = None
            new_flex_curve_2 = None
            if new_style == BondStyle.FLEX:
                if flex_curve_1 is not None or flex_curve_2 is not None:
                    new_flex_curve_1 = float(flex_curve_1) if flex_curve_1 is not None else None
                    new_flex_curve_2 = float(flex_curve_2) if flex_curve_2 is not None else None
                elif getattr(existing, "style", None) == BondStyle.FLEX:
                    new_flex_curve_1 = getattr(existing, "flex_curve_1", None)
                    new_flex_curve_2 = getattr(existing, "flex_curve_2", None)
            cmd = ChangeBondCommand(
                self.model,
                self,
                existing.id,
                new_order=new_order,
                new_style=new_style,
                new_stereo=new_stereo,
                new_is_aromatic=new_is_aromatic,
                new_donor_atom_id=new_donor_atom_id,
                new_flex_curve_1=new_flex_curve_1,
                new_flex_curve_2=new_flex_curve_2,
            )
            self.undo_stack.push(cmd)
            if new_is_aromatic:
                self._kekulize_aromatic_bonds(seed_atoms={a1_id, a2_id})
            return

        donor_atom_id = None
        if style == BondStyle.COORDINATION:
            donor_atom_id = self._infer_coordination_donor_atom(a1_id, a2_id, preferred=a1_id)
        cmd = AddBondCommand(
            self.model,
            self,
            a1_id,
            a2_id,
            order,
            style,
            stereo,
            is_aromatic=is_aromatic,
            donor_atom_id=donor_atom_id,
            flex_curve_1=flex_curve_1,
            flex_curve_2=flex_curve_2,
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
                new_donor_atom_id=None,
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
        if bond_id not in self.model.bonds:
            return
        bond = self.model.get_bond(bond_id)
        order = self.state.active_bond_order
        style = self.state.active_bond_style
        stereo = self.state.active_bond_stereo
        requested_aromatic = self.state.active_bond_aromatic
        is_aromatic = requested_aromatic
        if bond.is_aromatic and not requested_aromatic and style != BondStyle.PLAIN:
            is_aromatic = True
        if style != BondStyle.PLAIN:
            order = 1
        if requested_aromatic:
            order = 1
            style = BondStyle.PLAIN
            stereo = BondStereo.NONE
        donor_atom_id = None
        if style == BondStyle.COORDINATION:
            donor_atom_id = self._infer_coordination_donor_atom(
                bond.a1_id,
                bond.a2_id,
                preferred=getattr(bond, "donor_atom_id", None),
            )
        cmd = ChangeBondCommand(
            self.model,
            self,
            bond_id,
            new_order=order,
            new_style=style,
            new_stereo=stereo,
            new_is_aromatic=is_aromatic,
            new_donor_atom_id=donor_atom_id,
        )
        self.undo_stack.push(cmd)
        if is_aromatic:
            self._kekulize_aromatic_bonds(seed_atoms={bond.a1_id, bond.a2_id})

    def _create_first_bond(self, scene_pos: QPointF, modifiers: Qt.KeyboardModifiers) -> None:
        del modifiers
        self.undo_stack.beginMacro("Add bond")
        atom_cmd = AddAtomCommand(
            self.model,
            self,
            self.state.default_element,
            scene_pos.x(),
            scene_pos.y(),
            expected_bonds=1,
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
            near_anchor = dist < self.state.bond_length * 0.6
            aligned_with_existing = False
            if existing_angles:
                nearest = min(self._angle_distance(angle_from_click, angle) for angle in existing_angles)
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
        new_x = anchor.x + self.state.bond_length * math.cos(angle)
        new_y = anchor.y + self.state.bond_length * math.sin(angle)
        atom_cmd = AddAtomCommand(
            self.model,
            self,
            self.state.default_element,
            new_x,
            new_y,
            expected_bonds=1,
        )
        self.undo_stack.push(atom_cmd)
        new_atom_id = atom_cmd.atom_id
        if new_atom_id is None:
            return
        self._add_bond_between(anchor_id, new_atom_id)
        self._record_bond_angle(angle)
        self._set_bond_anchor(new_atom_id, reset_angle=False)

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
