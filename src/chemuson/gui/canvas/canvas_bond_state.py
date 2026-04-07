from __future__ import annotations

import math
from typing import Iterable, List, Optional

from chemuson.core.model import BondStyle, bond_is_structural
from chemuson.gui.geom import SP3_BOND_ANGLE_DEG, snap_angle_deg

from .canvas_constants import BOND_LAST_ANGLE_TOLERANCE_DEG


class CanvasBondStateMixin:
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
        best = min(candidates, key=lambda value: self._angle_distance(value, angle))
        return self._normalize_angle(best)

    def _snap_angle(self, angle: float, step: float) -> float:
        if step <= 0:
            return angle
        return self._normalize_angle(round(angle / step) * step)

    def _bond_environment_step(self, anchor_id: int) -> float:
        max_order = max(1, self.state.active_bond_order)
        for bond in self.model.bonds.values():
            if not bond_is_structural(bond):
                continue
            if bond.a1_id == anchor_id or bond.a2_id == anchor_id:
                order = 2 if bond.is_aromatic else bond.order
                max_order = max(max_order, order)
        if max_order >= 3:
            return math.pi
        if max_order == 2:
            return 2 * math.pi / 3
        if self.state.fixed_angles:
            step_deg = self._angle_snap_step_deg()
            if step_deg > 0:
                return math.radians(snap_angle_deg(SP3_BOND_ANGLE_DEG, step_deg))
        return math.radians(self._sp3_display_angle_deg())

    def _get_anchor_bond_angles(self, anchor_id: int) -> List[float]:
        angles = []
        anchor = self.model.get_atom(anchor_id)
        for bond in self.model.bonds.values():
            if not bond_is_structural(bond):
                continue
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
            sep = min(self._angle_distance(candidate, angle) for angle in existing)
            if sep > best_sep:
                best_sep = sep
                best = candidate
        return self._normalize_angle(best) if best is not None else 0.0

    def _normalize_angle(self, angle: float) -> float:
        return (angle + math.pi * 2) % (math.pi * 2)

    def _angle_distance(self, a: float, b: float) -> float:
        diff = (a - b + math.pi) % (2 * math.pi) - math.pi
        return abs(diff)

    def _angle_snap_step_deg(self) -> float:
        if not self.state.fixed_angles:
            return 0.0
        step = float(self.state.angle_step_deg)
        return step if step > 0 else 0.0

    def _snap_angles_to_grid(self, angles: Iterable[float]) -> list[float]:
        step = self._angle_snap_step_deg()
        if step <= 0:
            return list(angles)
        snapped = [snap_angle_deg(angle, step) for angle in angles]
        seen: set[float] = set()
        deduped: list[float] = []
        for angle in snapped:
            key = round(angle, 6)
            if key in seen:
                continue
            seen.add(key)
            deduped.append(angle)
        return deduped

    def _sp3_display_angle_deg(self) -> float:
        return SP3_BOND_ANGLE_DEG

    def _current_bond_geometry(self) -> tuple[int, bool]:
        order = self.state.active_bond_order
        is_aromatic = self.state.active_bond_aromatic
        if self.state.active_bond_style != BondStyle.PLAIN or is_aromatic:
            order = 1
        return order, is_aromatic

    def _peek_default_bond_angle(self, anchor_id: Optional[int]) -> float:
        if anchor_id is None:
            return 0.0
        existing = self._get_anchor_bond_angles(anchor_id)
        step = self._bond_environment_step(anchor_id)
        if self._bond_last_angle is not None and existing:
            tolerance = math.radians(BOND_LAST_ANGLE_TOLERANCE_DEG)
            if min(self._angle_distance(self._bond_last_angle, angle) for angle in existing) <= tolerance:
                angle = self._bond_last_angle + step * self._bond_zigzag_sign
                return self._normalize_angle(angle)
        if not existing:
            return 0.0
        candidates = []
        for base in existing:
            candidates.append(base + step)
            candidates.append(base - step)
        return self._best_separated_angle(candidates, existing)

    def _update_bond_angle_state(
        self,
        p0,
        p1,
        used_default: bool,
        anchor_id: Optional[int],
    ) -> None:
        angle = math.atan2(p1.y() - p0.y(), p1.x() - p0.x())
        self._bond_last_angle = self._normalize_angle(angle)
        if not used_default:
            self._bond_zigzag_sign = 1
            return
        step = (
            self._bond_environment_step(anchor_id)
            if anchor_id is not None
            else math.radians(self._sp3_display_angle_deg())
        )
        if step < math.pi:
            self._bond_zigzag_sign *= -1
