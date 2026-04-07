from __future__ import annotations

import math
from typing import Iterable, Optional

from PyQt6.QtCore import QPointF, Qt

from chemuson.core.model import bond_is_structural
from chemuson.gui.geom import (
    angle_deg,
    angle_distance_deg,
    candidate_directions_deg,
    closest_atom,
    endpoint_from_angle_len,
    filter_occupied_angles_deg,
    geometry_for_bond,
    pick_closest_direction_deg,
    segment_min_distance,
    segments_intersect,
)

from .canvas_constants import (
    ANGLE_OCCUPIED_TOLERANCE_DEG,
    ATOM_HIT_RADIUS,
    BOND_DRAG_THRESHOLD_PX,
    COLLISION_LENGTH_BOOST,
    MIN_ATOM_DIST_SCALE,
    MIN_BOND_DIST_SCALE,
    SP3_OCCUPIED_TOLERANCE_DEG,
)


class CanvasBondGeometryMixin:
    def _get_anchor_bond_angles_deg(self, anchor_id: int) -> list[float]:
        angles = []
        anchor = self.model.get_atom(anchor_id)
        p0 = QPointF(anchor.x, anchor.y)
        for bond in self.model.bonds.values():
            if not bond_is_structural(bond):
                continue
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
            if not bond_is_structural(bond):
                continue
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
        atom_candidates = self._interactive_atom_candidates()
        target_id = closest_atom(p1, atom_candidates, ATOM_HIT_RADIUS)
        if target_id is not None:
            excluded_atoms.add(target_id)

        for atom_id, atom_x, atom_y in atom_candidates:
            if atom_id in excluded_atoms:
                continue
            dist = math.hypot(p1.x() - atom_x, p1.y() - atom_y)
            min_atom_dist = min(min_atom_dist, dist)
            if dist < atom_threshold:
                intersections += 1

        for bond in self.model.bonds.values():
            if bond.a1_id in excluded_atoms or bond.a2_id in excluded_atoms:
                continue
            a1 = self.model.get_atom(bond.a1_id)
            a2 = self.model.get_atom(bond.a2_id)
            p_a = QPointF(a1.x, a1.y)
            p_b = QPointF(a2.x, a2.y)
            if segments_intersect(p0, p1, p_a, p_b):
                intersections += 1
            bond_dist = segment_min_distance(p0, p1, p_a, p_b)
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

    def _sp3_congested_directions_deg(self, existing_angles_deg: Iterable[float]) -> list[float]:
        angles = sorted((angle % 360.0) for angle in existing_angles_deg)
        if len(angles) < 2:
            return []
        gap_midpoints: list[float] = []
        for idx, start in enumerate(angles):
            end = angles[(idx + 1) % len(angles)]
            gap = (end - start) % 360.0
            midpoint = (start + 0.5 * gap) % 360.0
            gap_midpoints.append(midpoint)
        seen: set[float] = set()
        deduped: list[float] = []
        for angle in gap_midpoints:
            key = round(angle, 6)
            if key in seen:
                continue
            seen.add(key)
            deduped.append(angle)
        return deduped

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
        sp3_exact_mode = geometry == "sp3" and anchor_id is not None and len(existing_angles) >= 3
        if sp3_exact_mode:
            candidates = self._sp3_congested_directions_deg(existing_angles)
        else:
            candidates = candidate_directions_deg(geometry, existing_angles, mouse_angle_deg)
            if self.state.fixed_angles:
                candidates = self._snap_angles_to_grid(candidates)
        occupied_tolerance = (
            SP3_OCCUPIED_TOLERANCE_DEG if sp3_exact_mode else ANGLE_OCCUPIED_TOLERANCE_DEG
        )
        candidates = filter_occupied_angles_deg(candidates, existing_angles, occupied_tolerance)
        if not candidates:
            if sp3_exact_mode:
                candidates = self._sp3_congested_directions_deg(existing_angles)
            if not candidates:
                candidates = candidate_directions_deg(geometry, [], mouse_angle_deg)
                if self.state.fixed_angles and not sp3_exact_mode:
                    candidates = self._snap_angles_to_grid(candidates)

        preferred: list[float] = []
        if geometry == "sp3" and anchor_id is not None:
            incoming = self._incoming_angle_deg(anchor_id)
            if incoming is not None:
                sp3_angle = self._sp3_display_angle_deg()
                preferred = [(incoming + sp3_angle) % 360.0, (incoming - sp3_angle) % 360.0]
                if self.state.fixed_angles and not sp3_exact_mode:
                    preferred = self._snap_angles_to_grid(preferred)

        if not apply_collisions:
            picked = pick_closest_direction_deg(candidates, mouse_angle_deg, preferred)
            return (picked if picked is not None else mouse_angle_deg), length

        length_candidates = (
            (length, length * COLLISION_LENGTH_BOOST) if allow_length_boost else (length,)
        )
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

    def _bond_drag_distance(self, scene_pos: QPointF) -> float:
        if self._bond_drag_start is None:
            return 0.0
        return (scene_pos - self._bond_drag_start).manhattanLength()

    def _update_bond_drag_state(self, scene_pos: QPointF) -> None:
        self._drag_free_orientation = (
            self._bond_drag_start is not None
            and self._bond_drag_distance(scene_pos) >= BOND_DRAG_THRESHOLD_PX
        )

    def _should_use_default_bond_angle(
        self, modifiers, scene_pos: Optional[QPointF] = None
    ) -> bool:
        if self._drag_mode != "place_bond":
            return False
        if modifiers & Qt.KeyboardModifier.AltModifier:
            return False
        if scene_pos is None:
            scene_pos = self._last_scene_pos
        if self._bond_drag_distance(scene_pos) >= BOND_DRAG_THRESHOLD_PX:
            return False
        return True

    def _compute_default_bond_endpoint(
        self, anchor: QPointF, anchor_id: Optional[int]
    ) -> QPointF:
        order, aromatic = self._current_bond_geometry()
        default_angle_deg = math.degrees(self._peek_default_bond_angle(anchor_id))
        theta, final_length = self._pick_bond_direction_deg(
            anchor,
            anchor_id,
            default_angle_deg,
            order,
            aromatic,
            self.state.bond_length,
            apply_collisions=True,
            allow_length_boost=False,
        )
        p1 = endpoint_from_angle_len(anchor, theta, final_length)
        snap_id = closest_atom(p1, self._interactive_atom_candidates(), ATOM_HIT_RADIUS)
        if snap_id is not None and (anchor_id is None or snap_id != anchor_id):
            if anchor_id is not None and self.model.find_bond_between(anchor_id, snap_id) is not None:
                return p1
            atom = self.model.get_atom(snap_id)
            return QPointF(atom.x, atom.y)
        return p1
