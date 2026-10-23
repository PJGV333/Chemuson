from __future__ import annotations

import math
from typing import Iterable, Optional

from PyQt6.QtCore import QPointF, Qt

from chemuson.core.model import bond_affects_valence, bond_is_structural
from chemuson.gui.geom import (
    SP3_SKELETAL_ANGLE_2D_DEG,
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

    def _nitrogen_geometry_context(self, atom_id: int) -> str | None:
        """Classifica N por conectividade covalente, ordens e aromaticidade."""
        atom = self.model.get_atom(atom_id)
        if atom is None or atom.element != "N":
            return None
        if bool(getattr(atom, "is_coordination_center", False)):
            return "coordination"

        bonds = [
            bond
            for bond in self.model.bonds.values()
            if bond_affects_valence(bond) and atom_id in (bond.a1_id, bond.a2_id)
        ]
        if int(getattr(atom, "formal_charge", 0) or 0) != 0:
            return "other"
        if len(bonds) > 3:
            return "other"
        if any(bond.is_aromatic for bond in bonds):
            return "aromatic"

        carbonyl_carbons: set[int] = set()
        for bond in bonds:
            if bond.order != 1 or bond.is_aromatic:
                continue
            neighbor_id = bond.a2_id if bond.a1_id == atom_id else bond.a1_id
            neighbor = self.model.get_atom(neighbor_id)
            if neighbor is None or neighbor.element != "C":
                continue
            has_carbonyl = False
            for candidate in self.model.bonds.values():
                if (
                    candidate.a1_id != neighbor_id
                    and candidate.a2_id != neighbor_id
                ) or not bond_affects_valence(candidate):
                    continue
                hetero_id = (
                    candidate.a2_id
                    if candidate.a1_id == neighbor_id
                    else candidate.a1_id
                )
                hetero = self.model.get_atom(hetero_id)
                if (
                    candidate.order >= 2
                    and not candidate.is_aromatic
                    and hetero is not None
                    and hetero.element in {"O", "S"}
                ):
                    has_carbonyl = True
                    break
            if has_carbonyl:
                carbonyl_carbons.add(neighbor_id)

        for carbonyl_id in carbonyl_carbons:
            attached_nitrogens = 0
            for candidate in self.model.bonds.values():
                if (
                    candidate.a1_id != carbonyl_id
                    and candidate.a2_id != carbonyl_id
                ):
                    continue
                if (
                    candidate.order != 1
                    or candidate.is_aromatic
                    or not bond_affects_valence(candidate)
                ):
                    continue
                neighbor_id = candidate.a2_id if candidate.a1_id == carbonyl_id else candidate.a1_id
                neighbor = self.model.get_atom(neighbor_id)
                if neighbor is not None and neighbor.element == "N":
                    attached_nitrogens += 1
            if attached_nitrogens >= 2:
                return "urea"
            return "amide"

        for bond in bonds:
            neighbor_id = bond.a2_id if bond.a1_id == atom_id else bond.a1_id
            neighbor = self.model.get_atom(neighbor_id)
            if neighbor is not None and neighbor.element not in {"C", "H"}:
                return "other"
        return "amine"

    def _nitrogen_context_dependent_refresh_ids(
        self, endpoint_a_id: int, endpoint_b_id: int
    ) -> list[int]:
        """N adyacentes a los extremos de un enlace cuyo contexto puede cambiar.

        Al crear, actualizar o eliminar un enlace, la clasificación
        amida/urea/amina de los N unidos a esos extremos puede variar sin
        que se refresquen sus etiquetas; estos IDs complementan el refresco
        de los propios extremos.
        """
        affected: list[int] = []
        for endpoint_id in (endpoint_a_id, endpoint_b_id):
            for other_bond in self.model.bonds.values():
                if not bond_affects_valence(other_bond):
                    continue
                if endpoint_id not in (other_bond.a1_id, other_bond.a2_id):
                    continue
                other_id = (
                    other_bond.a2_id
                    if other_bond.a1_id == endpoint_id
                    else other_bond.a1_id
                )
                other_atom = self.model.get_atom(other_id)
                if other_atom is not None and other_atom.element == "N":
                    affected.append(other_id)
        return list(dict.fromkeys(affected))

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
        nitrogen_context = self._nitrogen_geometry_context(anchor_id)
        if nitrogen_context in {"amide", "urea", "aromatic"}:
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
        apply_angle_grid: bool = True,
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
            candidates = candidate_directions_deg(
                geometry,
                existing_angles,
                mouse_angle_deg,
                sp3_angle_deg=SP3_SKELETAL_ANGLE_2D_DEG,
            )
            if self.state.fixed_angles and apply_angle_grid:
                candidates = self._snap_angles_to_grid(candidates)
        occupied_tolerance = (
            SP3_OCCUPIED_TOLERANCE_DEG if sp3_exact_mode else ANGLE_OCCUPIED_TOLERANCE_DEG
        )
        candidates = filter_occupied_angles_deg(candidates, existing_angles, occupied_tolerance)
        if not candidates:
            if sp3_exact_mode:
                candidates = self._sp3_congested_directions_deg(existing_angles)
            if not candidates:
                candidates = candidate_directions_deg(
                    geometry,
                    [],
                    mouse_angle_deg,
                    sp3_angle_deg=SP3_SKELETAL_ANGLE_2D_DEG,
                )
                if self.state.fixed_angles and not sp3_exact_mode and apply_angle_grid:
                    candidates = self._snap_angles_to_grid(candidates)

        preferred: list[float] = []
        if geometry == "sp3" and anchor_id is not None:
            incoming = self._incoming_angle_deg(anchor_id)
            if incoming is not None:
                sp3_angle = self._sp3_skeletal_angle_deg()
                preferred = [(incoming + sp3_angle) % 360.0, (incoming - sp3_angle) % 360.0]
                if self.state.fixed_angles and not sp3_exact_mode and apply_angle_grid:
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
            apply_angle_grid=False,
        )
        p1 = endpoint_from_angle_len(anchor, theta, final_length)
        snap_id = closest_atom(p1, self._interactive_atom_candidates(), ATOM_HIT_RADIUS)
        if snap_id is not None and (anchor_id is None or snap_id != anchor_id):
            if anchor_id is not None and self.model.find_bond_between(anchor_id, snap_id) is not None:
                return p1
            atom = self.model.get_atom(snap_id)
            return QPointF(atom.x, atom.y)
        return p1
