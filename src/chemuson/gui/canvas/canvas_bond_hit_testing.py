from __future__ import annotations

import math
from typing import List, Optional, Tuple

from PyQt6.QtCore import QPointF

from chemuson.core.model import BondStereo, BondStyle
from chemuson.gui.commands import AddBondCommand
from chemuson.gui.geom import (
    angle_deg,
    closest_atom,
    closest_bond,
    endpoint_from_angle_len,
    segment_min_distance,
    segments_nearly_equal,
)

from .canvas_constants import (
    ATOM_HIT_RADIUS,
    BOND_OVERLAP_TOLERANCE_PX,
    HOVER_ATOM_RADIUS,
    HOVER_BOND_DISTANCE,
    IMPLICIT_H_OVERLAY_ANCHOR_ROLE,
    IMPLICIT_H_OVERLAY_ANGLE_ROLE,
)


class CanvasBondHitTestingMixin:
    def _pick_implicit_h_overlay(self, scene_pos: QPointF) -> tuple[Optional[int], Optional[float]]:
        best_anchor_id: Optional[int] = None
        best_angle: Optional[float] = None
        best_score = float("inf")
        line_tolerance = max(5.0, self.drawing_style.stroke_px + 3.0)

        def _parse_anchor(value) -> Optional[int]:
            try:
                return int(value)
            except (TypeError, ValueError):
                return None

        def _parse_angle(value) -> Optional[float]:
            try:
                return float(value)
            except (TypeError, ValueError):
                return None

        for _atom_id, overlays in self._implicit_h_overlays.items():
            for line_item, text_item in overlays:
                if text_item.scene() is self.scene:
                    text_rect = text_item.mapRectToScene(text_item.boundingRect()).adjusted(
                        -2.0,
                        -2.0,
                        2.0,
                        2.0,
                    )
                    if text_rect.contains(scene_pos):
                        anchor_id = _parse_anchor(text_item.data(IMPLICIT_H_OVERLAY_ANCHOR_ROLE))
                        if anchor_id is None or anchor_id not in self.model.atoms:
                            continue
                        text_center = text_item.mapToScene(text_item.boundingRect().center())
                        score = math.hypot(
                            scene_pos.x() - text_center.x(),
                            scene_pos.y() - text_center.y(),
                        )
                        if score <= best_score:
                            best_score = score
                            best_anchor_id = anchor_id
                            best_angle = _parse_angle(text_item.data(IMPLICIT_H_OVERLAY_ANGLE_ROLE))

                if line_item.scene() is not self.scene:
                    continue
                line = line_item.line()
                p1 = line_item.mapToScene(line.p1())
                p2 = line_item.mapToScene(line.p2())
                line_dist = segment_min_distance(scene_pos, scene_pos, p1, p2)
                if line_dist > line_tolerance or line_dist > best_score:
                    continue
                anchor_id = _parse_anchor(line_item.data(IMPLICIT_H_OVERLAY_ANCHOR_ROLE))
                if anchor_id is None or anchor_id not in self.model.atoms:
                    continue
                best_score = line_dist
                best_anchor_id = anchor_id
                line_angle = _parse_angle(line_item.data(IMPLICIT_H_OVERLAY_ANGLE_ROLE))
                if line_angle is None:
                    anchor_atom = self.model.get_atom(anchor_id)
                    line_angle = angle_deg(QPointF(anchor_atom.x, anchor_atom.y), p2)
                best_angle = line_angle

        return best_anchor_id, best_angle

    def _substitute_implicit_hydrogen(
        self,
        anchor_id: int,
        element: str,
        angle_value: Optional[float],
    ) -> bool:
        if anchor_id not in self.model.atoms:
            return False
        anchor_atom = self.model.get_atom(anchor_id)
        if anchor_atom.element == "H":
            return False
        if self._implicit_hydrogen_count(anchor_id, anchor_atom.element) <= 0:
            return False

        anchor_point = QPointF(anchor_atom.x, anchor_atom.y)
        if angle_value is None:
            direction = self._label_open_direction(anchor_id)
            if direction.manhattanLength() <= 1e-6:
                direction = QPointF(1.0, 0.0)
            target_probe = QPointF(anchor_point.x() + direction.x(), anchor_point.y() + direction.y())
            angle_value = angle_deg(anchor_point, target_probe)

        bond_length = max(1.0, self._local_bond_length(anchor_id))
        new_pos = endpoint_from_angle_len(anchor_point, angle_value, bond_length)
        cmd = AddBondCommand(
            self.model,
            self,
            anchor_id,
            None,
            1,
            BondStyle.PLAIN,
            BondStereo.NONE,
            is_aromatic=False,
            new_atom_element=element,
            new_atom_pos=(new_pos.x(), new_pos.y()),
        )
        self.undo_stack.push(cmd)
        return True

    def _pick_hover_target(self, scene_pos: QPointF) -> tuple[Optional[int], Optional[int]]:
        atoms = self._interactive_atom_candidates()
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
            (atom_id, atom_x, atom_y)
            for atom_id, atom_x, atom_y in self._interactive_atom_candidates()
            if atom_id not in excluded
        ]
        return closest_atom(pos, candidates, ATOM_HIT_RADIUS)

    def _build_ring_vertex_defs(
        self,
        vertices: List[QPointF],
        anchor_type: str,
        anchor_id: Optional[int],
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

        for vertex in remaining_vertices:
            snapped_id = self._snap_ring_vertex(vertex, used_ids)
            if snapped_id is not None:
                used_ids.add(snapped_id)
                vertex_defs.append((snapped_id, vertex.x(), vertex.y()))
            else:
                vertex_defs.append((None, vertex.x(), vertex.y()))

        return vertex_defs
