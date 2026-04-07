from __future__ import annotations

import math
from typing import List, Optional, Tuple

from PyQt6.QtCore import QPointF, Qt

from chemuson.core.model import BondStereo, BondStyle
from chemuson.gui.commands import AddChainCommand, AddRingCommand
from chemuson.gui.geom import (
    angle_deg,
    bond_side,
    choose_optimal_direction,
    endpoint_from_angle_len,
    snap_angle_deg,
)

from .canvas_constants import CHAIN_MAX_BONDS, OPTIMIZE_ZONE_SCALE

class CanvasToolsRingsChainsMixin:
    def _begin_place_ring(self, anchor_type: str, anchor_id: Optional[int], scene_pos: QPointF) -> None:
        """Método auxiliar para  begin place ring.

        Args:
            anchor_type: Descripción del parámetro.
            anchor_id: Descripción del parámetro.
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  update ring preview.

        Args:
            scene_pos: Descripción del parámetro.
            modifiers: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        vertices = self._compute_ring_vertices(scene_pos, modifiers)
        self._ring_last_vertices = vertices
        self._preview_ring_item.update_polygon(vertices)

    def _finalize_ring(self, modifiers: Qt.KeyboardModifiers) -> None:
        """Método auxiliar para  finalize ring.

        Args:
            modifiers: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  compute ring vertices.

        Args:
            scene_pos: Descripción del parámetro.
            modifiers: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  regular polygon from edge.

        Args:
            p1: Descripción del parámetro.
            p2: Descripción del parámetro.
            ring_size: Descripción del parámetro.
            direction: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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

    def _begin_place_chain(self, anchor_id: Optional[int], scene_pos: QPointF) -> None:
        """Método auxiliar para  begin place chain.

        Args:
            anchor_id: Descripción del parámetro.
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._drag_mode = "place_chain"
        if anchor_id is None:
            self._drag_anchor = {"type": "free", "id": None, "pos": QPointF(scene_pos)}
            self._optimize_zone.hide_zone()
        else:
            anchor = self.model.get_atom(anchor_id)
            self._drag_anchor = {"type": "atom", "id": anchor_id, "pos": QPointF(anchor.x, anchor.y)}
            radius = self.state.bond_length * OPTIMIZE_ZONE_SCALE
            self._optimize_zone.set_radius(radius)
            self._optimize_zone.update_center(anchor.x, anchor.y)
        self._update_chain_preview(scene_pos, Qt.KeyboardModifier.NoModifier)

    def _update_chain_preview(self, scene_pos: QPointF, modifiers: Qt.KeyboardModifiers) -> None:
        """Método auxiliar para  update chain preview.

        Args:
            scene_pos: Descripción del parámetro.
            modifiers: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._drag_anchor is None:
            return
        anchor_id = self._drag_anchor["id"]
        if anchor_id is None:
            p0 = QPointF(self._drag_anchor["pos"])
        else:
            anchor = self.model.get_atom(anchor_id)
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

        if not self.state.fixed_angles or use_alt:
            first_angle = raw_angle
        else:
            cursor_angle = raw_angle
            if anchor_id is None:
                # En cadena libre, usar el eje del cursor para permitir un zig-zag
                # "recto" (estilo ChemDraw) en cualquier dirección.
                step_deg = float(self.state.angle_step_deg)
                first_angle = snap_angle_deg(cursor_angle, step_deg) if step_deg > 0 and not use_shift else cursor_angle
            else:
                if use_optimize:
                    cursor_angle = self._select_preferred_angle(anchor_id, cursor_angle, 1, False)
                first_angle, _ = self._pick_bond_direction_deg(
                    p0,
                    anchor_id,
                    cursor_angle,
                    1,
                    False,
                    self.state.bond_length,
                    apply_collisions=True,
                    allow_length_boost=self.state.fixed_lengths and not use_shift,
                )

        n = max(1, min(CHAIN_MAX_BONDS, int(round(dist / self.state.bond_length))))
        points = [p0]
        # La herramienta de cadena siempre debe zigzaguear (si ángulos fijos están activos),
        # incluso cuando el ancla está en centros sp2/aromáticos.
        zigzag = not use_alt and self.state.fixed_angles
        sp3_angle = self._sp3_display_angle_deg()
        zigzag_turn = 180.0 - sp3_angle
        turn_pref = ((raw_angle - first_angle + 180.0) % 360.0) - 180.0
        zigzag_sign = 1.0 if turn_pref >= 0.0 else -1.0
        first_step_angle = first_angle
        second_step_angle = first_angle + zigzag_sign * zigzag_turn
        if zigzag and anchor_id is None:
            # Cadena libre: alternar simétricamente alrededor del eje pedido
            # por el cursor para formar una línea zig-zag recta.
            half_turn = zigzag_turn * 0.5
            axis_angle = first_angle
            first_step_angle = axis_angle - zigzag_sign * half_turn
            second_step_angle = axis_angle + zigzag_sign * half_turn
        current = p0
        for i in range(1, n + 1):
            if zigzag and i % 2 == 0:
                step_angle = second_step_angle
            else:
                step_angle = first_step_angle
            next_point = endpoint_from_angle_len(current, step_angle, self.state.bond_length)
            points.append(next_point)
            current = next_point

        self._chain_last_points = points
        self._preview_chain_item.update_polyline(points)
        self._preview_chain_label.update_label(str(n), points[-1] + QPointF(0, -10))

    def _finalize_chain(self, modifiers: Qt.KeyboardModifiers) -> None:
        """Método auxiliar para  finalize chain.

        Args:
            modifiers: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if not self._chain_last_points or self._drag_anchor is None:
            self._cancel_drag()
            return
        anchor_id = self._drag_anchor["id"]
        new_positions = [(p.x(), p.y()) for p in self._chain_last_points[1:]]
        if not new_positions:
            self._cancel_drag()
            return
        anchor_pos = None
        if anchor_id is None:
            anchor_point = QPointF(self._drag_anchor["pos"])
            anchor_pos = (anchor_point.x(), anchor_point.y())
        cmd = AddChainCommand(
            self.model,
            self,
            anchor_id,
            new_positions,
            element=self.state.default_element,
            anchor_position=anchor_pos,
        )
        self.undo_stack.push(cmd)
        self._cancel_drag()

    def _create_chain_hotkey(self, anchor_id: int, n: int) -> None:
        """Método auxiliar para  create chain hotkey.

        Args:
            anchor_id: Descripción del parámetro.
            n: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  create ring hotkey.

        Args:
            anchor_id: Descripción del parámetro.
            ring_size: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        anchor = self.model.get_atom(anchor_id)
        p0 = QPointF(anchor.x, anchor.y)
        angles = self._get_anchor_bond_angles_deg(anchor_id)
        theta = choose_optimal_direction(angles)
        p1 = endpoint_from_angle_len(p0, theta, self.state.bond_length)
        direction = bond_side(p0, p1, self._last_scene_pos)
        vertices = self._regular_polygon_from_edge(p0, p1, ring_size, direction)
        self._commit_ring_vertices(vertices, anchor_type="atom", anchor_id=anchor_id)

    def _create_ring_from_bond_hotkey(self, bond_id: int, ring_size: int) -> None:
        """Método auxiliar para  create ring from bond hotkey.

        Args:
            bond_id: Descripción del parámetro.
            ring_size: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  commit ring vertices.

        Args:
            vertices: Descripción del parámetro.
            anchor_type: Descripción del parámetro.
            anchor_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
