from __future__ import annotations

from ._shared import (
    ANGLE_OCCUPIED_TOLERANCE_DEG,
    ATOM_HIT_RADIUS,
    AddAtomCommand,
    AddBondCommand,
    BOND_DRAG_THRESHOLD_PX,
    BOND_LAST_ANGLE_TOLERANCE_DEG,
    BOND_OVERLAP_TOLERANCE_PX,
    BondStereo,
    BondStyle,
    COLLISION_LENGTH_BOOST,
    ChangeBondCommand,
    HOVER_ATOM_RADIUS,
    HOVER_BOND_DISTANCE,
    IMPLICIT_H_OVERLAY_ANCHOR_ROLE,
    IMPLICIT_H_OVERLAY_ANGLE_ROLE,
    Iterable,
    List,
    MIN_ATOM_DIST_SCALE,
    MIN_BOND_DIST_SCALE,
    OPTIMIZE_ZONE_SCALE,
    Optional,
    QPointF,
    Qt,
    SP3_BOND_ANGLE_DEG,
    SP3_OCCUPIED_TOLERANCE_DEG,
    Tuple,
    angle_deg,
    angle_distance_deg,
    bond_is_structural,
    candidate_directions_deg,
    closest_atom,
    closest_bond,
    endpoint_from_angle_len,
    filter_occupied_angles_deg,
    geometry_for_bond,
    math,
    pick_closest_direction_deg,
    segment_min_distance,
    segments_intersect,
    segments_nearly_equal,
    snap_angle_deg,
)

class CanvasToolsBondingMixin:
    def _set_bond_anchor(self, atom_id: int, reset_angle: bool = False) -> None:
        """Método auxiliar para  set bond anchor.

        Args:
            atom_id: Descripción del parámetro.
            reset_angle: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  clear bond anchor.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self.bond_anchor_id is not None and self.bond_anchor_id in self.atom_items:
            if self.bond_anchor_id not in self.state.selected_atoms:
                self.atom_items[self.bond_anchor_id].set_selected(False)
        self.bond_anchor_id = None
        self._bond_last_angle = None
        self._bond_zigzag_sign = 1

    def _infer_coordination_donor_atom(
        self,
        a1_id: int,
        a2_id: int,
        preferred: Optional[int] = None,
    ) -> int:
        """Infere el átomo donador para enlaces coordinativos."""
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
        # Si no hay centro claro, conservar dirección del primer click.
        return a1_id

    @staticmethod
    def _parse_bond_style_payload(bond_data: dict) -> BondStyle:
        """Resuelve estilo de enlace para payloads nuevos y legados."""
        raw_style = bond_data.get("style")
        if raw_style is None:
            raw_style = bond_data.get("type", BondStyle.PLAIN.value)
        try:
            return BondStyle(raw_style)
        except Exception:
            return BondStyle.PLAIN

    @staticmethod
    def _parse_bond_stereo_payload(bond_data: dict) -> BondStereo:
        """Resuelve estereo de enlace para payloads con fallback seguro."""
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
        """Método auxiliar para  create or update bond.

        Args:
            a1_id: Descripción del parámetro.
            a2_id: Descripción del parámetro.
            order: Descripción del parámetro.
            style: Descripción del parámetro.
            stereo: Descripción del parámetro.
            is_aromatic: Descripción del parámetro.
            flex_curve_1: Curvatura normalizada de control 1 (estilo FLEX).
            flex_curve_2: Curvatura normalizada de control 2 (estilo FLEX).

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
                    # Mantener aromaticidad al aplicar estilos visuales no planos
                    # (p. ej. bold para proyecciones tipo Haworth).
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
                    new_flex_curve_1 = (
                        float(flex_curve_1) if flex_curve_1 is not None else None
                    )
                    new_flex_curve_2 = (
                        float(flex_curve_2) if flex_curve_2 is not None else None
                    )
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
        """Método auxiliar para  add bond between.

        Args:
            a1_id: Descripción del parámetro.
            a2_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        order = self.state.active_bond_order
        style = self.state.active_bond_style
        stereo = self.state.active_bond_stereo
        is_aromatic = self.state.active_bond_aromatic
        if style != BondStyle.PLAIN or is_aromatic:
            order = 1
        self._create_or_update_bond(a1_id, a2_id, order, style, stereo, is_aromatic)

    def _cycle_bond_order(self, bond_id: int) -> None:
        """Método auxiliar para  cycle bond order.

        Args:
            bond_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  apply bond style.

        Args:
            bond_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if bond_id not in self.model.bonds:
            return
        bond = self.model.get_bond(bond_id)
        order = self.state.active_bond_order
        style = self.state.active_bond_style
        stereo = self.state.active_bond_stereo
        requested_aromatic = self.state.active_bond_aromatic
        is_aromatic = requested_aromatic
        if bond.is_aromatic and not requested_aromatic and style != BondStyle.PLAIN:
            # Preservar aromaticidad si solo se cambia el estilo visual.
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
        """Método auxiliar para  create first bond.

        Args:
            scene_pos: Descripción del parámetro.
            modifiers: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  extend bond from anchor.

        Args:
            scene_pos: Descripción del parámetro.
            modifiers: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  add bond with new atom.

        Args:
            anchor_id: Descripción del parámetro.
            angle: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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

    def _record_bond_angle_between(self, a1_id: int, a2_id: int) -> None:
        """Método auxiliar para  record bond angle between.

        Args:
            a1_id: Descripción del parámetro.
            a2_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        a1 = self.model.get_atom(a1_id)
        a2 = self.model.get_atom(a2_id)
        angle = math.atan2(a2.y - a1.y, a2.x - a1.x)
        self._record_bond_angle(angle)

    def _record_bond_angle(self, angle: float) -> None:
        """Método auxiliar para  record bond angle.

        Args:
            angle: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._bond_last_angle = self._normalize_angle(angle)

    def _default_bond_angle(self, anchor_id: int) -> float:
        """Método auxiliar para  default bond angle.

        Args:
            anchor_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  reset zigzag.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._bond_zigzag_sign = 1

    def _snap_angle_to_environment(self, angle: float, anchor_id: int, step: float) -> float:
        """Método auxiliar para  snap angle to environment.

        Args:
            angle: Descripción del parámetro.
            anchor_id: Descripción del parámetro.
            step: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  snap angle.

        Args:
            angle: Descripción del parámetro.
            step: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if step <= 0:
            return angle
        return self._normalize_angle(round(angle / step) * step)

    def _bond_environment_step(self, anchor_id: int) -> float:
        """Método auxiliar para  bond environment step.

        Args:
            anchor_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  get anchor bond angles.

        Args:
            anchor_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  best separated angle.

        Args:
            candidates: Descripción del parámetro.
            existing: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        best = None
        best_sep = -1.0
        for candidate in candidates:
            sep = min(self._angle_distance(candidate, a) for a in existing)
            if sep > best_sep:
                best_sep = sep
                best = candidate
        return self._normalize_angle(best) if best is not None else 0.0

    def _normalize_angle(self, angle: float) -> float:
        """Método auxiliar para  normalize angle.

        Args:
            angle: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        return (angle + math.pi * 2) % (math.pi * 2)

    def _angle_distance(self, a: float, b: float) -> float:
        """Método auxiliar para  angle distance.

        Args:
            a: Descripción del parámetro.
            b: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        diff = (a - b + math.pi) % (2 * math.pi) - math.pi
        return abs(diff)

    def _angle_snap_step_deg(self) -> float:
        """Método auxiliar para  angle snap step deg.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if not self.state.fixed_angles:
            return 0.0
        step = float(self.state.angle_step_deg)
        return step if step > 0 else 0.0

    def _snap_angles_to_grid(self, angles: Iterable[float]) -> list[float]:
        """Método auxiliar para  snap angles to grid.

        Args:
            angles: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  sp3 display angle deg.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        # Mantener 109.5° real para geometría sp3.
        # Si se ajusta a una retícula de 30°, se vuelve 120° (trigonal) y
        # termina bloqueando el 4º enlace en carbonos tetravalentes.
        return SP3_BOND_ANGLE_DEG

    def _pick_implicit_h_overlay(self, scene_pos: QPointF) -> tuple[Optional[int], Optional[float]]:
        """Devuelve el H implÃ­cito dibujado bajo el cursor (si existe)."""
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
                        anchor_id = _parse_anchor(
                            text_item.data(IMPLICIT_H_OVERLAY_ANCHOR_ROLE)
                        )
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
                            best_angle = _parse_angle(
                                text_item.data(IMPLICIT_H_OVERLAY_ANGLE_ROLE)
                            )

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
        """Convierte un H implÃ­cito visual en un Ã¡tomo real enlazado."""
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
        """Método auxiliar para  pick hover target.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  update hover.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  cancel drag.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  on undo stack changed.

        Args:
            _index: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._cancel_drag()
        if self._force_full_scene_sync_on_undo_change:
            snapshot = self._selection_snapshot()
            self._force_full_scene_sync_on_undo_change = False
            self._sync_scene_with_model()
            self._restore_selection_snapshot(snapshot)
        self._update_hover(self._last_scene_pos)
        self._update_selection_overlay()

    def _get_anchor_bond_angles_deg(self, anchor_id: int) -> list[float]:
        """Método auxiliar para  get anchor bond angles deg.

        Args:
            anchor_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  bond geometry.

        Args:
            anchor_id: Descripción del parámetro.
            bond_order: Descripción del parámetro.
            is_aromatic: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  incoming angle deg.

        Args:
            anchor_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        angles = self._get_anchor_bond_angles_deg(anchor_id)
        if len(angles) == 1:
            return angles[0]
        return None

    def _direction_collision_metrics(
        self, anchor_id: Optional[int], p0: QPointF, p1: QPointF
    ) -> tuple[int, float, float]:
        """Método auxiliar para  direction collision metrics.

        Args:
            anchor_id: Descripción del parámetro.
            p0: Descripción del parámetro.
            p1: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        intersections = 0
        min_atom_dist = float("inf")
        min_bond_dist = float("inf")
        atom_threshold = self.state.bond_length * MIN_ATOM_DIST_SCALE
        bond_threshold = self.state.bond_length * MIN_BOND_DIST_SCALE
        excluded_atoms: set[int] = set()
        if anchor_id is not None:
            excluded_atoms.add(anchor_id)
        atom_candidates = self._interactive_atom_candidates()
        target_id = closest_atom(
            p1,
            atom_candidates,
            ATOM_HIT_RADIUS,
        )
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
        """Método auxiliar para  direction score.

        Args:
            anchor_id: Descripción del parámetro.
            p0: Descripción del parámetro.
            p1: Descripción del parámetro.
            angle_deg_value: Descripción del parámetro.
            mouse_angle_deg: Descripción del parámetro.
            preferred_deg: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  direction is valid.

        Args:
            anchor_id: Descripción del parámetro.
            p0: Descripción del parámetro.
            p1: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Propone direcciones para centros sp3 congestionados (>=3 enlaces).

        Con ±109.5° sobre cada enlace existente se generan candidatos casi
        paralelos (separaciones ~10.5°). Aquí usamos bisectrices de todos los
        huecos angulares y dejamos que el score (cursor/colisiones) elija.
        """
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
        """Método auxiliar para  pick bond direction deg.

        Args:
            anchor: Descripción del parámetro.
            anchor_id: Descripción del parámetro.
            mouse_angle_deg: Descripción del parámetro.
            bond_order: Descripción del parámetro.
            is_aromatic: Descripción del parámetro.
            length: Descripción del parámetro.
            apply_collisions: Descripción del parámetro.
            allow_length_boost: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        existing_angles = self._get_anchor_bond_angles_deg(anchor_id) if anchor_id else []
        geometry = (
            self._bond_geometry(anchor_id, bond_order, is_aromatic)
            if anchor_id is not None
            else geometry_for_bond(bond_order, is_aromatic, [])
        )
        # Solo forzamos geometría tetraédrica exacta (109.5°) cuando el átomo
        # ya está congestionado (p. ej., intentando crear el 4º enlace).
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
        candidates = filter_occupied_angles_deg(
            candidates, existing_angles, occupied_tolerance
        )
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
                preferred = [
                    (incoming + sp3_angle) % 360.0,
                    (incoming - sp3_angle) % 360.0,
                ]
                if self.state.fixed_angles and not sp3_exact_mode:
                    preferred = self._snap_angles_to_grid(preferred)

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
        """Método auxiliar para  select preferred angle.

        Args:
            anchor_id: Descripción del parámetro.
            cursor_angle: Descripción del parámetro.
            bond_order: Descripción del parámetro.
            is_aromatic: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  find overlapping bond.

        Args:
            p0: Descripción del parámetro.
            p1: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        for bond in self.model.bonds.values():
            a1 = self.model.get_atom(bond.a1_id)
            a2 = self.model.get_atom(bond.a2_id)
            b0 = QPointF(a1.x, a1.y)
            b1 = QPointF(a2.x, a2.y)
            if segments_nearly_equal(p0, p1, b0, b1, BOND_OVERLAP_TOLERANCE_PX):
                return bond.id
        return None

    def _snap_ring_vertex(self, pos: QPointF, excluded: set[int]) -> Optional[int]:
        """Método auxiliar para  snap ring vertex.

        Args:
            pos: Descripción del parámetro.
            excluded: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        candidates = [
            (atom_id, atom_x, atom_y)
            for atom_id, atom_x, atom_y in self._interactive_atom_candidates()
            if atom_id not in excluded
        ]
        return closest_atom(pos, candidates, ATOM_HIT_RADIUS)

    def _build_ring_vertex_defs(
        self, vertices: List[QPointF], anchor_type: str, anchor_id: Optional[int]
    ) -> List[Tuple[Optional[int], float, float]]:
        """Método auxiliar para  build ring vertex defs.

        Args:
            vertices: Descripción del parámetro.
            anchor_type: Descripción del parámetro.
            anchor_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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

    def _bond_drag_distance(self, scene_pos: QPointF) -> float:
        """Método auxiliar para  bond drag distance.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._bond_drag_start is None:
            return 0.0
        return (scene_pos - self._bond_drag_start).manhattanLength()

    def _update_bond_drag_state(self, scene_pos: QPointF) -> None:
        """Método auxiliar para  update bond drag state.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._drag_free_orientation = (
            self._bond_drag_start is not None
            and self._bond_drag_distance(scene_pos) >= BOND_DRAG_THRESHOLD_PX
        )

    def _should_use_default_bond_angle(
        self, modifiers: Qt.KeyboardModifiers, scene_pos: Optional[QPointF] = None
    ) -> bool:
        """Método auxiliar para  should use default bond angle.

        Args:
            modifiers: Descripción del parámetro.
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._drag_mode != "place_bond":
            return False
        if modifiers & Qt.KeyboardModifier.AltModifier:
            return False
        if scene_pos is None:
            scene_pos = self._last_scene_pos
        if self._bond_drag_distance(scene_pos) >= BOND_DRAG_THRESHOLD_PX:
            return False
        return True

    def _current_bond_geometry(self) -> tuple[int, bool]:
        """Método auxiliar para  current bond geometry.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        order = self.state.active_bond_order
        is_aromatic = self.state.active_bond_aromatic
        if self.state.active_bond_style != BondStyle.PLAIN or is_aromatic:
            order = 1
        return order, is_aromatic

    def _peek_default_bond_angle(self, anchor_id: Optional[int]) -> float:
        """Método auxiliar para  peek default bond angle.

        Args:
            anchor_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if anchor_id is None:
            return 0.0
        existing = self._get_anchor_bond_angles(anchor_id)
        step = self._bond_environment_step(anchor_id)
        if self._bond_last_angle is not None and existing:
            tolerance = math.radians(BOND_LAST_ANGLE_TOLERANCE_DEG)
            if min(self._angle_distance(self._bond_last_angle, a) for a in existing) <= tolerance:
                angle = self._bond_last_angle + step * self._bond_zigzag_sign
                return self._normalize_angle(angle)
        if not existing:
            return 0.0
        candidates = []
        for base in existing:
            candidates.append(base + step)
            candidates.append(base - step)
        return self._best_separated_angle(candidates, existing)

    def _compute_default_bond_endpoint(
        self, anchor: QPointF, anchor_id: Optional[int]
    ) -> QPointF:
        """Método auxiliar para  compute default bond endpoint.

        Args:
            anchor: Descripción del parámetro.
            anchor_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        order, aromatic = self._current_bond_geometry()
        default_angle_deg = math.degrees(self._peek_default_bond_angle(anchor_id))
        length = self.state.bond_length
        theta, final_length = self._pick_bond_direction_deg(
            anchor,
            anchor_id,
            default_angle_deg,
            order,
            aromatic,
            length,
            apply_collisions=True,
            allow_length_boost=False,
        )
        p1 = endpoint_from_angle_len(anchor, theta, final_length)
        snap_id = closest_atom(
            p1,
            self._interactive_atom_candidates(),
            ATOM_HIT_RADIUS,
        )
        if snap_id is not None and (anchor_id is None or snap_id != anchor_id):
            # Evita auto-snap al átomo que ya está enlazado con el ancla:
            # en sp3 congestionado (3 enlaces), ese snap fuerza superposición
            # y termina ciclando el enlace existente en vez de crear el 4º.
            if anchor_id is not None and self.model.find_bond_between(anchor_id, snap_id) is not None:
                return p1
            atom = self.model.get_atom(snap_id)
            return QPointF(atom.x, atom.y)
        return p1

    def _update_bond_angle_state(
        self,
        p0: QPointF,
        p1: QPointF,
        used_default: bool,
        anchor_id: Optional[int],
    ) -> None:
        """Método auxiliar para  update bond angle state.

        Args:
            p0: Descripción del parámetro.
            p1: Descripción del parámetro.
            used_default: Descripción del parámetro.
            anchor_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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

    def _begin_place_bond(self, anchor_id: Optional[int], scene_pos: QPointF) -> None:
        """Método auxiliar para  begin place bond.

        Args:
            anchor_id: Descripción del parámetro.
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Método auxiliar para  update bond preview.

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
        """Calcula curvaturas normalizadas para enlaces FLEX."""
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
        # Shift activa curva compleja en "S" (controles opuestos).
        if modifiers & Qt.KeyboardModifier.ShiftModifier:
            c2 = -c1
        else:
            c2 = c1
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
        """Activa ajuste interactivo de curvatura para enlace FLEX."""
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
        """Actualiza previsualización durante ajuste de enlace FLEX."""
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
        """Confirma creación del enlace FLEX con la curvatura previsualizada."""
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
        """Método auxiliar para  finalize bond.

        Args:
            modifiers: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
            self._begin_flex_adjust_mode(
                p0,
                final_p1,
                anchor_id,
                target_id,
                use_default,
                modifiers,
            )
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
        """Método auxiliar para  compute bond endpoint.

        Args:
            anchor: Descripción del parámetro.
            scene_pos: Descripción del parámetro.
            modifiers: Descripción del parámetro.
            bond_order: Descripción del parámetro.
            is_aromatic: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
            snap_id = closest_atom(
                p1,
                self._interactive_atom_candidates(),
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
            allow_length_boost=False,
        )
        p1 = endpoint_from_angle_len(anchor, theta, final_length)
        snap_id = closest_atom(
            p1,
            self._interactive_atom_candidates(),
            ATOM_HIT_RADIUS,
        )
        if snap_id is not None and (self._drag_anchor is None or snap_id != self._drag_anchor["id"]):
            atom = self.model.get_atom(snap_id)
            return QPointF(atom.x, atom.y)
        return p1

    def _set_bond_order_hotkey(self, bond_id: int, order: int) -> None:
        """Método auxiliar para  set bond order hotkey.

        Args:
            bond_id: Descripción del parámetro.
            order: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
