from __future__ import annotations

from ._shared import (
    ATOM_HIT_RADIUS,
    AddArrowCommand,
    AddAtomCommand,
    AddBracketCommand,
    ArrowItem,
    AtomItem,
    BondItem,
    BondStereo,
    BondStyle,
    BracketItem,
    ChangeAtomCommand,
    ChangeBondLengthCommand,
    ChangeChargeCommand,
    CompositeDiagramItem,
    EnergyDiagramItem,
    GelBandItem,
    GelElectrophoresisItem,
    ImageAnnotationItem,
    MoveArrowItemsCommand,
    MoveAtomsCommand,
    MoveBracketItemsCommand,
    MovePlateItemsCommand,
    MoveTextItemsCommand,
    Optional,
    OrbitalAnnotationItem,
    QGraphicsItem,
    QGraphicsTextItem,
    QGraphicsView,
    QInputDialog,
    QPointF,
    QRectF,
    Qt,
    SYMBOL_TEXT_TOOLS,
    SetCoordinationCenterCommand,
    TLCPlateItem,
    TLCSpotItem,
    TextAnnotationItem,
    TransformEnergyDiagramItemsCommand,
    TransformImageItemsCommand,
    TransformOrbitalItemsCommand,
    WavyAnchorItem,
    closest_atom,
    energy_diagram_kind_from_tool_id,
    orbital_kind_from_tool_id,
)

class CanvasSelectionInputMixin:
    def set_current_tool(self, tool_id: str) -> None:
        """Actualiza el estado de current tool.

        Args:
            tool_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        previous_tool = self.current_tool
        tool_id = tool_id or "tool_none"
        if tool_id.startswith("atom_"):
            self.state.default_element = tool_id.split("_", 1)[1]
            tool_id = "tool_atom"
        elif tool_id.startswith("coord_"):
            tool_id = "tool_coordination_center"
        elif tool_id.startswith("bond_"):
            tool_id = "tool_bond"
        elif tool_id == "tool_coordination_bond":
            self.state.active_bond_order = 1
            self.state.active_bond_style = BondStyle.COORDINATION
            self.state.active_bond_stereo = BondStereo.NONE
            self.state.active_bond_mode = "set"
            self.state.active_bond_aromatic = False
            tool_id = "tool_bond"
        elif tool_id.startswith("tool_brackets_"):
            bracket_key = tool_id.split("tool_brackets_", 1)[1]
            mapping = {
                "round": "()",
                "square": "[]",
                "square_left": "[",
                "square_right": "]",
                "corner": "corner",
                "curly": "{}",
                "curly_left": "{",
                "curly_right": "}",
                "frame": "frame",
                "frame_rounded": "rounded_frame",
            }
            self.state.active_bracket_type = mapping.get(bracket_key, "[]")
            tool_id = "tool_brackets"
        elif tool_id.startswith("tool_orbital_"):
            orbital_kind = orbital_kind_from_tool_id(tool_id)
            if orbital_kind:
                self.state.active_orbital_kind = orbital_kind
            tool_id = "tool_orbital"
        elif tool_id.startswith("tool_energy_diagram_"):
            diagram_kind = energy_diagram_kind_from_tool_id(tool_id)
            if diagram_kind:
                self.state.active_energy_diagram_kind = diagram_kind
            tool_id = "tool_energy_diagram"

        self.current_tool = tool_id
        self.state.active_tool = tool_id
        if self._pending_template_graph is not None:
            self.cancel_template_insert_mode()
        self._clear_bond_anchor()
        self._cancel_drag()
        if tool_id != "tool_text":
            self._finish_active_text_editing()
        selection_tools = {"tool_select", "tool_select_lasso", "tool_rotate_3d_precise"}
        if tool_id not in selection_tools and previous_tool in selection_tools:
            self.scene.clearSelection()
            self._sync_selection_from_scene()
        if previous_tool == "tool_text" and tool_id != "tool_text":
            self.scene.clearSelection()
            self._sync_selection_from_scene()
        self._arrow_start_pos = None
        self._arrow_end_pos = None
        self._arrow_curve_adjust_mode = False
        self._arrow_curve_factor = ArrowItem.CURVE_FACTOR_DEFAULT
        if hasattr(self, "_preview_arrow_item") and self._preview_arrow_item is not None:
            self._preview_arrow_item.hide_preview()
        self._clear_bracket_preview()

        if tool_id in {"tool_select", "tool_select_lasso", "tool_rotate_3d_precise", "tool_none"}:
            self.setCursor(Qt.CursorShape.ArrowCursor)
            self.setDragMode(QGraphicsView.DragMode.NoDrag)
        else:
            self.setCursor(Qt.CursorShape.CrossCursor)
            self.setDragMode(QGraphicsView.DragMode.NoDrag)

    def set_active_bond(self, bond_spec: dict) -> None:
        """Actualiza el estado de active bond.

        Args:
            bond_spec: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
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
        """Actualiza el estado de active ring.

        Args:
            ring_spec: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self.state.active_ring_size = ring_spec.get("size", 6)
        self.state.active_ring_aromatic = ring_spec.get("aromatic", False)
        self.state.active_ring_template = ring_spec.get("template")
        self.state.active_ring_anomeric = ring_spec.get("anomeric")
        self.set_current_tool("tool_ring")

    def set_active_element(self, element: str) -> None:
        """Actualiza el estado de active element.

        Args:
            element: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self.state.default_element = element
        self.set_current_tool("tool_atom")

    def mousePressEvent(self, event) -> None:
        """Método auxiliar para mousePressEvent.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if event.button() == Qt.MouseButton.MiddleButton or (
            self._space_panning and event.button() == Qt.MouseButton.LeftButton
        ):
            self._panning = True
            self._pan_last_pos = event.pos()
            self.setCursor(Qt.CursorShape.ClosedHandCursor)
            event.accept()
            return
        scene_pos = self.mapToScene(event.pos())
        self._last_scene_pos = scene_pos
        clicked_item = self._selected_annotation_item_at(scene_pos) or self._resolve_click_item(scene_pos)
        active_text_item = self._active_text_edit_item()
        if event.button() == Qt.MouseButton.LeftButton and active_text_item is not None:
            if clicked_item is active_text_item:
                super().mousePressEvent(event)
                return
            self._finish_active_text_editing(clear_scene_selection=True)
            if self.current_tool == "tool_text":
                self._accept_input_event(event)
                return

        if self._drag_mode == "flex_adjust":
            if event.button() == Qt.MouseButton.LeftButton:
                self._commit_flex_bond()
                return
            if event.button() == Qt.MouseButton.RightButton:
                self._cancel_drag()
                return

        if (
            not self._is_on_paper(scene_pos.x(), scene_pos.y())
            and self.current_tool not in {"tool_select", "tool_select_lasso", "tool_rotate_3d_precise"}
        ):
            super().mousePressEvent(event)
            return

        if event.button() == Qt.MouseButton.RightButton:
            if isinstance(clicked_item, (TLCSpotItem, GelBandItem)):
                super().mousePressEvent(event)
                return
            if isinstance(
                clicked_item,
                (
                    TLCPlateItem,
                    GelElectrophoresisItem,
                ),
            ):
                if not (event.modifiers() & (Qt.KeyboardModifier.ShiftModifier | Qt.KeyboardModifier.ControlModifier)):
                    self.scene.clearSelection()
                clicked_item.setSelected(True)
                self._sync_selection_from_scene()
                self._show_context_menu(
                    scene_pos,
                    event.globalPosition().toPoint(),
                    clicked_item=clicked_item,
                )
                return
            if isinstance(
                clicked_item,
                (
                    AtomItem,
                    BondItem,
                    ArrowItem,
                    BracketItem,
                    TextAnnotationItem,
                    EnergyDiagramItem,
                    CompositeDiagramItem,
                    OrbitalAnnotationItem,
                    ImageAnnotationItem,
                ),
            ):
                if not (event.modifiers() & (Qt.KeyboardModifier.ShiftModifier | Qt.KeyboardModifier.ControlModifier)):
                    self.scene.clearSelection()
                clicked_item.setSelected(True)
                self._sync_selection_from_scene()
            self._show_context_menu(
                scene_pos,
                event.globalPosition().toPoint(),
                clicked_item=clicked_item,
            )
            return

        if event.button() != Qt.MouseButton.LeftButton:
            super().mousePressEvent(event)
            return

        if self._maybe_begin_selected_annotation_transform(scene_pos, event):
            return

        if self.current_tool == "tool_none":
            return

        if self._pending_template_graph is not None:
            clicked_atom_id, _clicked_bond_id = self._pick_hover_target(scene_pos)
            if not self._is_on_paper(scene_pos.x(), scene_pos.y()):
                return
            self._insert_molgraph_at(
                self._pending_template_graph,
                scene_pos,
                attach_to_atom_id=clicked_atom_id,
            )
            self.cancel_template_insert_mode()
            return

        if self.current_tool in {"tool_select", "tool_select_lasso", "tool_rotate_3d_precise"}:
            if self._space_panning:
                return
            precise_mode = self.current_tool == "tool_rotate_3d_precise"
            if not precise_mode:
                blocked_modifiers = (
                    Qt.KeyboardModifier.ShiftModifier
                    | Qt.KeyboardModifier.ControlModifier
                    | Qt.KeyboardModifier.MetaModifier
                )
                if (
                    event.modifiers() & Qt.KeyboardModifier.AltModifier
                    and not (event.modifiers() & blocked_modifiers)
                    and self._begin_3d_rotation_drag(
                        scene_pos,
                        clicked_item.atom_id if isinstance(clicked_item, AtomItem) else None,
                    )
                ):
                    return
                if (
                    (event.modifiers() & Qt.KeyboardModifier.AltModifier)
                    and isinstance(clicked_item, AtomItem)
                    and self._cycle_anchor_override(clicked_item.atom_id)
                ):
                    return
                if event.button() == Qt.MouseButton.LeftButton:
                    handle_hit = self._selection_handle_hit_kind(scene_pos)
                    if handle_hit == "scale":
                        self._begin_scale_drag(scene_pos)
                        self._accept_input_event(event)
                        return
                    if handle_hit == "rotate":
                        self._begin_rotation_drag(scene_pos)
                        self._accept_input_event(event)
                        return
                    if handle_hit == "move":
                        if self._selection_move_handle is not None:
                            self._selection_move_handle.setCursor(Qt.CursorShape.ClosedHandCursor)
                        self._begin_drag(scene_pos)
                        self._accept_input_event(event)
                        return
            if clicked_item is None:
                additive = bool(
                    event.modifiers()
                    & (Qt.KeyboardModifier.ShiftModifier | Qt.KeyboardModifier.ControlModifier)
                )
                free_select = precise_mode or self.current_tool == "tool_select_lasso" or bool(
                    event.modifiers() & Qt.KeyboardModifier.ControlModifier
                )
                self._begin_selection_drag(scene_pos, free_select, additive)
                self._accept_input_event(event)
                return

            clicked_atom_id = None
            clicked_bond_id = None
            if precise_mode:
                clicked_atom_id, clicked_bond_id = self._pick_hover_target(scene_pos)

            if event.modifiers() & (Qt.KeyboardModifier.ShiftModifier | Qt.KeyboardModifier.ControlModifier):
                clicked_item.setSelected(not clicked_item.isSelected())
                self._sync_selection_from_scene()
                return
            else:
                if not clicked_item.isSelected():
                    self.scene.clearSelection()
                    clicked_item.setSelected(True)
            self._sync_selection_from_scene()
            if isinstance(clicked_item, (TLCSpotItem, GelBandItem)):
                super().mousePressEvent(event)
                return
            if precise_mode:
                self._prompt_precise_3d_rotation(
                    clicked_atom_id=clicked_atom_id,
                    clicked_bond_id=clicked_bond_id,
                )
                return
            if clicked_item.isSelected():
                if isinstance(clicked_item, TextAnnotationItem) and (
                    clicked_item.textInteractionFlags()
                    & Qt.TextInteractionFlag.TextEditorInteraction
                ):
                    return
                if isinstance(clicked_item, (TLCSpotItem, GelBandItem)):
                    super().mousePressEvent(event)
                    return
                if isinstance(
                    clicked_item,
                    (
                        AtomItem,
                        TextAnnotationItem,
                        ArrowItem,
                        BracketItem,
                        EnergyDiagramItem,
                        CompositeDiagramItem,
                        OrbitalAnnotationItem,
                        ImageAnnotationItem,
                        TLCPlateItem,
                        GelElectrophoresisItem,
                    ),
                ):
                    self._begin_drag(scene_pos)
                    self._accept_input_event(event)
            return

        clicked_atom_id, clicked_bond_id = self._pick_hover_target(scene_pos)

        arrow_tools = {
            "tool_arrow_line": "line",
            "tool_arrow_line_dashed": "line_dashed",
            "tool_arrow_forward": "forward",
            "tool_arrow_forward_open": "forward_open",
            "tool_arrow_forward_dashed": "forward_dashed",
            "tool_arrow_retro": "retro",
            "tool_arrow_retro_open": "retro_open",
            "tool_arrow_retro_dashed": "retro_dashed",
            "tool_arrow_both": "both",
            "tool_arrow_both_open": "both_open",
            "tool_arrow_both_dashed": "both_dashed",
            "tool_arrow_equilibrium": "equilibrium",
            "tool_arrow_equilibrium_dashed": "equilibrium_dashed",
            "tool_arrow_retrosynthetic": "retrosynthetic",
            "tool_arrow_curved": "curved",
            "tool_arrow_curved_fishhook": "curved_fishhook",
        }
        if self.current_tool in arrow_tools:
            arrow_kind = arrow_tools[self.current_tool]
            curved_tools = {"tool_arrow_curved", "tool_arrow_curved_fishhook"}
            is_curved_tool = self.current_tool in curved_tools
            constrained_pos = (
                self._constrain_annotation_endpoint(self._arrow_start_pos, scene_pos, event.modifiers())
                if self._arrow_start_pos is not None
                else scene_pos
            )
            if self._arrow_start_pos is None:
                self._arrow_start_pos = scene_pos
                self._arrow_end_pos = None
                self._arrow_curve_adjust_mode = False
                self._arrow_curve_factor = ArrowItem.CURVE_FACTOR_DEFAULT
                return
            if is_curved_tool:
                if self._arrow_end_pos is None:
                    if (scene_pos - self._arrow_start_pos).manhattanLength() < 2.0:
                        self._arrow_start_pos = None
                        self._arrow_end_pos = None
                        self._arrow_curve_adjust_mode = False
                        self._arrow_curve_factor = ArrowItem.CURVE_FACTOR_DEFAULT
                        self._preview_arrow_item.hide_preview()
                        return
                    self._arrow_end_pos = scene_pos
                    self._arrow_curve_adjust_mode = True
                    self._preview_arrow_item.update_preview(
                        self._arrow_start_pos,
                        self._arrow_end_pos,
                        arrow_kind,
                        curve_factor=self._arrow_curve_factor,
                    )
                    return
                cmd = AddArrowCommand(
                    self,
                    self._arrow_start_pos,
                    self._arrow_end_pos,
                    arrow_kind,
                    curve_factor=self._arrow_curve_factor,
                )
                self.undo_stack.push(cmd)
                self._arrow_start_pos = None
                self._arrow_end_pos = None
                self._arrow_curve_adjust_mode = False
                self._arrow_curve_factor = ArrowItem.CURVE_FACTOR_DEFAULT
                self._preview_arrow_item.hide_preview()
                return
            if (constrained_pos - self._arrow_start_pos).manhattanLength() < 2.0:
                self._arrow_start_pos = None
                self._arrow_end_pos = None
                self._arrow_curve_adjust_mode = False
                self._arrow_curve_factor = ArrowItem.CURVE_FACTOR_DEFAULT
                self._preview_arrow_item.hide_preview()
                return
            cmd = AddArrowCommand(self, self._arrow_start_pos, constrained_pos, arrow_kind)
            self.undo_stack.push(cmd)
            self._arrow_start_pos = None
            self._arrow_end_pos = None
            self._arrow_curve_adjust_mode = False
            self._arrow_curve_factor = ArrowItem.CURVE_FACTOR_DEFAULT
            self._preview_arrow_item.hide_preview()
            return

        if self.current_tool == "tool_brackets":
            self._bracket_drag_start = scene_pos
            preview = self._ensure_bracket_preview()
            preview.setRect(QRectF(scene_pos, scene_pos))
            return

        if self.current_tool in {"tool_charge", "tool_charge_plus", "tool_charge_minus"}:
            if clicked_atom_id is None:
                item = self._get_item_at(scene_pos)
                if isinstance(item, AtomItem):
                    clicked_atom_id = item.atom_id
            if clicked_atom_id is None:
                if self.current_tool == "tool_charge":
                    self._insert_symbol_item("±", scene_pos, None, 1.0, False)
                return
            atom = self.model.get_atom(clicked_atom_id)
            if self.current_tool == "tool_charge_plus":
                new_charge = int(atom.charge) + 1
            elif self.current_tool == "tool_charge_minus":
                new_charge = int(atom.charge) - 1
            else:
                if atom.charge == 0:
                    new_charge = 1
                elif atom.charge > 0:
                    new_charge = -1
                else:
                    new_charge = 0
            cmd = ChangeChargeCommand(self.model, self, clicked_atom_id, new_charge)
            self.undo_stack.push(cmd)
            return

        if self.current_tool == "tool_symbol_wavy_anchor":
            anchor_atom_id = clicked_atom_id
            if anchor_atom_id is None:
                item = self._get_item_at(scene_pos)
                if isinstance(item, AtomItem):
                    anchor_atom_id = item.atom_id
            if anchor_atom_id is None:
                return
            self._insert_wavy_anchor(scene_pos, anchor_atom_id)
            return

        symbol_spec = SYMBOL_TEXT_TOOLS.get(self.current_tool)
        if symbol_spec is not None:
            anchor_atom_id = clicked_atom_id
            if anchor_atom_id is None:
                item = self._get_item_at(scene_pos)
                if isinstance(item, AtomItem):
                    anchor_atom_id = item.atom_id
            if anchor_atom_id is None and self.model.atoms:
                atoms = [(atom.id, atom.x, atom.y) for atom in self.model.atoms.values()]
                pick_radius = max(ATOM_HIT_RADIUS * 3, self.state.bond_length * 1.5)
                anchor_atom_id = closest_atom(scene_pos, atoms, pick_radius)
            if symbol_spec.get("electrons"):
                self._insert_electron_dots(
                    scene_pos,
                    anchor_atom_id,
                    symbol_spec["electrons"],
                    symbol_spec.get("scale", 1.0),
                    spread=None,
                    mode="lone_pair" if self.current_tool == "tool_symbol_lone_pair" else "default",
                )
                return
            self._insert_symbol_item(
                symbol_spec["text"],
                scene_pos,
                anchor_atom_id,
                symbol_spec.get("scale", 1.0),
                symbol_spec.get("anchor", True),
                symbol_spec.get("rotate", False),
            )
            return

        if self.current_tool == "tool_energy_diagram":
            self._insert_energy_diagram_item(scene_pos)
            return

        if self.current_tool == "tool_orbital":
            target_pos = scene_pos
            if clicked_atom_id is not None:
                atom = self.model.get_atom(clicked_atom_id)
                if atom is not None:
                    target_pos = QPointF(atom.x, atom.y)
            self._insert_orbital_item(target_pos)
            return

        if self.current_tool == "tool_tlc":
            self._insert_tlc_plate_item(scene_pos)
            return

        if self.current_tool == "tool_electrophoresis":
            self._insert_gel_electrophoresis_item(scene_pos)
            return

        if self.current_tool == "tool_atom":
            element = self.state.default_element
            implicit_anchor_id, implicit_angle = self._pick_implicit_h_overlay(scene_pos)
            if implicit_anchor_id is not None:
                if self._substitute_implicit_hydrogen(implicit_anchor_id, element, implicit_angle):
                    return
            if clicked_atom_id is None:
                is_explicit = None
                if element == "C" and not self.state.show_implicit_carbons:
                    is_explicit = True
                elif element == "H" and not self.state.show_implicit_hydrogens:
                    is_explicit = True
                cmd = AddAtomCommand(
                    self.model,
                    self,
                    element,
                    scene_pos.x(),
                    scene_pos.y(),
                    is_explicit=is_explicit,
                )
                self.undo_stack.push(cmd)
            else:
                cmd = ChangeAtomCommand(self.model, self, clicked_atom_id, element)
                self.undo_stack.push(cmd)
            return

        if self.current_tool == "tool_coordination_center":
            element = self.state.default_element or "Pd"
            if element == "C":
                element = "Pd"
            if clicked_atom_id is None:
                cmd = AddAtomCommand(
                    self.model,
                    self,
                    element,
                    scene_pos.x(),
                    scene_pos.y(),
                    is_explicit=True,
                    is_coordination_center=True,
                    sphere_radius=16.0,
                    sphere_filled=True,
                )
                self.undo_stack.push(cmd)
            else:
                atom = self.model.get_atom(clicked_atom_id)
                needs_center_flag = not getattr(atom, "is_coordination_center", False)
                if not needs_center_flag:
                    return
                self.undo_stack.push(
                    SetCoordinationCenterCommand(
                        self.model,
                        self,
                        clicked_atom_id,
                        True,
                    )
                )
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
            if self.state.active_ring_template:
                self._insert_ring_template(scene_pos)
                return
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
            if clicked_bond_id is None:
                self._begin_place_chain(None, scene_pos)
                return

        if self.current_tool == "tool_erase":
            if clicked_atom_id is not None:
                self._delete_selection({clicked_atom_id}, set())
            elif clicked_bond_id is not None:
                self._delete_selection(set(), {clicked_bond_id})
            return

        elif self.current_tool == "tool_text":
            # Create text annotation at position
            item = TextAnnotationItem(
                "",
                scene_pos.x(),
                scene_pos.y()
            )
            self._apply_text_settings(item)
            # Enter edit mode immediately
            item.setTextInteractionFlags(Qt.TextInteractionFlag.TextEditorInteraction)
            self.scene.addItem(item)
            item.setFocus()
            # Position cursor at end to ensure it's visible
            cursor = item.textCursor()
            cursor.movePosition(cursor.MoveOperation.End)
            item.setTextCursor(cursor)
            
            self.scene.clearSelection()
            item.setSelected(True)
            return

        super().mousePressEvent(event)

    def mouseMoveEvent(self, event) -> None:
        """Método auxiliar para mouseMoveEvent.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._panning and self._pan_last_pos is not None:
            delta = event.pos() - self._pan_last_pos
            self._pan_last_pos = event.pos()
            hbar = self.horizontalScrollBar()
            vbar = self.verticalScrollBar()
            hbar.setValue(hbar.value() - delta.x())
            vbar.setValue(vbar.value() - delta.y())
            event.accept()
            return
        scene_pos = self.mapToScene(event.pos())
        self._last_scene_pos = scene_pos

        if self._is_rotating_3d:
            self._update_3d_rotation_drag(scene_pos)
            return

        if self._scale_dragging:
            self._update_scale_drag(scene_pos)
            self._accept_input_event(event)
            return
        if self._rotation_dragging:
            self._update_rotation_drag(scene_pos)
            self._accept_input_event(event)
            return

        if self._select_drag_mode is not None:
            self._update_selection_drag(scene_pos)
            return

        if self._drag_mode == "flex_adjust":
            self._update_flex_bond_preview(scene_pos, event.modifiers())
            return

        arrow_tools = {
            "tool_arrow_line": "line",
            "tool_arrow_line_dashed": "line_dashed",
            "tool_arrow_forward": "forward",
            "tool_arrow_forward_open": "forward_open",
            "tool_arrow_forward_dashed": "forward_dashed",
            "tool_arrow_retro": "retro",
            "tool_arrow_retro_open": "retro_open",
            "tool_arrow_retro_dashed": "retro_dashed",
            "tool_arrow_both": "both",
            "tool_arrow_both_open": "both_open",
            "tool_arrow_both_dashed": "both_dashed",
            "tool_arrow_equilibrium": "equilibrium",
            "tool_arrow_equilibrium_dashed": "equilibrium_dashed",
            "tool_arrow_retrosynthetic": "retrosynthetic",
            "tool_arrow_curved": "curved",
            "tool_arrow_curved_fishhook": "curved_fishhook",
        }
        if self.current_tool in arrow_tools and self._arrow_start_pos is not None:
            arrow_kind = arrow_tools[self.current_tool]
            curved_tools = {"tool_arrow_curved", "tool_arrow_curved_fishhook"}
            is_curved_tool = self.current_tool in curved_tools
            constrained_pos = self._constrain_annotation_endpoint(
                self._arrow_start_pos,
                scene_pos,
                event.modifiers(),
            )
            if is_curved_tool and self._arrow_end_pos is not None and self._arrow_curve_adjust_mode:
                curve_factor = self._curve_factor_from_pointer(
                    self._arrow_start_pos,
                    self._arrow_end_pos,
                    scene_pos,
                )
                if event.modifiers() & Qt.KeyboardModifier.ShiftModifier:
                    curve_factor = -curve_factor
                self._arrow_curve_factor = curve_factor
                self._preview_arrow_item.update_preview(
                    self._arrow_start_pos,
                    self._arrow_end_pos,
                    arrow_kind,
                    curve_factor=self._arrow_curve_factor,
                )
            else:
                preview_end = self._arrow_end_pos if (is_curved_tool and self._arrow_end_pos is not None) else constrained_pos
                preview_curve = self._arrow_curve_factor if is_curved_tool else None
                self._preview_arrow_item.update_preview(
                    self._arrow_start_pos,
                    preview_end,
                    arrow_kind,
                    curve_factor=preview_curve,
                )
            return

        if self._drag_mode == "place_bond":
            self._update_bond_drag_state(scene_pos)
            self._update_bond_preview(scene_pos, event.modifiers())
            return
        if self._drag_mode == "place_ring":
            self._update_ring_preview(scene_pos, event.modifiers())
            return
        if self._drag_mode == "place_chain":
            self._update_chain_preview(scene_pos, event.modifiers())
            return

        if self._bracket_drag_start is not None and self.current_tool == "tool_brackets":
            preview = self._ensure_bracket_preview()
            preview.setRect(QRectF(self._bracket_drag_start, scene_pos).normalized())
            return

        if self._dragging_selection and self._drag_start_pos is not None:
            delta = scene_pos - self._drag_start_pos
            if delta.manhattanLength() > 0:
                self._drag_has_moved = True

                # Move atoms
                for atom_id, (x, y) in self._drag_start_positions.items():
                    nx = x + delta.x()
                    ny = y + delta.y()
                    self.model.update_atom_position(atom_id, nx, ny)
                    self.update_atom_item(atom_id, nx, ny)
                self._update_live_drag_bond_geometry(set(self._drag_start_positions.keys()))
                
                # Move text items
                if hasattr(self, "_drag_start_text_positions"):
                    for item, (pos, rot) in self._drag_start_text_positions.items():
                        new_pos = pos + delta
                        item.setPos(new_pos)

                # Move arrows
                if hasattr(self, "_drag_start_arrow_positions"):
                    for item, (start, end) in self._drag_start_arrow_positions.items():
                        item.update_positions(start + delta, end + delta)

                # Move brackets
                if hasattr(self, "_drag_start_bracket_rects"):
                    for item, rect in self._drag_start_bracket_rects.items():
                        item.set_rect(rect.translated(delta))

                if hasattr(self, "_drag_start_energy_diagram_snapshots"):
                    for item, (pos, width, height, rotation) in self._drag_start_energy_diagram_snapshots.items():
                        item.set_display_rect(
                            QRectF(
                                pos.x() + delta.x(),
                                pos.y() + delta.y(),
                                width,
                                height,
                            )
                        )
                        item.setRotation(rotation)

                if hasattr(self, "_drag_start_semantic_diagram_snapshots"):
                    for item, (pos, width, height, rotation) in self._drag_start_semantic_diagram_snapshots.items():
                        item.set_display_rect(
                            QRectF(
                                pos.x() + delta.x(),
                                pos.y() + delta.y(),
                                width,
                                height,
                            )
                        )
                        item.setRotation(rotation)

                if hasattr(self, "_drag_start_orbital_snapshots"):
                    for item, (anchor0, anchor1) in self._drag_start_orbital_snapshots.items():
                        item.set_anchors(anchor0 + delta, anchor1 + delta)

                if hasattr(self, "_drag_start_image_snapshots"):
                    for item, (pos, width, height, rotation) in self._drag_start_image_snapshots.items():
                        item.set_display_rect(
                            QRectF(
                                pos.x() + delta.x(),
                                pos.y() + delta.y(),
                                width,
                                height,
                            )
                        )
                        item.setRotation(rotation)

                if hasattr(self, "_drag_start_plate_snapshots"):
                    for item, (pos, rotation) in self._drag_start_plate_snapshots.items():
                        item.setPos(pos + delta)
                        item.setRotation(rotation)

                self._update_drag_selection_overlay(delta)
            self._accept_input_event(event)
            return

        self._update_hover(scene_pos)

        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event) -> None:
        """Método auxiliar para mouseReleaseEvent.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._panning:
            self._panning = False
            self._pan_last_pos = None
            self.setCursor(Qt.CursorShape.ArrowCursor)
            event.accept()
            return
        if self._is_rotating_3d:
            self._finalize_3d_rotation_drag()
            return
        if self._scale_dragging:
            self._finalize_scale_drag()
            self._accept_input_event(event)
            return
        if self._rotation_dragging:
            self._finalize_rotation_drag()
            self._accept_input_event(event)
            return
        if self._select_drag_mode is not None:
            self._finalize_selection_drag()
            if self.current_tool == "tool_rotate_3d_precise":
                self._prompt_precise_3d_rotation()
            return
        if self._drag_mode == "flex_adjust":
            return
        if self._drag_mode == "place_bond":
            self._finalize_bond(event.modifiers())
            return
        if self._drag_mode == "place_ring":
            self._finalize_ring(event.modifiers())
            return
        if self._drag_mode == "place_chain":
            self._finalize_chain(event.modifiers())
            return

        if self._bracket_drag_start is not None and self.current_tool == "tool_brackets":
            rect = QRectF(self._bracket_drag_start, self._last_scene_pos).normalized()
            self._clear_bracket_preview()
            if rect.width() < 4.0 or rect.height() < 4.0:
                return
            kind = self.state.active_bracket_type
            pair = self._split_bracket_kind(kind)
            if pair:
                self.undo_stack.beginMacro("Add brackets")
                for side in pair:
                    self.undo_stack.push(AddBracketCommand(self, rect, side, padding=0.0))
                self.undo_stack.endMacro()
            else:
                cmd = AddBracketCommand(self, rect, kind, padding=0.0)
                self.undo_stack.push(cmd)
            return

        if self._dragging_selection:
            if self._selection_move_handle is not None:
                self._selection_move_handle.setCursor(Qt.CursorShape.OpenHandCursor)
            had_moved = bool(self._drag_has_moved)
            atom_before = dict(self._drag_start_positions)
            text_before = dict(getattr(self, "_drag_start_text_positions", {}))
            arrow_before = dict(getattr(self, "_drag_start_arrow_positions", {}))
            bracket_before = dict(getattr(self, "_drag_start_bracket_rects", {}))
            energy_before = dict(getattr(self, "_drag_start_energy_diagram_snapshots", {}))
            semantic_before = dict(getattr(self, "_drag_start_semantic_diagram_snapshots", {}))
            orbital_before = dict(getattr(self, "_drag_start_orbital_snapshots", {}))
            image_before = dict(getattr(self, "_drag_start_image_snapshots", {}))
            plate_before = dict(getattr(self, "_drag_start_plate_snapshots", {}))
            atom_after = (
                {
                    atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
                    for atom_id in atom_before
                    if atom_id in self.model.atoms
                }
                if had_moved and atom_before
                else {}
            )
            text_after = (
                {
                    item: (item.pos(), item.rotation())
                    for item in text_before.keys()
                }
                if had_moved and text_before
                else {}
            )
            arrow_after = (
                {
                    item: (item.start_point(), item.end_point())
                    for item in arrow_before.keys()
                }
                if had_moved and arrow_before
                else {}
            )
            bracket_after = (
                {
                    item: item.base_rect()
                    for item in bracket_before.keys()
                }
                if had_moved and bracket_before
                else {}
            )
            energy_after = (
                {
                    item: self._energy_diagram_transform_snapshot(item)
                    for item in energy_before.keys()
                }
                if had_moved and energy_before
                else {}
            )
            semantic_after = (
                {
                    item: self._semantic_diagram_transform_snapshot(item)
                    for item in semantic_before.keys()
                }
                if had_moved and semantic_before
                else {}
            )
            orbital_after = (
                {
                    item: self._orbital_transform_snapshot(item)
                    for item in orbital_before.keys()
                }
                if had_moved and orbital_before
                else {}
            )
            image_after = (
                {
                    item: self._image_transform_snapshot(item)
                    for item in image_before.keys()
                }
                if had_moved and image_before
                else {}
            )
            plate_after = (
                {
                    item: (item.pos(), item.rotation())
                    for item in plate_before.keys()
                }
                if had_moved and plate_before
                else {}
            )
            moved_atom_ids = set(atom_after.keys())

            self._dragging_selection = False
            self._drag_start_pos = None
            self._drag_start_positions = {}
            self._drag_start_text_positions = {}
            self._drag_start_arrow_positions = {}
            self._drag_start_bracket_rects = {}
            self._drag_start_energy_diagram_snapshots = {}
            self._drag_start_semantic_diagram_snapshots = {}
            self._drag_start_orbital_snapshots = {}
            self._drag_start_image_snapshots = {}
            self._drag_start_plate_snapshots = {}
            self._drag_start_selection_bbox = None
            self._drag_affected_bond_ids = set()
            self._drag_affects_ring_centers = False
            self._drag_has_moved = False

            if had_moved:
                if moved_atom_ids:
                    self.refresh_ring_centers()
                    self.update_bond_items_for_atoms(moved_atom_ids)
                self._update_selection_overlay()
                move_atoms = bool(atom_before and atom_after)
                move_text = bool(text_before and text_after)
                move_arrows = bool(arrow_before and arrow_after)
                move_brackets = bool(bracket_before and bracket_after)
                move_energy_diagrams = bool(energy_before and energy_after)
                move_semantic_diagrams = bool(semantic_before and semantic_after)
                move_orbitals = bool(orbital_before and orbital_after)
                move_images = bool(image_before and image_after)
                move_plates = bool(plate_before and plate_after)
                move_count = sum(
                    [
                        move_atoms,
                        move_text,
                        move_arrows,
                        move_brackets,
                        move_energy_diagrams,
                        move_semantic_diagrams,
                        move_orbitals,
                        move_images,
                        move_plates,
                    ]
                )

                if move_count > 1:
                    self.undo_stack.beginMacro("Move selection")

                if move_atoms:
                    self.undo_stack.push(
                        MoveAtomsCommand(
                            self.model,
                            self,
                            atom_before,
                            atom_after,
                            skip_first_redo=True,
                        )
                    )

                if move_text:
                    cmd = MoveTextItemsCommand(self, text_before, text_after)
                    # QUndoStack.push() llama redo(); en este caso la posición ya es la final.
                    self.undo_stack.push(cmd)

                if move_arrows:
                    self.undo_stack.push(MoveArrowItemsCommand(self, arrow_before, arrow_after))

                if move_brackets:
                    self.undo_stack.push(MoveBracketItemsCommand(self, bracket_before, bracket_after))

                if move_energy_diagrams:
                    self.undo_stack.push(
                        TransformEnergyDiagramItemsCommand(
                            self,
                            energy_before,
                            energy_after,
                            "Move energy diagrams",
                        )
                    )

                if move_semantic_diagrams:
                    self.undo_stack.push(
                        TransformImageItemsCommand(
                            self,
                            semantic_before,
                            semantic_after,
                            "Move semantic diagrams",
                        )
                    )

                if move_orbitals:
                    self.undo_stack.push(
                        TransformOrbitalItemsCommand(self, orbital_before, orbital_after, "Move orbitals")
                    )

                if move_images:
                    self.undo_stack.push(TransformImageItemsCommand(self, image_before, image_after, "Move images"))
                    
                if move_plates:
                    self.undo_stack.push(MovePlateItemsCommand(self, plate_before, plate_after))

                if move_count > 1:
                    self.undo_stack.endMacro()
            self._release_interaction_mouse()
            self._accept_input_event(event)
            return

        super().mouseReleaseEvent(event)

    def _get_item_at(self, scene_pos: QPointF):
        """Método auxiliar para  get item at.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        for item in self.scene.items(scene_pos):
            semantic_parent = self._semantic_diagram_parent(item)
            if semantic_parent is not None:
                return semantic_parent
            if isinstance(item, AtomItem) and self._is_disposable_orphan_atom(item.atom_id):
                continue
            if isinstance(
                item,
                (
                    AtomItem,
                    BondItem,
                    ArrowItem,
                    BracketItem,
                    TextAnnotationItem,
                    EnergyDiagramItem,
                    CompositeDiagramItem,
                    OrbitalAnnotationItem,
                    ImageAnnotationItem,
                    WavyAnchorItem,
                    TLCSpotItem,
                    GelBandItem,
                ),
            ):
                return item
            # If we clicked a label/text child of an atom, return the atom.
            if isinstance(item, QGraphicsTextItem):
                parent = item.parentItem()
                if isinstance(parent, AtomItem):
                    if self._is_disposable_orphan_atom(parent.atom_id):
                        continue
                    return parent
                semantic_parent = self._semantic_diagram_parent(parent)
                if semantic_parent is not None:
                    return semantic_parent
        return None

    @staticmethod
    def _semantic_diagram_parent(item: Optional[QGraphicsItem]) -> Optional[CompositeDiagramItem]:
        """Promueve hijos internos de un diagrama semántico a su item raíz."""
        current = item
        while current is not None:
            if isinstance(current, CompositeDiagramItem):
                return current
            current = current.parentItem()
        return None

    def _resolve_click_item(self, scene_pos: QPointF):
        """Resuelve el objetivo de clic priorizando el pick geométrico de átomos."""
        item = self._get_item_at(scene_pos)
        if isinstance(
            item,
            (
                TextAnnotationItem,
                EnergyDiagramItem,
                CompositeDiagramItem,
                OrbitalAnnotationItem,
                ImageAnnotationItem,
                ArrowItem,
                BracketItem,
                WavyAnchorItem,
            ),
        ):
            return item

        atom_id, bond_id = self._pick_hover_target(scene_pos)
        if atom_id is not None:
            atom_item = self.atom_items.get(atom_id)
            if atom_item is not None:
                return atom_item

        if item is not None:
            return item

        if bond_id is not None:
            return self.bond_items.get(bond_id)
        return None

    def _selected_annotation_item_at(self, scene_pos: QPointF) -> Optional[QGraphicsItem]:
        """Prioriza un orbital/imagen ya seleccionado si el puntero cae sobre él.

        Esto evita que objetos químicos superpuestos secuestren el gesto de
        mover/redimensionar/rotar cuando la intención del usuario es seguir
        manipulando una anotación activa.
        """
        best_item: Optional[QGraphicsItem] = None
        best_z = float("-inf")
        for item in [
            *self._selected_energy_diagram_items(),
            *self._selected_semantic_diagram_items(),
            *self._selected_orbital_items(),
            *self._selected_image_items(),
            *self._selected_plate_items(),
        ]:
            if item.scene() is not self.scene:
                continue
            try:
                if not item.isVisible():
                    continue
                if not item.sceneBoundingRect().contains(scene_pos):
                    continue
                local_pos = item.mapFromScene(scene_pos)
                if not item.contains(local_pos):
                    continue
                z_value = float(item.zValue())
            except RuntimeError:
                continue
            if best_item is None or z_value >= best_z:
                best_item = item
                best_z = z_value
        return best_item

    def _maybe_begin_selected_annotation_transform(self, scene_pos: QPointF, event) -> bool:
        """Permite transformar una anotación seleccionada incluso fuera de tool_select."""
        if (
            not self._selected_energy_diagram_items()
            and not self._selected_semantic_diagram_items()
            and not self._selected_orbital_items()
            and not self._selected_image_items()
            and not self._selected_plate_items()
        ):
            return False

        handle_hit = self._selection_handle_hit_kind(scene_pos)
        if handle_hit == "scale":
            self._begin_scale_drag(scene_pos)
            self._accept_input_event(event)
            return True
        if handle_hit == "rotate":
            self._begin_rotation_drag(scene_pos)
            self._accept_input_event(event)
            return True
        if handle_hit == "move":
            if self._selection_move_handle is not None:
                self._selection_move_handle.setCursor(Qt.CursorShape.ClosedHandCursor)
            self._begin_drag(scene_pos)
            self._accept_input_event(event)
            return True

        clicked_selected_annotation = self._selected_annotation_item_at(scene_pos)
        if clicked_selected_annotation is None:
            return False

        selection_tools = {"tool_select", "tool_select_lasso", "tool_rotate_3d_precise"}
        if self.current_tool not in selection_tools:
            self._begin_drag(scene_pos)
            self._accept_input_event(event)
            return True
        return False

    def _grab_interaction_mouse(self) -> None:
        """Fija el mouse a la vista durante transformaciones interactivas."""
        if self._interaction_mouse_grabbed:
            return
        viewport = self.viewport()
        try:
            viewport.grabMouse()
            self._interaction_mouse_grabbed = True
        except Exception:
            self._interaction_mouse_grabbed = False

    def _release_interaction_mouse(self) -> None:
        """Libera el mouse previamente fijado a la vista."""
        if not self._interaction_mouse_grabbed:
            return
        viewport = self.viewport()
        try:
            viewport.releaseMouse()
        except Exception:
            pass
        self._interaction_mouse_grabbed = False

    @staticmethod
    def _accept_input_event(event) -> None:
        """Acepta el evento solo si expone la API esperada por Qt."""
        accept = getattr(event, "accept", None)
        if callable(accept):
            accept()

    def mouseDoubleClickEvent(self, event) -> None:
        """Método auxiliar para mouseDoubleClickEvent.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        scene_pos = self.mapToScene(event.pos())
        item = self._resolve_click_item(scene_pos)
        if isinstance(item, BondItem):
            bond = self.model.get_bond(item.bond_id)
            current = bond.length_px if bond.length_px is not None else self.state.bond_length
            value, ok = QInputDialog.getDouble(
                self,
                "Longitud de enlace",
                "Longitud (px, 0 para usar global):",
                float(current),
                0.0,
                999.0,
                1,
            )
            if ok:
                new_length = None if value <= 0 else float(value)
                cmd = ChangeBondLengthCommand(self.model, self, bond.id, new_length)
                self.undo_stack.push(cmd)
            return
        if isinstance(item, AtomItem):
            atom = self.model.get_atom(item.atom_id)
            current_label = self._editable_atom_label(atom)
            value, anchor = self._prompt_atom_label(
                current_label, atom_id=atom.id, initial=current_label
            )
            if value is None:
                return
            cmd = ChangeAtomCommand(
                self.model,
                self,
                atom.id,
                value,
                anchor_override=anchor,
            )
            self.undo_stack.push(cmd)
            return
        super().mouseDoubleClickEvent(event)

    def _begin_drag(self, scene_pos: QPointF) -> None:
        """Método auxiliar para  begin drag.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if (
            not self._selected_atom_ids_for_transform()
            and not self._selected_text_items()
            and not self._selected_arrow_items()
            and not self._selected_bracket_items()
            and not self._selected_energy_diagram_items()
            and not self._selected_semantic_diagram_items()
            and not self._selected_orbital_items()
            and not self._selected_image_items()
            and not self._selected_plate_items()
        ):
            return
        self._dragging_selection = True
        self._grab_interaction_mouse()
        self._drag_start_pos = scene_pos
        
        # Capture atoms
        atom_ids = self._selected_atom_ids_for_transform()
        self._drag_start_positions = {
            atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
            for atom_id in atom_ids
            if atom_id in self.model.atoms
        }
        
        # Capture text items
        self._drag_start_text_positions = {}
        for item in self._selected_text_items():
            self._drag_start_text_positions[item] = (item.pos(), item.rotation())

        # Capture arrows
        self._drag_start_arrow_positions = {}
        for item in self._selected_arrow_items():
            self._drag_start_arrow_positions[item] = (item.start_point(), item.end_point())

        # Capture brackets
        self._drag_start_bracket_rects = {}
        for item in self._selected_bracket_items():
            self._drag_start_bracket_rects[item] = item.base_rect()

        self._drag_start_energy_diagram_snapshots = {}
        for item in self._selected_energy_diagram_items():
            self._drag_start_energy_diagram_snapshots[item] = self._energy_diagram_transform_snapshot(item)

        self._drag_start_semantic_diagram_snapshots = {}
        for item in self._selected_semantic_diagram_items():
            self._drag_start_semantic_diagram_snapshots[item] = self._semantic_diagram_transform_snapshot(item)

        self._drag_start_orbital_snapshots = {}
        for item in self._selected_orbital_items():
            self._drag_start_orbital_snapshots[item] = self._orbital_transform_snapshot(item)

        self._drag_start_image_snapshots = {}
        for item in self._selected_image_items():
            self._drag_start_image_snapshots[item] = self._image_transform_snapshot(item)
            
        self._drag_start_plate_snapshots = {}
        for item in self._selected_plate_items():
            self._drag_start_plate_snapshots[item] = (item.pos(), item.rotation())

        self._drag_start_selection_bbox = self._selected_items_bbox()
        self._drag_affected_bond_ids = {
            bond.id
            for bond in self.model.bonds.values()
            if bond.a1_id in atom_ids or bond.a2_id in atom_ids
        }
        self._drag_affects_ring_centers = any(
            self.model.bonds[bond_id].ring_id is not None
            for bond_id in self._drag_affected_bond_ids
            if bond_id in self.model.bonds
        )

        self._drag_has_moved = False

    def _update_live_drag_bond_geometry(self, atom_ids: set[int]) -> None:
        """Actualiza sólo la geometría indispensable durante un drag."""
        if not atom_ids or not self._drag_affected_bond_ids:
            return
        self._clear_trackball_reference_if_desynced(atom_ids)
        if self._drag_affects_ring_centers:
            self.refresh_ring_centers()

        for bond_id in self._drag_affected_bond_ids:
            bond = self.model.bonds.get(bond_id)
            item = self.bond_items.get(bond_id)
            if bond is None or item is None:
                continue
            atom1 = self.model.atoms.get(bond.a1_id)
            atom2 = self.model.atoms.get(bond.a2_id)
            if atom1 is None or atom2 is None:
                continue

            ring_center = self._aromatic_ring_center_for_bond(bond) if bond.is_aromatic else None
            if ring_center is None and bond.ring_id is not None:
                ring_center = self._ring_centers.get(bond.ring_id)

            item.set_ring_context(ring_center)
            item.set_offset_sign(self._bond_offset_sign(bond))
            item.update_positions(atom1, atom2)

        self._update_aromatic_circles_for_atoms(atom_ids)
