from __future__ import annotations

from ._shared import *

class CanvasInputMixin:
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

    @staticmethod
    def _bond_effective_order_for_render(bond: Bond) -> int:
        """Devuelve el orden efectivo usado para dibujar el enlace."""
        if bond.is_aromatic and bond.display_order is not None:
            return int(bond.display_order)
        return int(bond.order)

    def _selected_double_bond_for_pi_toggle(self) -> tuple[Bond, BondItem] | None:
        """Obtiene un único doble enlace seleccionado apto para invertir orientación."""
        if len(self.state.selected_bonds) != 1:
            return None
        if self.state.selected_atoms:
            return None
        bond_id = next(iter(self.state.selected_bonds))
        if bond_id not in self.model.bonds:
            return None
        bond = self.model.get_bond(bond_id)
        if self._bond_effective_order_for_render(bond) != 2:
            return None
        item = self.bond_items.get(bond_id)
        if item is None:
            return None
        return bond, item

    def toggle_selected_double_bond_orientation(self) -> bool:
        """Invierte orientación de un doble enlace seleccionado usando undo/redo."""
        selected = self._selected_double_bond_for_pi_toggle()
        if selected is None:
            return False
        bond, item = selected
        old_sign = getattr(bond, "pi_offset_sign", None)
        new_sign = item.next_manual_pi_offset()
        cmd = ChangeDoubleBondOrientationCommand(
            self.model,
            self,
            bond.id,
            old_sign=old_sign,
            new_sign=new_sign,
        )
        self.undo_stack.push(cmd)
        return True

    def wheelEvent(self, event: QWheelEvent) -> None:
        """Método auxiliar para wheelEvent.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        modifiers = event.modifiers()
        # Atajo para invertir línea pi en dobles enlaces:
        # Ctrl+Alt+Rueda evita conflicto con Shift+Rueda (scroll horizontal).
        if (
            (modifiers & Qt.KeyboardModifier.ControlModifier)
            and (modifiers & Qt.KeyboardModifier.AltModifier)
        ):
            if self.toggle_selected_double_bond_orientation():
                event.accept()
                return
        if modifiers & Qt.KeyboardModifier.ShiftModifier:
            hbar = self.horizontalScrollBar()
            delta = event.angleDelta().y()
            if delta == 0:
                delta = event.angleDelta().x()
            if delta != 0:
                hbar.setValue(hbar.value() - delta)
                event.accept()
                return
        if event.angleDelta().y() > 0:
            self.zoom_in()
        else:
            self.zoom_out()

    def keyPressEvent(self, event) -> None:
        """Método auxiliar para keyPressEvent.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        focus_item = self.scene.focusItem()
        if isinstance(focus_item, TextAnnotationItem) and (
            focus_item.textInteractionFlags()
            & Qt.TextInteractionFlag.TextEditorInteraction
        ):
            # Let the text editor handle Delete/Backspace while editing text.
            super().keyPressEvent(event)
            return
        if isinstance(focus_item, EnergyDiagramItem) and focus_item.is_editing():
            super().keyPressEvent(event)
            return
        if event.key() == Qt.Key.Key_Space and not self._space_panning:
            self._space_panning = True
            self.setCursor(Qt.CursorShape.OpenHandCursor)
            event.accept()
            return
        modifiers = event.modifiers()
        ctrl = modifiers == Qt.KeyboardModifier.ControlModifier
        if ctrl and event.key() == Qt.Key.Key_X:
            self.cut_to_clipboard()
            event.accept()
            return
        if ctrl and event.key() == Qt.Key.Key_D:
            self.duplicate_selection()
            event.accept()
            return
        if event.key() in (Qt.Key.Key_Delete, Qt.Key.Key_Backspace):
            has_selected_arrows = any(
                isinstance(item, ArrowItem) for item in self.scene.selectedItems()
            )
            has_selected_brackets = any(
                isinstance(item, BracketItem) for item in self.scene.selectedItems()
            )
            has_selected_text = any(
                isinstance(item, TextAnnotationItem) for item in self.scene.selectedItems()
            )
            has_selected_orbitals = any(
                isinstance(item, OrbitalAnnotationItem) for item in self.scene.selectedItems()
            )
            has_selected_energy_diagrams = any(
                isinstance(item, EnergyDiagramItem) for item in self.scene.selectedItems()
            )
            has_selected_semantic_diagrams = any(
                isinstance(item, CompositeDiagramItem) for item in self.scene.selectedItems()
            )
            has_selected_images = any(
                isinstance(item, ImageAnnotationItem) for item in self.scene.selectedItems()
            )
            has_selected_plates = any(
                isinstance(item, (TLCPlateItem, GelElectrophoresisItem))
                for item in self.scene.selectedItems()
            )
            has_selected_spots_bands = any(
                isinstance(item, (TLCSpotItem, GelBandItem))
                for item in self.scene.selectedItems()
            )
            if (
                self.state.selected_atoms
                or self.state.selected_bonds
                or has_selected_arrows
                or has_selected_brackets
                or has_selected_text
                or has_selected_energy_diagrams
                or has_selected_semantic_diagrams
                or has_selected_orbitals
                or has_selected_images
                or has_selected_plates
                or has_selected_spots_bands
            ):
                self.delete_selection()
                return
            if self._delete_hovered():
                return
        if event.key() == Qt.Key.Key_Escape and self._pending_template_graph is not None:
            self.cancel_template_insert_mode()
            event.accept()
            return
        if self._handle_atom_text_entry(event):
            return
        if self._handle_select_all(event):
            return
        if self._handle_nudge(event):
            return
        if self._handle_hotkeys(event):
            return
        super().keyPressEvent(event)

    def remember_text_edit_item(self, item: TextAnnotationItem | None) -> None:
        """Método auxiliar para remember text edit item.

        Args:
            item: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._last_text_edit_item = item

    def restore_text_edit_focus(self) -> None:
        """Método auxiliar para restore text edit focus.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = getattr(self, "_last_text_edit_item", None)
        if item is None:
            return
        if item.scene() is not self.scene:
            self._last_text_edit_item = None
            return
        item.setTextInteractionFlags(Qt.TextInteractionFlag.TextEditorInteraction)
        item.setFocus()

    def keyReleaseEvent(self, event) -> None:
        """Método auxiliar para keyReleaseEvent.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if event.key() == Qt.Key.Key_Space and self._space_panning:
            self._space_panning = False
            if not self._panning:
                self.setCursor(Qt.CursorShape.ArrowCursor)
            event.accept()
            return
        super().keyReleaseEvent(event)

    def _handle_atom_text_entry(self, event) -> bool:
        """Método auxiliar para  handle atom text entry.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if event.modifiers() & (
            Qt.KeyboardModifier.ControlModifier
            | Qt.KeyboardModifier.AltModifier
            | Qt.KeyboardModifier.MetaModifier
        ):
            return False
        if not self.state.selected_atoms or self.state.selected_bonds:
            return False
        if len(self.state.selected_atoms) != 1:
            return False
        typed = event.text()
        if not typed or not typed.isalpha():
            return False

        atom_id = next(iter(self.state.selected_atoms))
        atom = self.model.get_atom(atom_id)
        seed = self._normalize_atom_label(typed) or typed
        value, anchor = self._prompt_atom_label(
            self._editable_atom_label(atom), atom_id=atom_id, initial=seed
        )
        if value is None:
            return True
        cmd = ChangeAtomCommand(
            self.model,
            self,
            atom_id,
            value,
            anchor_override=anchor,
        )
        self.undo_stack.push(cmd)
        return True

    def _handle_nudge(self, event) -> bool:
        """Método auxiliar para  handle nudge.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        key = event.key()
        if key not in (
            Qt.Key.Key_Left,
            Qt.Key.Key_Right,
            Qt.Key.Key_Up,
            Qt.Key.Key_Down,
        ):
            return False

        focus_item = self.scene.focusItem()
        if isinstance(focus_item, TextAnnotationItem) and (
            focus_item.textInteractionFlags()
            & Qt.TextInteractionFlag.TextEditorInteraction
        ):
            return False

        selected_atom_ids = self._selected_atom_ids_for_transform()
        selected_text_items = self._selected_text_items()
        selected_arrows = self._selected_arrow_items()
        selected_brackets = self._selected_bracket_items()
        selected_energy_diagrams = self._selected_energy_diagram_items()
        selected_semantic_diagrams = self._selected_semantic_diagram_items()
        selected_orbitals = self._selected_orbital_items()
        selected_images = self._selected_image_items()
        if (
            not selected_atom_ids
            and not selected_text_items
            and not selected_arrows
            and not selected_brackets
            and not selected_energy_diagrams
            and not selected_semantic_diagrams
            and not selected_orbitals
            and not selected_images
        ):
            return False

        step = 1.0
        if event.modifiers() & Qt.KeyboardModifier.ShiftModifier:
            step = 10.0

        dx = 0.0
        dy = 0.0
        if key == Qt.Key.Key_Left:
            dx = -step
        elif key == Qt.Key.Key_Right:
            dx = step
        elif key == Qt.Key.Key_Up:
            dy = -step
        elif key == Qt.Key.Key_Down:
            dy = step

        if selected_atom_ids:
            before = {
                atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
                for atom_id in selected_atom_ids
                if atom_id in self.model.atoms
            }
            after = {
                atom_id: (x + dx, y + dy)
                for atom_id, (x, y) in before.items()
            }
            self.undo_stack.push(MoveAtomsCommand(self.model, self, before, after))

        if selected_text_items:
            before = {item: (item.pos(), item.rotation()) for item in selected_text_items}
            after = {
                item: (QPointF(pos.x() + dx, pos.y() + dy), rot)
                for item, (pos, rot) in before.items()
            }
            self.undo_stack.push(MoveTextItemsCommand(self, before, after))

        if selected_arrows:
            before = {item: (item.start_point(), item.end_point()) for item in selected_arrows}
            after = {
                item: (start + QPointF(dx, dy), end + QPointF(dx, dy))
                for item, (start, end) in before.items()
            }
            self.undo_stack.push(MoveArrowItemsCommand(self, before, after))

        if selected_brackets:
            before = {item: item.base_rect() for item in selected_brackets}
            after = {item: rect.translated(dx, dy) for item, rect in before.items()}
            self.undo_stack.push(MoveBracketItemsCommand(self, before, after))

        if selected_energy_diagrams:
            before = {item: self._energy_diagram_transform_snapshot(item) for item in selected_energy_diagrams}
            after = {
                item: (QPointF(pos.x() + dx, pos.y() + dy), width, height, rotation)
                for item, (pos, width, height, rotation) in before.items()
            }
            self.undo_stack.push(
                TransformEnergyDiagramItemsCommand(self, before, after, "Move energy diagrams")
            )

        if selected_semantic_diagrams:
            before = {item: self._semantic_diagram_transform_snapshot(item) for item in selected_semantic_diagrams}
            after = {
                item: (QPointF(pos.x() + dx, pos.y() + dy), width, height, rotation)
                for item, (pos, width, height, rotation) in before.items()
            }
            self.undo_stack.push(
                TransformImageItemsCommand(self, before, after, "Move semantic diagrams")
            )

        if selected_orbitals:
            before = {item: self._orbital_transform_snapshot(item) for item in selected_orbitals}
            after = {
                item: (anchor0 + QPointF(dx, dy), anchor1 + QPointF(dx, dy))
                for item, (anchor0, anchor1) in before.items()
            }
            self.undo_stack.push(TransformOrbitalItemsCommand(self, before, after, "Move orbitals"))

        if selected_images:
            before = {item: self._image_transform_snapshot(item) for item in selected_images}
            after = {
                item: (QPointF(pos.x() + dx, pos.y() + dy), width, height, rotation)
                for item, (pos, width, height, rotation) in before.items()
            }
            self.undo_stack.push(TransformImageItemsCommand(self, before, after, "Move images"))

        self._update_selection_overlay()
        return True

    def _select_all_items(self) -> None:
        """Método auxiliar para  select all items.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self.scene.clearSelection()
        for item in self.scene.items():
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
                ),
            ) and item.isVisible():
                item.setSelected(True)
        self._sync_selection_from_scene()

    def _handle_select_all(self, event) -> bool:
        """Método auxiliar para  handle select all.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if not (event.modifiers() & (Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.MetaModifier)):
            return False
        if event.key() not in (Qt.Key.Key_A, Qt.Key.Key_E):
            return False
        focus_item = self.scene.focusItem()
        if isinstance(focus_item, TextAnnotationItem) and (
            focus_item.textInteractionFlags()
            & Qt.TextInteractionFlag.TextEditorInteraction
        ):
            return False
        self._select_all_items()
        return True

    def _show_context_menu(
        self,
        scene_pos: QPointF,
        global_pos,
        clicked_item: Optional[QGraphicsItem] = None,
    ) -> None:
        """Método auxiliar para  show context menu.

        Args:
            scene_pos: Descripción del parámetro.
            global_pos: Descripción del parámetro.
            clicked_item: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        menu = QMenu(self)
        has_selection = bool(
            self.state.selected_atoms
            or self.state.selected_bonds
            or self._selected_text_items()
            or any(
                isinstance(
                    item,
                    (
                        ArrowItem,
                        BracketItem,
                        EnergyDiagramItem,
                        CompositeDiagramItem,
                        OrbitalAnnotationItem,
                        ImageAnnotationItem,
                        TLCPlateItem,
                        GelElectrophoresisItem,
                    ),
                )
                for item in self.scene.selectedItems()
            )
        )
        has_bond_selection = bool(self.state.selected_bonds)
        has_stroke_selection = bool(
            has_bond_selection or self._selected_arrow_items() or self._selected_bracket_items()
        )
        selected_line_arrow_items = self._selected_line_arrow_items()
        line_arrow_target = (
            clicked_item
            if self._is_uniform_line_arrow(clicked_item)
            else selected_line_arrow_items[0] if len(selected_line_arrow_items) == 1 else None
        )
        coordination_bond_ids: list[int] = []
        if isinstance(clicked_item, BondItem):
            coordination_bond_ids = [clicked_item.bond_id]
        elif has_bond_selection:
            coordination_bond_ids = sorted(
                bond_id for bond_id in self.state.selected_bonds if bond_id in self.model.bonds
            )
        coordination_bond_enable = False
        if coordination_bond_ids:
            coordination_bond_enable = not all(
                self.model.get_bond(bond_id).style == BondStyle.COORDINATION
                for bond_id in coordination_bond_ids
            )
        coordination_atom_ids: list[int] = []
        if isinstance(clicked_item, AtomItem):
            coordination_atom_ids = [clicked_item.atom_id]
        else:
            coordination_atom_ids = sorted(
                atom_id
                for atom_id in self._selected_atom_ids_for_transform()
                if atom_id in self.model.atoms
            )
        coordination_enable = False
        if coordination_atom_ids:
            coordination_enable = not all(
                getattr(self.model.get_atom(atom_id), "is_coordination_center", False)
                for atom_id in coordination_atom_ids
            )
        styled_coordination_ids = [
            atom_id
            for atom_id in coordination_atom_ids
            if bool(getattr(self.model.get_atom(atom_id), "is_coordination_center", False))
        ]
        coordination_all_filled = bool(styled_coordination_ids) and all(
            bool(getattr(self.model.get_atom(atom_id), "sphere_filled", True))
            for atom_id in styled_coordination_ids
        )
        coordination_all_transparent = bool(styled_coordination_ids) and all(
            bool(getattr(self.model.get_atom(atom_id), "sphere_transparent", False))
            for atom_id in styled_coordination_ids
        )
        coordination_has_custom_color = any(
            bool(getattr(self.model.get_atom(atom_id), "sphere_color", None))
            for atom_id in styled_coordination_ids
        )
        selected_energy_diagram_items = self._selected_energy_diagram_items()
        selected_row_energy_diagram_items = [
            item for item in selected_energy_diagram_items if item.family() == "row"
        ]
        energy_diagram_target = self._single_energy_diagram_target(clicked_item)
        energy_fill_all_visible = bool(selected_energy_diagram_items) and all(
            bool(item.effective_style().get("fill_visible", True))
            for item in selected_energy_diagram_items
        )
        energy_box_outline_all_visible = bool(selected_energy_diagram_items) and all(
            bool(item.effective_style().get("box_stroke_visible", True))
            for item in selected_energy_diagram_items
        )
        energy_box_base_all_visible = bool(selected_energy_diagram_items) and all(
            bool(item.effective_style().get("box_base_visible", False))
            for item in selected_energy_diagram_items
        )
        energy_has_custom_style = any(
            bool(item.style_payload()) for item in selected_energy_diagram_items
        )
        selected_semantic_diagram_items = self._selected_semantic_diagram_items()
        semantic_diagram_target = self._single_semantic_diagram_target(clicked_item)
        selected_orbital_items = self._selected_orbital_items()
        orbital_target = self._single_orbital_target(clicked_item)
        charge_atom_ids: list[int] = []
        if isinstance(clicked_item, AtomItem):
            charge_atom_ids = [clicked_item.atom_id]
        elif self.state.selected_atoms:
            charge_atom_ids = sorted(
                atom_id for atom_id in self.state.selected_atoms if atom_id in self.model.atoms
            )
        charge_all_no_implicit = bool(charge_atom_ids) and all(
            bool(getattr(self.model.get_atom(atom_id), "no_implicit", False))
            for atom_id in charge_atom_ids
        )

        act_cut = menu.addAction("Cortar")
        act_copy = menu.addAction("Copiar")
        act_paste = menu.addAction("Pegar")
        menu.addSeparator()
        act_delete = menu.addAction("Eliminar")
        act_scale_selection = menu.addAction("Redimensionar selección...") if has_selection else None
        act_thicker = None
        act_thinner = None
        act_reset_thickness = None
        act_line_equalize = None
        act_line_set_length = None
        act_color = None
        act_reset_color = None
        branch_action_handlers: dict = {}
        if has_stroke_selection:
            menu.addSeparator()
            thickness_menu = menu.addMenu("Grosor de enlace/flecha/corchete")
            act_thicker = thickness_menu.addAction("Incrementar grosor")
            act_thinner = thickness_menu.addAction("Disminuir grosor")
            act_reset_thickness = thickness_menu.addAction("Restablecer grosor")
        if selected_line_arrow_items:
            line_menu = menu.addMenu("Longitud de lineas")
            if len(selected_line_arrow_items) >= 2:
                equalize_label = (
                    "Igualar a la linea activa"
                    if line_arrow_target is not None
                    else "Igualar longitudes"
                )
                act_line_equalize = line_menu.addAction(equalize_label)
            act_line_set_length = line_menu.addAction("Establecer longitud...")
        branch_bond_id = None
        if isinstance(clicked_item, BondItem):
            branch_bond_id = clicked_item.bond_id
        else:
            branch_bond_id = self._selected_single_bond_id()
        if branch_bond_id is not None and branch_bond_id in self.model.bonds:
            menu.addSeparator()
            branch_menu = menu.addMenu("Reordenar rama")
            bond = self.model.get_bond(branch_bond_id)
            branch_contexts = [
                self._branch_rotation_context(branch_bond_id, fixed_atom_id=bond.a1_id),
                self._branch_rotation_context(branch_bond_id, fixed_atom_id=bond.a2_id),
            ]
            branch_contexts = [context for context in branch_contexts if context is not None]
            if branch_contexts:
                for context in branch_contexts:
                    side_menu = branch_menu.addMenu(self._branch_context_title(context))
                    act_branch_minus = side_menu.addAction(
                        f"Girar -{int(BRANCH_ROTATION_STEP_DEG)}°"
                    )
                    act_branch_plus = side_menu.addAction(
                        f"Girar +{int(BRANCH_ROTATION_STEP_DEG)}°"
                    )
                    act_branch_flip = side_menu.addAction("Invertir rama (180°)")
                    side_menu.addSeparator()
                    act_branch_auto = side_menu.addAction("Autoacomodar rama")
                    branch_action_handlers[act_branch_minus] = (
                        lambda ctx=context: self.rotate_branch_side_degrees(
                            int(ctx["bond_id"]),
                            int(ctx["fixed_atom_id"]),
                            -BRANCH_ROTATION_STEP_DEG,
                        )
                    )
                    branch_action_handlers[act_branch_plus] = (
                        lambda ctx=context: self.rotate_branch_side_degrees(
                            int(ctx["bond_id"]),
                            int(ctx["fixed_atom_id"]),
                            BRANCH_ROTATION_STEP_DEG,
                        )
                    )
                    branch_action_handlers[act_branch_flip] = (
                        lambda ctx=context: self.invert_branch_side(
                            int(ctx["bond_id"]),
                            int(ctx["fixed_atom_id"]),
                        )
                    )
                    branch_action_handlers[act_branch_auto] = (
                        lambda ctx=context: self.auto_arrange_branch_side(
                            int(ctx["bond_id"]),
                            int(ctx["fixed_atom_id"]),
                        )
                    )
            else:
                act_branch_unavailable = branch_menu.addAction(
                    "No disponible en enlaces cíclicos o coordinativos"
                )
                act_branch_unavailable.setEnabled(False)
        if has_bond_selection:
            act_color = menu.addAction("Color de enlace...")
            act_reset_color = menu.addAction("Restablecer color de enlace")
        act_anchor = None
        act_charge_increment = None
        act_charge_decrement = None
        act_charge_set = None
        act_charge_reset = None
        act_no_implicit = None
        act_toggle_coord_sphere = None
        act_toggle_coord_bond = None
        act_coord_sphere_color = None
        act_coord_sphere_reset_color = None
        act_coord_sphere_toggle_fill = None
        act_coord_sphere_toggle_transparent = None
        act_coord_sphere_radius = None
        act_arrange_coordination = None
        act_energy_label = None
        act_energy_box_count = None
        act_energy_occupancies = None
        act_energy_label_left = None
        act_energy_label_right = None
        act_energy_stroke_color = None
        act_energy_fill_color = None
        act_energy_label_color = None
        act_energy_arrow_up_color = None
        act_energy_arrow_down_color = None
        act_energy_fill_toggle = None
        act_energy_box_outline_toggle = None
        act_energy_box_base_toggle = None
        act_energy_reset_style = None
        act_semantic_title = None
        act_semantic_edit_builder = None
        act_semantic_summary_toggle = None
        act_semantic_level_label = None
        act_semantic_lane_label = None
        act_orbital_color = None
        act_orbital_reset_color = None
        act_orbital_opacity = None
        act_orbital_lobe_color = None
        act_orbital_lobe_opacity = None
        act_orbital_lobe_move = None
        act_orbital_lobe_reset = None
        if isinstance(clicked_item, AtomItem):
            atom = self.model.get_atom(clicked_item.atom_id)
            if atom is not None and atom.element not in ELEMENT_SYMBOLS:
                candidates = self._label_anchor_candidates(atom.element)
                if candidates:
                    menu.addSeparator()
                    act_anchor = menu.addAction("Elegir átomo de unión...")
        if charge_atom_ids:
            menu.addSeparator()
            charge_menu = menu.addMenu("Carga formal")
            act_charge_increment = charge_menu.addAction("Incrementar (+1)")
            act_charge_decrement = charge_menu.addAction("Disminuir (-1)")
            act_charge_set = charge_menu.addAction("Establecer...")
            act_charge_reset = charge_menu.addAction("Restablecer a 0")
            charge_menu.addSeparator()
            act_no_implicit = charge_menu.addAction("No implícitos")
            act_no_implicit.setCheckable(True)
            act_no_implicit.setChecked(charge_all_no_implicit)
        if coordination_bond_ids:
            menu.addSeparator()
            bond_text = (
                "Convertir a enlace coordinativo"
                if coordination_bond_enable
                else "Convertir a enlace normal"
            )
            act_toggle_coord_bond = menu.addAction(bond_text)
        if coordination_atom_ids:
            menu.addSeparator()
            coord_text = (
                "Activar esfera de coordinación"
                if coordination_enable
                else "Desactivar esfera de coordinación"
            )
            act_toggle_coord_sphere = menu.addAction(coord_text)
        if styled_coordination_ids:
            sphere_style_menu = menu.addMenu("Estilo de esfera")
            transparent_text = (
                "Mostrar esfera"
                if coordination_all_transparent
                else "Hacer esfera transparente"
            )
            act_coord_sphere_toggle_transparent = sphere_style_menu.addAction(transparent_text)
            act_coord_sphere_color = sphere_style_menu.addAction("Color de esfera...")
            act_coord_sphere_reset_color = sphere_style_menu.addAction("Restablecer color")
            fill_text = (
                "Quitar fondo de esfera"
                if coordination_all_filled
                else "Mostrar fondo de esfera"
            )
            act_coord_sphere_toggle_fill = sphere_style_menu.addAction(fill_text)
            act_coord_sphere_radius = sphere_style_menu.addAction("Tamaño de esfera...")
            act_coord_sphere_reset_color.setEnabled(coordination_has_custom_color)
            act_coord_sphere_color.setEnabled(not coordination_all_transparent)
            act_coord_sphere_reset_color.setEnabled(
                coordination_has_custom_color and not coordination_all_transparent
            )
            act_coord_sphere_toggle_fill.setEnabled(not coordination_all_transparent)
        arrangement_center_ids = [
            atom_id
            for atom_id in coordination_atom_ids
            if bool(getattr(self.model.get_atom(atom_id), "is_coordination_center", False))
        ]
        if arrangement_center_ids:
            can_arrange = any(
                len(self._coordination_ligands_for_center(center_id)) >= 2
                for center_id in arrangement_center_ids
            )
            act_arrange_coordination = menu.addAction("Distribuir ligandos alrededor")
            act_arrange_coordination.setEnabled(can_arrange)
        if selected_energy_diagram_items:
            menu.addSeparator()
            energy_menu = menu.addMenu("Diagrama de energia")
            if energy_diagram_target is not None and energy_diagram_target.supports_free_label():
                act_energy_label = energy_menu.addAction("Etiqueta...")
                side_menu = energy_menu.addMenu("Posición de etiqueta")
                act_energy_label_left = side_menu.addAction("Izquierda")
                act_energy_label_right = side_menu.addAction("Derecha")
                energy_menu.addSeparator()
            if selected_row_energy_diagram_items:
                act_energy_box_count = energy_menu.addAction("Numero de cajas...")
            if energy_diagram_target is not None:
                act_energy_occupancies = energy_menu.addAction("Ocupación electrónica...")
                energy_menu.addSeparator()
            energy_style_menu = energy_menu.addMenu("Estilo")
            act_energy_stroke_color = energy_style_menu.addAction("Color de contorno...")
            act_energy_fill_color = energy_style_menu.addAction("Color de relleno...")
            act_energy_label_color = energy_style_menu.addAction("Color de etiqueta...")
            act_energy_arrow_up_color = energy_style_menu.addAction("Color de flecha ↑...")
            act_energy_arrow_down_color = energy_style_menu.addAction("Color de flecha ↓...")
            fill_text = "Quitar fondo" if energy_fill_all_visible else "Mostrar fondo"
            act_energy_fill_toggle = energy_style_menu.addAction(fill_text)
            outline_text = (
                "Quitar contorno de cajas"
                if energy_box_outline_all_visible
                else "Mostrar contorno de cajas"
            )
            act_energy_box_outline_toggle = energy_style_menu.addAction(outline_text)
            base_text = (
                "Quitar base de cajas"
                if energy_box_base_all_visible
                else "Mostrar base de cajas"
            )
            act_energy_box_base_toggle = energy_style_menu.addAction(base_text)
            energy_style_menu.addSeparator()
            act_energy_reset_style = energy_style_menu.addAction("Restablecer estilo")
            act_energy_reset_style.setEnabled(energy_has_custom_style)
        if selected_semantic_diagram_items and semantic_diagram_target is not None:
            menu.addSeparator()
            semantic_menu = menu.addMenu("Diagrama electrónico")
            act_semantic_edit_builder = semantic_menu.addAction("Edit Electronic Diagram...")
            summary_lines = [
                str(line).strip()
                for line in list(
                    semantic_diagram_target.semantic_diagram.metadata.get("summary_lines", []) or []
                )
                if str(line).strip()
            ]
            act_semantic_summary_toggle = semantic_menu.addAction("Mostrar resumen inferior")
            act_semantic_summary_toggle.setCheckable(True)
            act_semantic_summary_toggle.setChecked(
                bool(semantic_diagram_target.semantic_diagram.metadata.get("show_summary", True))
            )
            act_semantic_summary_toggle.setEnabled(bool(summary_lines))
            semantic_menu.addSeparator()
            act_semantic_title = semantic_menu.addAction("Editar título...")
            if semantic_diagram_target.semantic_diagram.levels:
                act_semantic_level_label = semantic_menu.addAction("Editar etiqueta de nivel...")
            if semantic_diagram_target.semantic_diagram.lanes:
                act_semantic_lane_label = semantic_menu.addAction("Editar etiqueta de carril...")
        if selected_orbital_items:
            menu.addSeparator()
            orbital_menu = menu.addMenu("Orbital")
            act_orbital_color = orbital_menu.addAction("Color de orbital...")
            act_orbital_reset_color = orbital_menu.addAction("Restablecer color de orbital")
            act_orbital_opacity = orbital_menu.addAction("Opacidad de orbital...")
            if orbital_target is not None:
                orbital_menu.addSeparator()
                lobe_menu = orbital_menu.addMenu("Editar lóbulo")
                act_orbital_lobe_color = lobe_menu.addAction("Color por lóbulo...")
                act_orbital_lobe_opacity = lobe_menu.addAction("Opacidad por lóbulo...")
                act_orbital_lobe_move = lobe_menu.addAction("Mover lóbulo...")
                act_orbital_lobe_reset = lobe_menu.addAction("Restablecer estilo de lóbulo")
        menu.addSeparator()
        act_select_all = menu.addAction("Seleccionar todo")
        menu.addSeparator()
        analysis_menu = menu.addMenu("Análisis")
        act_analysis_name = analysis_menu.addAction("Nombre (SMILES)")
        act_analysis_iupac = analysis_menu.addAction("Nombre (IUPAC)")
        act_analysis_formula = analysis_menu.addAction("Fórmula química")
        act_analysis_exact = analysis_menu.addAction("Masa exacta")
        act_analysis_weight = analysis_menu.addAction("Peso molecular")
        act_analysis_mz = analysis_menu.addAction("m/z")
        act_analysis_elemental = analysis_menu.addAction("Análisis elemental")
        analysis_menu.addSeparator()
        act_analysis_all = analysis_menu.addAction("Todo")

        act_cut.setEnabled(has_selection)
        act_copy.setEnabled(has_selection)
        act_paste.setEnabled(self.can_paste_from_clipboard())
        act_delete.setEnabled(has_selection)
        analysis_menu.setEnabled(bool(self.model.atoms))

        action = menu.exec(global_pos)
        if action is None:
            return
        if action in branch_action_handlers:
            ok, message = branch_action_handlers[action]()
            if not ok:
                self._show_status_message(message)
            return
        if act_thicker is not None and action == act_thicker:
            self._adjust_selected_bond_stroke(self._bond_stroke_step())
            return
        if act_thinner is not None and action == act_thinner:
            self._adjust_selected_bond_stroke(-self._bond_stroke_step())
            return
        if act_reset_thickness is not None and action == act_reset_thickness:
            self._reset_selected_bond_stroke()
            return
        if act_line_equalize is not None and action == act_line_equalize:
            self._equalize_selected_line_arrow_lengths(reference_item=line_arrow_target)
            return
        if act_line_set_length is not None and action == act_line_set_length:
            self._prompt_selected_line_arrow_length(reference_item=line_arrow_target)
            return
        if act_color is not None and action == act_color:
            self._prompt_bond_color()
            return
        if act_reset_color is not None and action == act_reset_color:
            self._set_selected_bond_color(None)
            return
        if act_orbital_color is not None and action == act_orbital_color:
            self._prompt_selected_orbital_color()
            return
        if act_orbital_reset_color is not None and action == act_orbital_reset_color:
            self._reset_selected_orbital_color()
            return
        if act_orbital_opacity is not None and action == act_orbital_opacity:
            self._prompt_selected_orbital_opacity()
            return
        if energy_diagram_target is not None and act_energy_label is not None and action == act_energy_label:
            self._prompt_energy_diagram_label(energy_diagram_target)
            return
        if act_energy_box_count is not None and action == act_energy_box_count:
            self._prompt_selected_energy_diagram_box_count()
            return
        if energy_diagram_target is not None and act_energy_occupancies is not None and action == act_energy_occupancies:
            self._prompt_energy_diagram_occupancies(energy_diagram_target)
            return
        if act_energy_label_left is not None and action == act_energy_label_left:
            self._set_selected_energy_diagram_label_side("left")
            return
        if act_energy_label_right is not None and action == act_energy_label_right:
            self._set_selected_energy_diagram_label_side("right")
            return
        if act_energy_stroke_color is not None and action == act_energy_stroke_color:
            self._prompt_selected_energy_diagram_color("stroke_color", "Color de contorno")
            return
        if act_energy_fill_color is not None and action == act_energy_fill_color:
            self._prompt_selected_energy_diagram_color("fill_color", "Color de relleno")
            return
        if act_energy_label_color is not None and action == act_energy_label_color:
            self._prompt_selected_energy_diagram_color("label_color", "Color de etiqueta")
            return
        if act_energy_arrow_up_color is not None and action == act_energy_arrow_up_color:
            self._prompt_selected_energy_diagram_color("arrow_up_color", "Color de flecha ↑")
            return
        if act_energy_arrow_down_color is not None and action == act_energy_arrow_down_color:
            self._prompt_selected_energy_diagram_color("arrow_down_color", "Color de flecha ↓")
            return
        if act_energy_fill_toggle is not None and action == act_energy_fill_toggle:
            self._set_selected_energy_diagram_fill_visible(not energy_fill_all_visible)
            return
        if act_energy_box_outline_toggle is not None and action == act_energy_box_outline_toggle:
            self._set_selected_energy_diagram_box_stroke_visible(
                not energy_box_outline_all_visible
            )
            return
        if act_energy_box_base_toggle is not None and action == act_energy_box_base_toggle:
            self._set_selected_energy_diagram_box_base_visible(
                not energy_box_base_all_visible
            )
            return
        if act_energy_reset_style is not None and action == act_energy_reset_style:
            self._reset_selected_energy_diagram_style()
            return
        if (
            semantic_diagram_target is not None
            and act_semantic_edit_builder is not None
            and action == act_semantic_edit_builder
        ):
            window = self.window()
            edit_handler = getattr(window, "_edit_selected_semantic_diagram", None)
            if callable(edit_handler):
                edit_handler()
            return
        if (
            semantic_diagram_target is not None
            and act_semantic_summary_toggle is not None
            and action == act_semantic_summary_toggle
        ):
            self._set_semantic_diagram_summary_visible(
                semantic_diagram_target,
                act_semantic_summary_toggle.isChecked(),
            )
            return
        if semantic_diagram_target is not None and act_semantic_title is not None and action == act_semantic_title:
            self._prompt_semantic_diagram_title(semantic_diagram_target)
            return
        if semantic_diagram_target is not None and act_semantic_level_label is not None and action == act_semantic_level_label:
            self._prompt_semantic_diagram_level_label(semantic_diagram_target)
            return
        if semantic_diagram_target is not None and act_semantic_lane_label is not None and action == act_semantic_lane_label:
            self._prompt_semantic_diagram_lane_title(semantic_diagram_target)
            return
        if orbital_target is not None and act_orbital_lobe_color is not None and action == act_orbital_lobe_color:
            self._prompt_orbital_part_color(orbital_target)
            return
        if orbital_target is not None and act_orbital_lobe_opacity is not None and action == act_orbital_lobe_opacity:
            self._prompt_orbital_part_opacity(orbital_target)
            return
        if orbital_target is not None and act_orbital_lobe_move is not None and action == act_orbital_lobe_move:
            self._prompt_orbital_part_offset(orbital_target)
            return
        if orbital_target is not None and act_orbital_lobe_reset is not None and action == act_orbital_lobe_reset:
            self._reset_orbital_part_style(orbital_target)
            return
        if act_anchor is not None and action == act_anchor and isinstance(clicked_item, AtomItem):
            self._prompt_anchor_for_atom(clicked_item.atom_id)
            return
        if act_charge_increment is not None and action == act_charge_increment:
            self._adjust_atom_formal_charge(charge_atom_ids, +1)
            return
        if act_charge_decrement is not None and action == act_charge_decrement:
            self._adjust_atom_formal_charge(charge_atom_ids, -1)
            return
        if act_charge_set is not None and action == act_charge_set:
            self._prompt_set_atom_formal_charge(charge_atom_ids)
            return
        if act_charge_reset is not None and action == act_charge_reset:
            self._set_atom_formal_charge(charge_atom_ids, 0)
            return
        if act_no_implicit is not None and action == act_no_implicit:
            self._set_atom_no_implicit(charge_atom_ids, not charge_all_no_implicit)
            return
        if act_toggle_coord_bond is not None and action == act_toggle_coord_bond:
            self._set_bond_coordination_style(coordination_bond_ids, coordination_bond_enable)
            return
        if act_toggle_coord_sphere is not None and action == act_toggle_coord_sphere:
            self._set_coordination_sphere(coordination_atom_ids, coordination_enable)
            return
        if (
            act_coord_sphere_toggle_transparent is not None
            and action == act_coord_sphere_toggle_transparent
        ):
            self._set_coordination_sphere_transparent(
                styled_coordination_ids,
                not coordination_all_transparent,
            )
            return
        if act_coord_sphere_color is not None and action == act_coord_sphere_color:
            self._prompt_coordination_sphere_color(styled_coordination_ids)
            return
        if act_coord_sphere_reset_color is not None and action == act_coord_sphere_reset_color:
            self._set_coordination_sphere_color(styled_coordination_ids, None)
            return
        if act_coord_sphere_toggle_fill is not None and action == act_coord_sphere_toggle_fill:
            self._set_coordination_sphere_fill(styled_coordination_ids, not coordination_all_filled)
            return
        if act_coord_sphere_radius is not None and action == act_coord_sphere_radius:
            self._prompt_coordination_sphere_radius(styled_coordination_ids)
            return
        if act_arrange_coordination is not None and action == act_arrange_coordination:
            self._auto_arrange_coordination_ligands(arrangement_center_ids)
            return
        if action == act_analysis_name:
            self._run_analysis_action("name", scene_pos)
            return
        if action == act_analysis_iupac:
            self._run_analysis_action("iupac", scene_pos)
            return
        if action == act_analysis_formula:
            self._run_analysis_action("formula", scene_pos)
            return
        if action == act_analysis_exact:
            self._run_analysis_action("exact", scene_pos)
            return
        if action == act_analysis_weight:
            self._run_analysis_action("weight", scene_pos)
            return
        if action == act_analysis_mz:
            self._run_analysis_action("mz", scene_pos)
            return
        if action == act_analysis_elemental:
            self._run_analysis_action("elemental", scene_pos)
            return
        if action == act_analysis_all:
            self._run_analysis_action("all", scene_pos)
            return
        if action == act_copy:
            self.copy_to_clipboard()
            return
        if action == act_paste:
            self.paste_from_clipboard()
            return
        if action == act_cut:
            self.copy_to_clipboard()
            self.delete_selection()
            return
        if action == act_delete:
            self.delete_selection()
            return
        if act_scale_selection is not None and action == act_scale_selection:
            self.open_selection_scale_dialog()
            return
        if action == act_select_all:
            self._select_all_items()

    def _adjust_atom_formal_charge(self, atom_ids: Iterable[int], delta: int) -> None:
        """Ajusta la carga formal de los átomos seleccionados."""
        valid_ids = [atom_id for atom_id in atom_ids if atom_id in self.model.atoms]
        if not valid_ids:
            return
        self.undo_stack.beginMacro("Adjust atom formal charge")
        for atom_id in valid_ids:
            atom = self.model.get_atom(atom_id)
            self.undo_stack.push(
                ChangeChargeCommand(self.model, self, atom_id, int(atom.charge) + int(delta))
            )
        self.undo_stack.endMacro()

    def _set_atom_formal_charge(self, atom_ids: Iterable[int], charge: int) -> None:
        """Establece una carga formal fija para una lista de átomos."""
        valid_ids = [atom_id for atom_id in atom_ids if atom_id in self.model.atoms]
        if not valid_ids:
            return
        normalized_charge = int(charge)
        self.undo_stack.beginMacro("Set atom formal charge")
        for atom_id in valid_ids:
            self.undo_stack.push(ChangeChargeCommand(self.model, self, atom_id, normalized_charge))
        self.undo_stack.endMacro()

    def _prompt_set_atom_formal_charge(self, atom_ids: Iterable[int]) -> None:
        """Solicita la carga formal y la aplica a los átomos seleccionados."""
        valid_ids = [atom_id for atom_id in atom_ids if atom_id in self.model.atoms]
        if not valid_ids:
            return
        first_charge = int(self.model.get_atom(valid_ids[0]).charge)
        value, ok = QInputDialog.getInt(
            self,
            "Carga formal",
            "Carga formal:",
            first_charge,
            -15,
            15,
            1,
        )
        if not ok:
            return
        self._set_atom_formal_charge(valid_ids, value)

    def _set_atom_no_implicit(self, atom_ids: Iterable[int], enabled: bool) -> None:
        """Activa/desactiva la bandera `no_implicit` para átomos seleccionados."""
        valid_ids = [atom_id for atom_id in atom_ids if atom_id in self.model.atoms]
        if not valid_ids:
            return
        enabled_flag = bool(enabled)
        self.undo_stack.beginMacro("Set no implicit hydrogens")
        for atom_id in valid_ids:
            self.undo_stack.push(
                ChangeNoImplicitCommand(self.model, self, atom_id, enabled_flag)
            )
        self.undo_stack.endMacro()

    def _set_coordination_sphere(self, atom_ids: Iterable[int], enabled: bool) -> None:
        """Activa/desactiva esfera de coordinación para una lista de átomos."""
        enabled_flag = bool(enabled)
        valid_ids = [
            atom_id
            for atom_id in atom_ids
            if atom_id in self.model.atoms
            and bool(getattr(self.model.get_atom(atom_id), "is_coordination_center", False))
            != enabled_flag
        ]
        if not valid_ids:
            return
        self.undo_stack.beginMacro(
            "Enable coordination spheres" if enabled_flag else "Disable coordination spheres"
        )
        for atom_id in valid_ids:
            self.undo_stack.push(
                SetCoordinationCenterCommand(
                    self.model,
                    self,
                    atom_id,
                    enabled_flag,
                )
            )
        self.undo_stack.endMacro()

    def _prompt_coordination_sphere_color(self, atom_ids: Iterable[int]) -> None:
        """Solicita color de esfera y lo aplica a los átomos seleccionados."""
        valid_ids = [atom_id for atom_id in atom_ids if atom_id in self.model.atoms]
        if not valid_ids:
            return
        initial = QColor("#D9DDE3")
        for atom_id in valid_ids:
            atom = self.model.get_atom(atom_id)
            color = getattr(atom, "sphere_color", None)
            if color:
                candidate = QColor(color)
                if candidate.isValid():
                    initial = candidate
                    break
        color = QColorDialog.getColor(initial, self, "Seleccionar color de esfera")
        if not color.isValid():
            return
        self._set_coordination_sphere_color(valid_ids, color.name())

    def _set_coordination_sphere_color(
        self, atom_ids: Iterable[int], color: Optional[str]
    ) -> None:
        """Aplica color de esfera a una lista de átomos."""
        valid_ids = [atom_id for atom_id in atom_ids if atom_id in self.model.atoms]
        if not valid_ids:
            return
        self.undo_stack.beginMacro("Set coordination sphere color")
        for atom_id in valid_ids:
            self.undo_stack.push(
                ChangeCoordinationSphereStyleCommand(
                    self.model,
                    self,
                    atom_id,
                    new_color=color,
                )
            )
        self.undo_stack.endMacro()

    def _set_coordination_sphere_transparent(
        self, atom_ids: Iterable[int], transparent: bool
    ) -> None:
        """Activa/desactiva esfera transparente (solo etiqueta)."""
        valid_ids = [atom_id for atom_id in atom_ids if atom_id in self.model.atoms]
        if not valid_ids:
            return
        self.undo_stack.beginMacro("Toggle coordination sphere transparency")
        for atom_id in valid_ids:
            self.undo_stack.push(
                ChangeCoordinationSphereStyleCommand(
                    self.model,
                    self,
                    atom_id,
                    new_transparent=bool(transparent),
                )
            )
        self.undo_stack.endMacro()

    def _set_coordination_sphere_fill(self, atom_ids: Iterable[int], filled: bool) -> None:
        """Activa o desactiva el relleno de la esfera de coordinación."""
        valid_ids = [atom_id for atom_id in atom_ids if atom_id in self.model.atoms]
        if not valid_ids:
            return
        self.undo_stack.beginMacro("Toggle coordination sphere fill")
        for atom_id in valid_ids:
            self.undo_stack.push(
                ChangeCoordinationSphereStyleCommand(
                    self.model,
                    self,
                    atom_id,
                    new_filled=bool(filled),
                )
            )
        self.undo_stack.endMacro()

    def _prompt_coordination_sphere_radius(self, atom_ids: Iterable[int]) -> None:
        """Solicita un radio de esfera y lo aplica a los átomos seleccionados."""
        valid_ids = [atom_id for atom_id in atom_ids if atom_id in self.model.atoms]
        if not valid_ids:
            return
        current = 16.0
        for atom_id in valid_ids:
            atom = self.model.get_atom(atom_id)
            radius = getattr(atom, "sphere_radius", None)
            if radius is not None:
                current = float(radius)
                break
        value, ok = QInputDialog.getDouble(
            self,
            "Tamaño de esfera",
            "Radio (px):",
            current,
            4.0,
            200.0,
            1,
        )
        if not ok:
            return
        self.undo_stack.beginMacro("Set coordination sphere radius")
        for atom_id in valid_ids:
            self.undo_stack.push(
                ChangeCoordinationSphereStyleCommand(
                    self.model,
                    self,
                    atom_id,
                    new_radius=float(value),
                )
            )
        self.undo_stack.endMacro()

    def _set_bond_coordination_style(self, bond_ids: Iterable[int], enabled: bool) -> None:
        """Convierte enlaces seleccionados entre normal y coordinativo."""
        target_style = BondStyle.COORDINATION if enabled else BondStyle.PLAIN
        valid_ids = [bond_id for bond_id in bond_ids if bond_id in self.model.bonds]
        if not valid_ids:
            return
        changed_ids = [
            bond_id
            for bond_id in valid_ids
            if self.model.get_bond(bond_id).style != target_style
        ]
        if not changed_ids:
            return
        self.undo_stack.beginMacro(
            "Set coordination bond style" if enabled else "Set normal bond style"
        )
        for bond_id in changed_ids:
            bond = self.model.get_bond(bond_id)
            donor_atom_id = None
            if enabled:
                donor_atom_id = self._infer_coordination_donor_atom(
                    bond.a1_id,
                    bond.a2_id,
                    preferred=getattr(bond, "donor_atom_id", None),
                )
            self.undo_stack.push(
                ChangeBondCommand(
                    self.model,
                    self,
                    bond_id,
                    new_order=1,
                    new_style=target_style,
                    new_stereo=BondStereo.NONE,
                    new_is_aromatic=False,
                    new_donor_atom_id=donor_atom_id,
                )
            )
        self.undo_stack.endMacro()

    def _coordination_ligands_for_center(self, center_atom_id: int) -> list[int]:
        """Lista ligandos conectados al centro mediante enlaces coordinativos."""
        if center_atom_id not in self.model.atoms:
            return []
        ligands: list[int] = []
        seen: set[int] = set()
        for bond in self.model.bonds.values():
            if bond.style != BondStyle.COORDINATION:
                continue
            other_id: Optional[int] = None
            if bond.a1_id == center_atom_id:
                other_id = bond.a2_id
            elif bond.a2_id == center_atom_id:
                other_id = bond.a1_id
            if other_id is None or other_id == center_atom_id:
                continue
            if other_id not in self.model.atoms or other_id in seen:
                continue
            seen.add(other_id)
            ligands.append(other_id)
        return ligands

    @staticmethod
    def _coordination_arrangement_step_deg(count: int) -> float:
        """Devuelve separación angular para distribuir ligandos."""
        if count <= 1:
            return 360.0
        if count == 2:
            return 180.0
        if count == 3:
            return 120.0
        if count == 4:
            return 90.0
        if count == 6:
            return 60.0
        return 360.0 / float(count)

    def _coordination_arrangement_radius(self, center_id: int, ligand_ids: list[int]) -> float:
        """Calcula radio objetivo para auto-arrange de ligandos."""
        center = self.model.get_atom(center_id)
        distances: list[float] = []
        for ligand_id in ligand_ids:
            ligand = self.model.atoms.get(ligand_id)
            if ligand is None:
                continue
            distances.append(math.hypot(ligand.x - center.x, ligand.y - center.y))
        base_length = max(18.0, float(self.state.bond_length))
        if not distances:
            return base_length
        avg = sum(distances) / float(len(distances))
        target = max(base_length * 0.8, avg)
        return min(target, base_length * 2.2)

    def _auto_arrange_coordination_ligands(self, center_atom_ids: Iterable[int]) -> None:
        """Distribuye ligandos coordinativos alrededor de uno o más centros."""
        centers = [
            atom_id
            for atom_id in center_atom_ids
            if atom_id in self.model.atoms
            and bool(getattr(self.model.get_atom(atom_id), "is_coordination_center", False))
        ]
        if not centers:
            return
        before: Dict[int, Tuple[float, float]] = {}
        after: Dict[int, Tuple[float, float]] = {}
        moved_ligands: set[int] = set()
        for center_id in centers:
            center = self.model.get_atom(center_id)
            ligand_ids = self._coordination_ligands_for_center(center_id)
            if len(ligand_ids) < 2:
                continue
            ordered_ligands = sorted(
                ligand_ids,
                key=lambda atom_id: math.atan2(
                    self.model.get_atom(atom_id).y - center.y,
                    self.model.get_atom(atom_id).x - center.x,
                ),
            )
            step_deg = self._coordination_arrangement_step_deg(len(ordered_ligands))
            first_atom = self.model.get_atom(ordered_ligands[0])
            start_angle_deg = math.degrees(
                math.atan2(first_atom.y - center.y, first_atom.x - center.x)
            )
            radius = self._coordination_arrangement_radius(center_id, ordered_ligands)

            for index, ligand_id in enumerate(ordered_ligands):
                if ligand_id in moved_ligands:
                    continue
                ligand = self.model.get_atom(ligand_id)
                angle_deg = start_angle_deg + step_deg * index
                angle_rad = math.radians(angle_deg)
                new_x = center.x + radius * math.cos(angle_rad)
                new_y = center.y + radius * math.sin(angle_rad)
                if abs(new_x - ligand.x) < 0.05 and abs(new_y - ligand.y) < 0.05:
                    continue
                before[ligand_id] = (ligand.x, ligand.y)
                after[ligand_id] = (new_x, new_y)
                moved_ligands.add(ligand_id)
        if not after:
            return
        self.undo_stack.push(MoveAtomsCommand(self.model, self, before, after))

    def _bond_stroke_step(self) -> float:
        """Método auxiliar para  bond stroke step.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        return max(0.6, self.drawing_style.stroke_px * 0.35)

    def increase_selected_bond_thickness(self) -> None:
        """Método auxiliar para increase selected bond thickness.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._adjust_selected_bond_stroke(self._bond_stroke_step())

    def decrease_selected_bond_thickness(self) -> None:
        """Método auxiliar para decrease selected bond thickness.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._adjust_selected_bond_stroke(-self._bond_stroke_step())

    def reset_selected_bond_thickness(self) -> None:
        """Método auxiliar para reset selected bond thickness.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._reset_selected_bond_stroke()

    def _prompt_bond_color(self) -> None:
        """Método auxiliar para  prompt bond color.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        bond_ids = list(self.state.selected_bonds)
        if not bond_ids:
            return
        initial = QColor(self.drawing_style.bond_color)
        for bond_id in bond_ids:
            bond = self.model.get_bond(bond_id)
            if bond.color:
                initial = QColor(bond.color)
                break
        color = QColorDialog.getColor(initial, self, "Seleccionar color de enlace")
        if not color.isValid():
            return
        self._set_selected_bond_color(color.name())

    def _set_selected_bond_color(self, color: Optional[str]) -> None:
        """Método auxiliar para  set selected bond color.

        Args:
            color: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        bond_ids = list(self.state.selected_bonds)
        if not bond_ids:
            return
        self.undo_stack.beginMacro("Change bond color")
        for bond_id in bond_ids:
            cmd = ChangeBondColorCommand(self.model, self, bond_id, color)
            self.undo_stack.push(cmd)
        self.undo_stack.endMacro()

    @staticmethod
    def _orbital_part_display_name(name: str) -> str:
        labels = {
            "sphere": "Esfera",
            "top": "Lóbulo superior",
            "bottom": "Lóbulo inferior",
            "left": "Lóbulo izquierdo",
            "right": "Lóbulo derecho",
            "major": "Lóbulo principal",
            "minor": "Lóbulo secundario",
            "minor_left": "Lóbulo menor izquierdo",
            "minor_right": "Lóbulo menor derecho",
            "bond": "Nube central",
            "upper": "Nube superior",
            "lower": "Nube inferior",
            "ring": "Toroide",
            "ring_upper": "Toroide superior",
            "ring_lower": "Toroide inferior",
        }
        if name.startswith("lobe_"):
            suffix = name.split("_", 1)[1]
            if suffix.isdigit():
                return f"Lóbulo {int(suffix) + 1}"
        return labels.get(name, name.replace("_", " "))

    def _single_orbital_target(
        self,
        clicked_item: Optional[QGraphicsItem] = None,
    ) -> Optional[OrbitalAnnotationItem]:
        if isinstance(clicked_item, OrbitalAnnotationItem):
            return clicked_item
        selected = self._selected_orbital_items()
        return selected[0] if len(selected) == 1 else None

    def _choose_orbital_part(
        self,
        item: OrbitalAnnotationItem,
        *,
        title: str = "Seleccionar lóbulo",
    ) -> Optional[str]:
        part_names = item.part_names()
        if not part_names:
            return None
        display_to_name = {
            self._orbital_part_display_name(name): name
            for name in part_names
        }
        value, ok = QInputDialog.getItem(
            self,
            title,
            "Lóbulo:",
            list(display_to_name.keys()),
            0,
            False,
        )
        if not ok:
            return None
        return display_to_name.get(value)

    def _push_orbital_style_change(
        self,
        updates: dict[OrbitalAnnotationItem, dict[str, dict[str, object]]],
        *,
        text: str,
    ) -> bool:
        before: dict[OrbitalAnnotationItem, dict[str, dict[str, object]]] = {}
        after: dict[OrbitalAnnotationItem, dict[str, dict[str, object]]] = {}
        for item, payload in updates.items():
            if item.scene() is not self.scene:
                continue
            current = item.part_styles()
            if current == payload:
                continue
            before[item] = current
            after[item] = payload
        if not after:
            return False
        self.undo_stack.push(StyleOrbitalItemsCommand(self, before, after, text))
        return True

    def _prompt_selected_orbital_color(self) -> None:
        orbitals = self._selected_orbital_items()
        if not orbitals:
            return
        initial = QColor("#111111")
        for item in orbitals:
            for name in item.part_names():
                color_value = item.effective_part_style(name).get("color")
                if color_value:
                    candidate = QColor(str(color_value))
                    if candidate.isValid():
                        initial = candidate
                        break
            else:
                continue
            break
        color = QColorDialog.getColor(initial, self, "Seleccionar color de orbital")
        if not color.isValid():
            return
        updates: dict[OrbitalAnnotationItem, dict[str, dict[str, object]]] = {}
        for item in orbitals:
            payload = item.part_styles()
            for name in item.part_names():
                part_payload = dict(payload.get(name, {}))
                part_payload["color"] = color.name(QColor.NameFormat.HexRgb)
                payload[name] = part_payload
            updates[item] = payload
        self._push_orbital_style_change(updates, text="Change orbital color")

    def _reset_selected_orbital_color(self) -> None:
        orbitals = self._selected_orbital_items()
        if not orbitals:
            return
        updates: dict[OrbitalAnnotationItem, dict[str, dict[str, object]]] = {}
        for item in orbitals:
            payload = item.part_styles()
            for name in list(payload.keys()):
                part_payload = dict(payload.get(name, {}))
                part_payload.pop("color", None)
                if part_payload:
                    payload[name] = part_payload
                else:
                    payload.pop(name, None)
            updates[item] = payload
        self._push_orbital_style_change(updates, text="Reset orbital color")

    def _prompt_selected_orbital_opacity(self) -> None:
        orbitals = self._selected_orbital_items()
        if not orbitals:
            return
        current = 100
        values = [self.effective_item_opacity(item) for item in orbitals if item.scene() is self.scene]
        if values:
            current = int(round(sum(values) / float(len(values)) * 100.0))
        value, ok = QInputDialog.getInt(
            self,
            "Opacidad de orbital",
            "Opacidad (%):",
            current,
            0,
            100,
            1,
        )
        if not ok:
            return
        target = max(0.0, min(1.0, float(value) / 100.0))
        item_values = {item: target for item in orbitals if item.scene() is self.scene}
        if not item_values:
            return
        self.undo_stack.push(
            ChangeCanvasOpacityCommand(
                self.model,
                self,
                item_values=item_values,
                text="Change orbital opacity",
            )
        )

    def _prompt_orbital_part_color(self, item: OrbitalAnnotationItem) -> None:
        part_name = self._choose_orbital_part(item, title="Color por lóbulo")
        if not part_name:
            return
        current = item.effective_part_style(part_name).get("color")
        initial = QColor(str(current)) if current else QColor("#111111")
        color = QColorDialog.getColor(initial, self, f"Color de {self._orbital_part_display_name(part_name)}")
        if not color.isValid():
            return
        payload = item.part_styles()
        part_payload = dict(payload.get(part_name, {}))
        part_payload["color"] = color.name(QColor.NameFormat.HexRgb)
        payload[part_name] = part_payload
        self._push_orbital_style_change({item: payload}, text="Change orbital lobe color")

    def _prompt_orbital_part_opacity(self, item: OrbitalAnnotationItem) -> None:
        part_name = self._choose_orbital_part(item, title="Opacidad por lóbulo")
        if not part_name:
            return
        current = item.effective_part_style(part_name).get("opacity")
        percent = 100 if current is None else int(round(float(current) * 100.0))
        value, ok = QInputDialog.getInt(
            self,
            "Opacidad por lóbulo",
            f"{self._orbital_part_display_name(part_name)} (%):",
            percent,
            0,
            100,
            1,
        )
        if not ok:
            return
        payload = item.part_styles()
        part_payload = dict(payload.get(part_name, {}))
        opacity = max(0.0, min(1.0, float(value) / 100.0))
        if abs(opacity - 1.0) <= 1e-6:
            part_payload.pop("opacity", None)
        else:
            part_payload["opacity"] = opacity
        if part_payload:
            payload[part_name] = part_payload
        else:
            payload.pop(part_name, None)
        self._push_orbital_style_change({item: payload}, text="Change orbital lobe opacity")

    def _prompt_orbital_part_offset(self, item: OrbitalAnnotationItem) -> None:
        part_name = self._choose_orbital_part(item, title="Mover lóbulo")
        if not part_name:
            return
        current = item.effective_part_style(part_name)
        current_x = float(current.get("offset_x", 0.0) or 0.0)
        current_y = float(current.get("offset_y", 0.0) or 0.0)
        offset_x, ok = QInputDialog.getDouble(
            self,
            "Mover lóbulo",
            f"{self._orbital_part_display_name(part_name)} desplazamiento X:",
            current_x,
            -80.0,
            80.0,
            1,
        )
        if not ok:
            return
        offset_y, ok = QInputDialog.getDouble(
            self,
            "Mover lóbulo",
            f"{self._orbital_part_display_name(part_name)} desplazamiento Y:",
            current_y,
            -80.0,
            80.0,
            1,
        )
        if not ok:
            return
        payload = item.part_styles()
        part_payload = dict(payload.get(part_name, {}))
        if abs(offset_x) <= 1e-6:
            part_payload.pop("offset_x", None)
        else:
            part_payload["offset_x"] = float(offset_x)
        if abs(offset_y) <= 1e-6:
            part_payload.pop("offset_y", None)
        else:
            part_payload["offset_y"] = float(offset_y)
        if part_payload:
            payload[part_name] = part_payload
        else:
            payload.pop(part_name, None)
        self._push_orbital_style_change({item: payload}, text="Move orbital lobe")

    def _reset_orbital_part_style(self, item: OrbitalAnnotationItem) -> None:
        part_name = self._choose_orbital_part(item, title="Restablecer lóbulo")
        if not part_name:
            return
        payload = item.part_styles()
        payload.pop(part_name, None)
        self._push_orbital_style_change({item: payload}, text="Reset orbital lobe style")

    def _single_energy_diagram_target(
        self,
        clicked_item: Optional[QGraphicsItem] = None,
    ) -> Optional[EnergyDiagramItem]:
        """Devuelve un diagrama único para acciones de edición puntual."""
        if isinstance(clicked_item, EnergyDiagramItem):
            return clicked_item
        selected = self._selected_energy_diagram_items()
        return selected[0] if len(selected) == 1 else None

    def _single_semantic_diagram_target(
        self,
        clicked_item: Optional[QGraphicsItem] = None,
    ) -> Optional[CompositeDiagramItem]:
        """Devuelve un diagrama semántico único para acciones de edición."""
        if isinstance(clicked_item, CompositeDiagramItem):
            return clicked_item
        selected = self._selected_semantic_diagram_items()
        return selected[0] if len(selected) == 1 else None

    def _set_semantic_diagram_summary_visible(
        self,
        item: CompositeDiagramItem,
        visible: bool,
    ) -> bool:
        """Alterna la visibilidad del resumen textual inferior con undo/redo."""
        if item is None or item.scene() is not self.scene:
            return False
        before_payload = item.to_json()
        before_payload["opacity"] = self.item_raw_opacity(item)
        current_visible = bool(
            item.semantic_diagram.metadata.get("show_summary", True)
        )
        if current_visible == bool(visible):
            return False
        after_payload = deepcopy(before_payload)
        semantic_payload = dict(after_payload.get("semantic_diagram", {}) or {})
        metadata = dict(semantic_payload.get("metadata", {}) or {})
        metadata["show_summary"] = bool(visible)
        semantic_payload["metadata"] = metadata
        after_payload["semantic_diagram"] = semantic_payload
        self.undo_stack.push(
            EditSemanticDiagramCommand(
                self,
                item,
                before_payload,
                after_payload,
                text="Toggle semantic diagram summary",
            )
        )
        return True

    @staticmethod
    def _plain_text_from_markup(value: object) -> str:
        raw = str(value or "")
        if not raw:
            return ""
        if not Qt.mightBeRichText(raw):
            return raw
        document = QTextDocument()
        document.setHtml(raw)
        return document.toPlainText()

    def _prompt_rich_text_value(
        self,
        title: str,
        label: str,
        initial_text: str,
    ) -> tuple[str, bool]:
        """Abre el editor enriquecido del main window o usa fallback plano."""
        window = self.window()
        prompt = getattr(window, "_open_rich_text_value_dialog", None)
        if callable(prompt):
            return prompt(title=title, label=label, initial_text=initial_text)
        value, ok = QInputDialog.getText(self, title, label, text=str(initial_text or ""))
        return str(value or ""), bool(ok)

    def _prompt_semantic_diagram_title(self, item: CompositeDiagramItem) -> None:
        """Solicita un nuevo título para el diagrama semántico."""
        value, ok = self._prompt_rich_text_value(
            "Título del diagrama",
            "Título:",
            str(item.semantic_diagram.title or ""),
        )
        if not ok:
            return
        if item.set_diagram_title(str(value or "")):
            self._sync_selection_from_scene()

    def _prompt_semantic_diagram_lane_title(
        self,
        item: CompositeDiagramItem,
        lane_id: Optional[str] = None,
    ) -> None:
        """Edita la etiqueta de un carril del diagrama semántico."""
        target_lane_id = str(lane_id or "")
        lane = None
        if target_lane_id:
            lane = next(
                (candidate for candidate in item.semantic_diagram.lanes if candidate.id == target_lane_id),
                None,
            )
        if lane is None:
            lanes = [candidate for candidate in item.semantic_diagram.lanes]
            if not lanes:
                return
            lane_labels = [
                self._plain_text_from_markup(candidate.title or candidate.id or f"Carril {index + 1}")
                for index, candidate in enumerate(lanes)
            ]
            selected_label, ok = QInputDialog.getItem(
                self,
                "Etiqueta de carril",
                "Carril:",
                lane_labels,
                0,
                False,
            )
            if not ok:
                return
            chosen_index = lane_labels.index(selected_label)
            lane = lanes[chosen_index]
        value, ok = self._prompt_rich_text_value(
            "Etiqueta de carril",
            "Etiqueta:",
            str(lane.title or ""),
        )
        if not ok:
            return
        if item.set_lane_title(lane.id, str(value or "")):
            self._sync_selection_from_scene()

    def _prompt_semantic_diagram_level_label(
        self,
        item: CompositeDiagramItem,
        level_id: Optional[str] = None,
    ) -> None:
        """Edita la etiqueta de un nivel del diagrama semántico."""
        target_level_id = str(level_id or "")
        level = None
        if target_level_id:
            level = next(
                (candidate for candidate in item.semantic_diagram.levels if candidate.id == target_level_id),
                None,
            )
        if level is None:
            levels = [candidate for candidate in item.semantic_diagram.levels]
            if not levels:
                return
            level_labels = [
                self._plain_text_from_markup(candidate.label or candidate.id or f"Nivel {index + 1}")
                for index, candidate in enumerate(levels)
            ]
            selected_label, ok = QInputDialog.getItem(
                self,
                "Etiqueta de nivel",
                "Nivel:",
                level_labels,
                0,
                False,
            )
            if not ok:
                return
            chosen_index = level_labels.index(selected_label)
            level = levels[chosen_index]
        value, ok = self._prompt_rich_text_value(
            "Etiqueta de nivel",
            "Etiqueta:",
            str(level.label or ""),
        )
        if not ok:
            return
        if item.set_level_label(level.id, str(value or "")):
            self._sync_selection_from_scene()

    def _push_energy_diagram_config_change(
        self,
        updates: dict[EnergyDiagramItem, dict[str, object]],
        *,
        text: str,
    ) -> bool:
        """Empaqueta cambios de configuración de diagramas en undo/redo."""
        before: dict[EnergyDiagramItem, dict[str, object]] = {}
        after: dict[EnergyDiagramItem, dict[str, object]] = {}
        for item, payload in updates.items():
            if item.scene() is not self.scene:
                continue
            current = item.config_payload()
            if current == payload:
                continue
            before[item] = current
            after[item] = payload
        if not after:
            return False
        self.undo_stack.push(ConfigureEnergyDiagramItemsCommand(self, before, after, text))
        return True

    def _prompt_energy_diagram_label(self, item: EnergyDiagramItem) -> None:
        """Solicita una nueva etiqueta para un diagrama."""
        if not item.supports_free_label():
            return
        value, ok = self._prompt_rich_text_value(
            "Etiqueta de diagrama",
            "Etiqueta:",
            item.label(),
        )
        if not ok:
            return
        payload = item.config_payload()
        payload["label"] = str(value or "")
        self._push_energy_diagram_config_change(
            {item: payload},
            text="Change energy diagram label",
        )

    def _prompt_energy_diagram_occupancies(self, item: EnergyDiagramItem) -> None:
        """Edita la ocupación electrónica caja por caja."""
        value, ok = QInputDialog.getText(
            self,
            "Ocupación electrónica",
            "Ocupaciones separadas por comas (empty, up, down, pair, upup, downdown):",
            text=", ".join(item.occupancies()),
        )
        if not ok:
            return
        payload = item.config_payload()
        payload["occupancies"] = list(
            normalize_energy_occupancies(
                value,
                kind=item.kind(),
                box_count=item.slot_count(),
            )
        )
        self._push_energy_diagram_config_change(
            {item: payload},
            text="Change energy diagram occupancies",
        )

    def _prompt_selected_energy_diagram_box_count(self) -> None:
        """Solicita un nuevo numero de cajas para filas de orbitales."""
        items = [item for item in self._selected_energy_diagram_items() if item.family() == "row"]
        if not items:
            return
        current = int(items[0].slot_count())
        value, ok = QInputDialog.getInt(
            self,
            "Numero de cajas",
            "Numero de cajas:",
            current,
            1,
            20,
        )
        if not ok:
            return
        updates: dict[EnergyDiagramItem, dict[str, object]] = {}
        for item in items:
            payload = item.config_payload()
            payload["slot_count"] = int(value)
            updates[item] = payload
        self._push_energy_diagram_config_change(
            updates,
            text="Change energy diagram box count",
        )

    def _set_selected_energy_diagram_label_side(self, side: str) -> None:
        """Aplica la posición de etiqueta a los diagramas seleccionados."""
        items = self._selected_energy_diagram_items()
        if not items:
            return
        updates: dict[EnergyDiagramItem, dict[str, object]] = {}
        for item in items:
            if not item.supports_free_label():
                continue
            payload = item.config_payload()
            payload["label_side"] = side
            updates[item] = payload
        self._push_energy_diagram_config_change(
            updates,
            text="Change energy diagram label side",
        )

    def _prompt_selected_energy_diagram_color(self, key: str, title: str) -> None:
        """Solicita un color para un canal visual de diagramas seleccionados."""
        items = self._selected_energy_diagram_items()
        if not items:
            return
        initial = QColor("#222222")
        for item in items:
            candidate = QColor(str(item.effective_style().get(key, "")))
            if candidate.isValid():
                initial = candidate
                break
        color = QColorDialog.getColor(initial, self, title)
        if not color.isValid():
            return
        updates: dict[EnergyDiagramItem, dict[str, object]] = {}
        for item in items:
            payload = item.config_payload()
            style_payload = dict(payload.get("style_payload", {}) or {})
            style_payload[key] = color.name(QColor.NameFormat.HexRgb)
            payload["style_payload"] = style_payload
            updates[item] = payload
        self._push_energy_diagram_config_change(updates, text=f"Change {key}")

    def _set_selected_energy_diagram_fill_visible(self, visible: bool) -> None:
        """Muestra u oculta el fondo de los diagramas seleccionados."""
        items = self._selected_energy_diagram_items()
        if not items:
            return
        updates: dict[EnergyDiagramItem, dict[str, object]] = {}
        for item in items:
            payload = item.config_payload()
            style_payload = dict(payload.get("style_payload", {}) or {})
            style_payload["fill_visible"] = bool(visible)
            payload["style_payload"] = style_payload
            updates[item] = payload
        self._push_energy_diagram_config_change(
            updates,
            text="Toggle energy diagram fill",
        )

    def _set_selected_energy_diagram_box_stroke_visible(self, visible: bool) -> None:
        """Muestra u oculta el contorno de cajas/niveles editables."""
        items = self._selected_energy_diagram_items()
        if not items:
            return
        updates: dict[EnergyDiagramItem, dict[str, object]] = {}
        for item in items:
            payload = item.config_payload()
            style_payload = dict(payload.get("style_payload", {}) or {})
            style_payload["box_stroke_visible"] = bool(visible)
            if visible:
                style_payload["box_base_visible"] = False
            payload["style_payload"] = style_payload
            updates[item] = payload
        self._push_energy_diagram_config_change(
            updates,
            text="Toggle energy diagram box outline",
        )

    def _set_selected_energy_diagram_box_base_visible(self, visible: bool) -> None:
        """Muestra u oculta la base horizontal de las cajas."""
        items = self._selected_energy_diagram_items()
        if not items:
            return
        updates: dict[EnergyDiagramItem, dict[str, object]] = {}
        for item in items:
            payload = item.config_payload()
            style_payload = dict(payload.get("style_payload", {}) or {})
            style_payload["box_base_visible"] = bool(visible)
            if visible:
                style_payload["box_stroke_visible"] = False
            payload["style_payload"] = style_payload
            updates[item] = payload
        self._push_energy_diagram_config_change(
            updates,
            text="Toggle energy diagram box base",
        )

    def _reset_selected_energy_diagram_style(self) -> None:
        """Elimina overrides de color/fondo de los diagramas seleccionados."""
        items = self._selected_energy_diagram_items()
        if not items:
            return
        updates: dict[EnergyDiagramItem, dict[str, object]] = {}
        for item in items:
            payload = item.config_payload()
            payload["style_payload"] = {}
            updates[item] = payload
        self._push_energy_diagram_config_change(
            updates,
            text="Reset energy diagram style",
        )

    def _adjust_selected_bond_stroke(self, delta: float) -> None:
        """Método auxiliar para  adjust selected bond stroke.

        Args:
            delta: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        bond_ids = list(self.state.selected_bonds)
        arrow_items = self._selected_arrow_items()
        bracket_items = self._selected_bracket_items()
        if not bond_ids and not arrow_items and not bracket_items:
            return
        default_stroke = self.drawing_style.stroke_px
        self.undo_stack.beginMacro("Change bond/arrow/bracket thickness")
        for bond_id in bond_ids:
            bond = self.model.get_bond(bond_id)
            current = bond.stroke_px if bond.stroke_px is not None else default_stroke
            new_value = max(0.6, current + delta)
            if abs(new_value - default_stroke) < 0.05:
                new_value = None
            cmd = ChangeBondStrokeCommand(self.model, self, bond_id, new_value)
            self.undo_stack.push(cmd)
        for item in arrow_items:
            current = item.stroke_px() if item.stroke_px() is not None else default_stroke
            new_value = max(0.6, current + delta)
            if abs(new_value - default_stroke) < 0.05:
                new_value = None
            self.undo_stack.push(ChangeArrowStrokeCommand(self, item, new_value))
        for item in bracket_items:
            current = item.stroke_px() if item.stroke_px() is not None else default_stroke
            new_value = max(0.6, current + delta)
            if abs(new_value - default_stroke) < 0.05:
                new_value = None
            self.undo_stack.push(ChangeBracketStrokeCommand(self, item, new_value))
        self.undo_stack.endMacro()

    def _reset_selected_bond_stroke(self) -> None:
        """Método auxiliar para  reset selected bond stroke.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        bond_ids = list(self.state.selected_bonds)
        arrow_items = self._selected_arrow_items()
        bracket_items = self._selected_bracket_items()
        if not bond_ids and not arrow_items and not bracket_items:
            return
        self.undo_stack.beginMacro("Reset bond/arrow/bracket thickness")
        for bond_id in bond_ids:
            cmd = ChangeBondStrokeCommand(self.model, self, bond_id, None)
            self.undo_stack.push(cmd)
        for item in arrow_items:
            self.undo_stack.push(ChangeArrowStrokeCommand(self, item, None))
        for item in bracket_items:
            self.undo_stack.push(ChangeBracketStrokeCommand(self, item, None))
        self.undo_stack.endMacro()

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

    def _handle_hotkeys(self, event) -> bool:
        """Método auxiliar para  handle hotkeys.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        key = event.key()
        modifiers = event.modifiers()

        if (
            key == Qt.Key.Key_0
            and (modifiers & Qt.KeyboardModifier.AltModifier)
            and not (
                modifiers
                & (
                    Qt.KeyboardModifier.ControlModifier
                    | Qt.KeyboardModifier.MetaModifier
                    | Qt.KeyboardModifier.ShiftModifier
                )
            )
        ):
            if self._is_rotating_3d:
                self._finalize_3d_rotation_drag()
            self._reset_3d_rotation_perspective()
            return True

        if Qt.Key.Key_0 <= key <= Qt.Key.Key_9:
            digit = key - Qt.Key.Key_0
            if digit == 0:
                return False
            if self.hovered_atom_id is not None:
                if modifiers & Qt.KeyboardModifier.ShiftModifier:
                    if digit >= 3:
                        self._create_ring_hotkey(self.hovered_atom_id, digit)
                        return True
                else:
                    self._create_chain_hotkey(self.hovered_atom_id, digit)
                    return True
            if self.hovered_bond_id is not None:
                if modifiers & Qt.KeyboardModifier.ShiftModifier:
                    if digit >= 3:
                        self._create_ring_from_bond_hotkey(self.hovered_bond_id, digit)
                        return True
                else:
                    self._set_bond_order_hotkey(self.hovered_bond_id, digit)
                    return True
        return False

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
