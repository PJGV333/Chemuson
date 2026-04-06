from __future__ import annotations

from ._shared import (
    ArrowItem,
    AtomItem,
    BRANCH_ROTATION_STEP_DEG,
    BondItem,
    BondStereo,
    BondStyle,
    BracketItem,
    ChangeArrowStrokeCommand,
    ChangeBondColorCommand,
    ChangeBondCommand,
    ChangeBondStrokeCommand,
    ChangeBracketStrokeCommand,
    ChangeChargeCommand,
    ChangeCoordinationSphereStyleCommand,
    ChangeNoImplicitCommand,
    CompositeDiagramItem,
    Dict,
    ELEMENT_SYMBOLS,
    EnergyDiagramItem,
    GelElectrophoresisItem,
    ImageAnnotationItem,
    Iterable,
    MoveAtomsCommand,
    Optional,
    OrbitalAnnotationItem,
    QColor,
    QColorDialog,
    QGraphicsItem,
    QInputDialog,
    QMenu,
    QPointF,
    SetCoordinationCenterCommand,
    TLCPlateItem,
    Tuple,
    angle_deg,
    math,
)

class CanvasContextMenuMixin:
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
