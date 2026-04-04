from __future__ import annotations

from ._shared import *

class CanvasStructureMixin:
    def get_persistence_data(self) -> dict:
        """Collects canvas settings and non-structural items for serialization."""
        data = {
            "settings": {
                "paper_width": self.paper_width,
                "paper_height": self.paper_height,
                "show_grid": self.show_grid,
                "show_rulers": self.show_rulers,
                "bond_length": self.state.bond_length,
                "use_aromatic_circles": self.state.use_aromatic_circles,
                "show_carbons": self.state.show_implicit_carbons,
                "show_hydrogens": self.state.show_implicit_hydrogens,
                "use_element_colors": self.state.use_element_colors,
                "font_family": self.state.label_font_family,
                "font_size": self.state.label_font_size,
                "font_bold": self.state.label_font_bold,
                "font_italic": self.state.label_font_italic,
                "font_underline": self.state.label_font_underline,
                "canvas_opacity": self.state.canvas_opacity,
                "name_advanced_enabled": bool(getattr(self, "name_advanced_enabled", True)),
                "name_rdkit_isolated": bool(getattr(self, "name_rdkit_isolated", True)),
            },
            "annotations": {
                "arrows": [],
                "brackets": [],
                "energy_diagrams": [],
                "semantic_diagrams": [],
                "orbitals": [],
                "images": [],
                "text_items": [],
                "wavy_anchors": [],
                "plates": []
            },
            "label_anchors": {
                str(atom_id): anchor for atom_id, anchor in self._group_anchor_overrides.items()
            },
        }

        for item in self.scene.items():
            parent_item = item.parentItem()
            if isinstance(parent_item, CompositeDiagramItem):
                continue
            if isinstance(item, ArrowItem):
                if isinstance(item, PreviewArrowItem):
                    continue
                data["annotations"]["arrows"].append({
                    "kind": item.kind(),
                    "start": {"x": item.start_point().x(), "y": item.start_point().y()},
                    "end": {"x": item.end_point().x(), "y": item.end_point().y()},
                    "curve_factor": item.curve_factor(),
                    "stroke_px": item.stroke_px(),
                    "opacity": self.item_raw_opacity(item),
                })
            elif isinstance(item, BracketItem):
                data["annotations"]["brackets"].append({
                    "kind": item._kind,
                    "rect": {
                        "x": item._base_rect.x(),
                        "y": item._base_rect.y(),
                        "w": item._base_rect.width(),
                        "h": item._base_rect.height()
                    },
                    "padding": item._padding,
                    "stroke_px": item.stroke_px(),
                    "opacity": self.item_raw_opacity(item),
                })
            elif isinstance(item, TextAnnotationItem):
                if bool(item.data(NUMBERING_TEXT_ROLE)):
                    continue
                data["annotations"]["text_items"].append({
                    "text": item.toPlainText(),
                    "html": item.toHtml(), # Store HTML for rich text (formatting, colors)
                    "x": item.pos().x(),
                    "y": item.pos().y(),
                    "rotation": item.rotation(),
                    "font": item.font().toString(),
                    "color": item.defaultTextColor().name(),
                    "text_width": item.textWidth(),
                    "opacity": self.item_raw_opacity(item),
                })
            elif isinstance(item, ImageAnnotationItem):
                rect = item.display_rect()
                data["annotations"]["images"].append({
                    "x": rect.x(),
                    "y": rect.y(),
                    "width": rect.width(),
                    "height": rect.height(),
                    "rotation": item.rotation(),
                    "z": item.zValue(),
                    "mime_type": item.mime_type(),
                    "data_b64": base64.b64encode(item.data_bytes()).decode("ascii"),
                    "source_name": item.source_name(),
                    "opacity": self.item_raw_opacity(item),
                })
            elif isinstance(item, (TLCPlateItem, GelElectrophoresisItem)):
                data["annotations"]["plates"].append(item.to_dict())
            elif isinstance(item, EnergyDiagramItem):
                rect = item.display_rect()
                data["annotations"]["energy_diagrams"].append({
                    "x": rect.x(),
                    "y": rect.y(),
                    "width": rect.width(),
                    "height": rect.height(),
                    "rotation": item.rotation(),
                    "z": item.zValue(),
                    "kind": item.kind(),
                    "slot_count": item.slot_count(),
                    "label": item.label(),
                    "label_side": item.label_side(),
                    "occupancies": list(item.occupancies()),
                    "style_payload": item.style_payload(),
                    "metadata": item.metadata(),
                    "opacity": self.item_raw_opacity(item),
                })
            elif isinstance(item, CompositeDiagramItem):
                payload = item.to_json()
                payload["opacity"] = self.item_raw_opacity(item)
                data["annotations"]["semantic_diagrams"].append(payload)
            elif isinstance(item, OrbitalAnnotationItem):
                anchor0 = item.anchor0()
                anchor1 = item.anchor1()
                data["annotations"]["orbitals"].append({
                    "kind": item.kind(),
                    "anchor0": {"x": anchor0.x(), "y": anchor0.y()},
                    "anchor1": {"x": anchor1.x(), "y": anchor1.y()},
                    "stroke_shaded_lobes": item.stroke_shaded_lobes(),
                    "part_styles": item.part_styles(),
                    "z": item.zValue(),
                    "opacity": self.item_raw_opacity(item),
                })
            elif isinstance(item, WavyAnchorItem):
                anchor_id = item.data(WAVY_ANCHOR_ROLE)
                angle = item.data(WAVY_ANCHOR_ANGLE_ROLE)
                length = item.data(WAVY_ANCHOR_LENGTH_ROLE)
                bond_id = item.data(WAVY_ANCHOR_BOND_ROLE)
                if anchor_id is None:
                    continue
                data["annotations"]["wavy_anchors"].append({
                    "anchor_id": int(anchor_id),
                    "angle": float(angle) if angle is not None else 0.0,
                    "length": float(length) if length is not None else self._wavy_anchor_length(),
                    "bond_id": int(bond_id) if bond_id is not None else None,
                    "opacity": self.item_raw_opacity(item),
                })
        return data

    def load_persistence_data(self, data: dict) -> None:
        """Restores canvas settings and non-structural items from a dictionary."""
        settings = data.get("settings", {})
        self.paper_width = settings.get("paper_width", DEFAULT_PAPER_WIDTH)
        self.paper_height = settings.get("paper_height", DEFAULT_PAPER_HEIGHT)
        self.set_show_grid(settings.get("show_grid", False))
        self.set_show_rulers(settings.get("show_rulers", False))
        self.state.bond_length = settings.get("bond_length", DEFAULT_BOND_LENGTH)
        self.state.use_aromatic_circles = settings.get("use_aromatic_circles", False)
        self.state.show_implicit_carbons = settings.get("show_carbons", False)
        self.state.show_implicit_hydrogens = settings.get("show_hydrogens", False)
        self.state.use_element_colors = settings.get("use_element_colors", True)
        self.state.label_font_family = settings.get("font_family", "Arial")
        self.state.label_font_size = settings.get("font_size", 11.0)
        self.state.label_font_bold = settings.get("font_bold", False)
        self.state.label_font_italic = settings.get("font_italic", False)
        self.state.label_font_underline = settings.get("font_underline", False)
        self.state.canvas_opacity = normalize_opacity(settings.get("canvas_opacity", 1.0))
        self.name_advanced_enabled = bool(settings.get("name_advanced_enabled", True))
        self.name_rdkit_isolated = bool(settings.get("name_rdkit_isolated", True))

        self._group_anchor_overrides = {
            int(atom_id): anchor
            for atom_id, anchor in data.get("label_anchors", {}).items()
            if anchor
        }

        # Re-apply paper size visual
        self._create_paper()

        # Restore Annotations
        annotations = data.get("annotations", {})
        
        for arrow_d in annotations.get("arrows", []):
            start = QPointF(arrow_d["start"]["x"], arrow_d["start"]["y"])
            end = QPointF(arrow_d["end"]["x"], arrow_d["end"]["y"])
            curve_factor = arrow_d.get("curve_factor")
            stroke_px = arrow_d.get("stroke_px")
            try:
                curve_factor_value = float(curve_factor) if curve_factor is not None else None
            except Exception:
                curve_factor_value = None
            try:
                stroke_px_value = float(stroke_px) if stroke_px is not None else None
            except Exception:
                stroke_px_value = None
            self.add_arrow_item(
                start,
                end,
                kind=arrow_d["kind"],
                curve_factor=curve_factor_value,
                stroke_px=stroke_px_value,
                opacity=arrow_d.get("opacity"),
            )

        for br_d in annotations.get("brackets", []):
            rect_d = br_d["rect"]
            rect = QRectF(rect_d["x"], rect_d["y"], rect_d["w"], rect_d["h"])
            kind = br_d.get("kind", "[]")
            padding = br_d.get("padding", 8.0)
            stroke_px = br_d.get("stroke_px")
            bracket = BracketItem(
                rect,
                kind=kind,
                padding=padding,
                stroke_px=stroke_px,
                style=self.drawing_style,
            )
            self.readd_bracket_item(
                bracket,
                rect,
                kind,
                padding=padding,
                stroke_px=stroke_px,
                opacity=br_d.get("opacity"),
            )

        for txt_d in annotations.get("text_items", []):
            text_item = TextAnnotationItem(txt_d["text"], txt_d["x"], txt_d["y"])
            if "html" in txt_d:
                text_item.setHtml(txt_d["html"])
            if "rotation" in txt_d:
                text_item.setRotation(txt_d["rotation"])
            if "font" in txt_d:
                font = QFont()
                font.fromString(txt_d["font"])
                text_item.setFont(font)
            if "color" in txt_d:
                text_item.setDefaultTextColor(QColor(txt_d["color"]))
            if "text_width" in txt_d:
                text_item.setTextWidth(float(txt_d["text_width"]))
            self.set_graphics_item_opacity(text_item, txt_d.get("opacity"))
            self.scene.addItem(text_item)

        for img_d in annotations.get("images", []):
            raw_b64 = img_d.get("data_b64")
            mime_type = img_d.get("mime_type", "image/png")
            if not raw_b64:
                continue
            try:
                image_data = base64.b64decode(raw_b64)
            except Exception:
                continue
            item = ImageAnnotationItem(
                image_data,
                mime_type,
                width=float(img_d.get("width", 1.0)),
                height=float(img_d.get("height", 1.0)),
                source_name=img_d.get("source_name"),
            )
            if not item.is_valid_image():
                continue
            item.set_display_rect(
                QRectF(
                    float(img_d.get("x", 0.0)),
                    float(img_d.get("y", 0.0)),
                    float(img_d.get("width", 1.0)),
                    float(img_d.get("height", 1.0)),
                )
            )
            item.setRotation(float(img_d.get("rotation", 0.0)))
            item.setZValue(float(img_d.get("z", 8.0)))
            self.set_graphics_item_opacity(item, img_d.get("opacity"))
            self.readd_image_item(item)

        for energy_d in annotations.get("energy_diagrams", []):
            item = EnergyDiagramItem(
                energy_d.get("kind", DEFAULT_ENERGY_DIAGRAM_KIND),
                label=energy_d.get("label"),
                label_side=energy_d.get("label_side"),
                occupancies=energy_d.get("occupancies"),
                slot_count=energy_d.get("slot_count"),
                style_payload=energy_d.get("style_payload"),
                metadata=energy_d.get("metadata"),
                width=float(energy_d.get("width", 1.0)),
                height=float(energy_d.get("height", 1.0)),
            )
            item.set_display_rect(
                QRectF(
                    float(energy_d.get("x", 0.0)),
                    float(energy_d.get("y", 0.0)),
                    float(energy_d.get("width", 1.0)),
                    float(energy_d.get("height", 1.0)),
                )
            )
            item.setRotation(float(energy_d.get("rotation", 0.0)))
            item.setZValue(float(energy_d.get("z", 44.0)))
            self.set_graphics_item_opacity(item, energy_d.get("opacity"))
            self.readd_energy_diagram_item(item)

        for semantic_d in annotations.get("semantic_diagrams", []):
            try:
                item = CompositeDiagramItem.from_json(semantic_d)
            except Exception:
                continue
            self.set_graphics_item_opacity(item, semantic_d.get("opacity"))
            self.readd_semantic_diagram_item(item)

        for orbital_d in annotations.get("orbitals", []):
            try:
                anchor0_d = orbital_d.get("anchor0", {})
                anchor1_d = orbital_d.get("anchor1", {})
                item = OrbitalAnnotationItem(
                    orbital_d.get("kind", DEFAULT_ORBITAL_KIND),
                    QPointF(float(anchor0_d.get("x", 0.0)), float(anchor0_d.get("y", 0.0))),
                    QPointF(float(anchor1_d.get("x", 0.0)), float(anchor1_d.get("y", 0.0))),
                    stroke_shaded_lobes=orbital_d.get("stroke_shaded_lobes"),
                    part_styles=orbital_d.get("part_styles"),
                )
            except Exception:
                continue
            item.setZValue(float(orbital_d.get("z", 44.0)))
            self.set_graphics_item_opacity(item, orbital_d.get("opacity"))
            self.readd_orbital_item(item)

        for anchor_d in annotations.get("wavy_anchors", []):
            anchor_id = anchor_d.get("anchor_id")
            if anchor_id is None or anchor_id not in self.model.atoms:
                continue
            angle = float(anchor_d.get("angle", 0.0))
            length = float(anchor_d.get("length", self._wavy_anchor_length()))
            bond_id = anchor_d.get("bond_id")
            item = WavyAnchorItem(QPointF(0.0, 0.0), QPointF(1.0, 0.0), style=self.drawing_style)
            item.setData(WAVY_ANCHOR_ROLE, int(anchor_id))
            item.setData(WAVY_ANCHOR_ANGLE_ROLE, angle)
            item.setData(WAVY_ANCHOR_LENGTH_ROLE, length)
            if bond_id is not None:
                item.setData(WAVY_ANCHOR_BOND_ROLE, int(bond_id))
            self.set_graphics_item_opacity(item, anchor_d.get("opacity"))
            self.readd_wavy_anchor_item(item)
            
        for plate_d in annotations.get("plates", []):
            item = PlateItem.from_json(plate_d)
            if item:
                # 1. Add plate to scene. This automatically brings in all child items (lanes, labels, etc.)
                # and ensures the relative coordinate system is ready.
                self.readd_plate_item(item)
                # 2. load_dict handles creating and positioning child spots/bands at their saved relative coordinates.
                item.load_dict(plate_d, scene=self.scene)

        # Full refresh to update atom visibility and circles
        self.refresh_atom_visibility()
        self.refresh_aromatic_circles()
        self.refresh_numbering_opacity()

    def clear_canvas(self) -> None:
        """Método auxiliar para clear canvas.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._overlays_ready = False
        self._cancel_drag()
        self.undo_stack.blockSignals(True)
        
        # Nullify Python references to scene items before scene.clear()
        # so that subsequent recreation logic doesn't try to access deleted C++ objects.
        self.paper = None
        self._grid_minor_item = None
        self._grid_major_item = None
        self._selection_box = None
        self._selection_handle = None
        self._selection_move_handle = None
        self._select_preview_path = None
        self._select_preview_rect = None
        self._bracket_preview = None
        self._hover_atom_indicator = None
        self._hover_bond_indicator = None
        self._optimize_zone = None
        self._preview_bond_item = None
        self._preview_ring_item = None
        self._preview_chain_item = None
        self._preview_chain_label = None
        self._preview_arrow_item = None
        self._electron_dots.clear()
        self._wavy_anchors.clear()
        self._numbering_overlay_items.clear()
        self._numbering_atom_items.clear()
        self._numbering_structure_items.clear()
        self._numbering_atom_base_pos.clear()
        self._numbering_structure_base_pos.clear()
        self._atom_numbering.clear()
        self._structure_numbering.clear()

        self.scene.clear()
        self.model.clear()
        self.undo_stack.clear()
        self.undo_stack.blockSignals(False)
        self.atom_items.clear()
        self.bond_items.clear()
        self.aromatic_circles.clear()
        self.arrow_items.clear()
        self.bracket_items.clear()
        self.energy_diagram_items.clear()
        self.semantic_diagram_items.clear()
        self.orbital_items.clear()
        self.image_items.clear()
        self._electron_dots.clear()
        self._implicit_h_overlays.clear()
        self._group_anchor_overrides.clear()
        self._ring_centers.clear()
        self._next_ring_id = 1
        self.state.selected_atoms.clear()
        self.state.selected_bonds.clear()
        self.bond_anchor_id = None
        self.hovered_atom_id = None
        self.hovered_bond_id = None
        
        self._create_paper()
        self._create_overlays()
        self.recompute_numbering()

    def _rebuild_items_from_model(self) -> None:
        """Create AtomItems and BondItems for everything currently in self.model."""
        self.atom_items.clear()
        self.bond_items.clear()
        self._suspend_numbering_refresh = True
        try:
            # 1. Create AtomItems
            for atom in self.model.atoms.values():
                item = AtomItem(
                    atom,
                    radius=ATOM_HIT_RADIUS,
                    show_carbon=self.state.show_implicit_carbons,
                    show_hydrogen=self.state.show_implicit_hydrogens,
                    label_font=self._label_font_for_atom(atom.id),
                    style=self.drawing_style,
                    use_element_colors=self.state.use_element_colors
                )
                item.setOpacity(self.effective_atom_opacity(atom))
                self.scene.addItem(item)
                self.atom_items[atom.id] = item

            # 2. Recalculate ring centers for double bond offsets
            self.refresh_ring_centers()

            # 3. Create BondItems with ring context
            ring_pairs = self._compute_ring_bond_pairs()
            for bond in self.model.bonds.values():
                a1 = self.model.atoms.get(bond.a1_id)
                a2 = self.model.atoms.get(bond.a2_id)
                if a1 and a2:
                    item = BondItem(
                        bond, a1, a2,
                        render_aromatic_as_circle=self.state.use_aromatic_circles,
                        style=self.drawing_style
                    )
                    ring_center = self._aromatic_ring_center_for_bond(bond) if bond.is_aromatic else None
                    if ring_center is None and bond.ring_id is not None:
                        ring_center = self._ring_centers.get(bond.ring_id)
                    item.set_ring_context(ring_center)
                    item.set_bond_in_ring(self._bond_in_ring_for_pairs(bond, ring_pairs))
                    item.set_offset_sign(self._bond_offset_sign(bond))
                    self._configure_bond_rendering(bond, item)
                    self._set_bond_item_join_context(bond, item)
                    trim_start = self._bond_endpoint_trim(bond, bond.a1_id)
                    trim_end = self._bond_endpoint_trim(bond, bond.a2_id)
                    item.set_endpoint_trim(trim_start, trim_end)
                    extend_start = self._bond_endpoint_extend(bond, bond.a1_id)
                    extend_end = self._bond_endpoint_extend(bond, bond.a2_id)
                    item.set_endpoint_extend(extend_start, extend_end)
                    # Explicitly update positions now that context and sign are set
                    item.update_positions(a1, a2)
                    item.setOpacity(self.effective_bond_opacity(bond))
                    self.scene.addItem(item)
                    self.bond_items[bond.id] = item

            # 4. Global refresh: labels, visibility, shrinks, and implicit H overlays
            self.refresh_atom_visibility()
            self.refresh_aromatic_circles()
            self.validate_structure()
        finally:
            self._suspend_numbering_refresh = False
        self.recompute_numbering()
        self.refresh_numbering_opacity()

    def refresh_ring_centers(self) -> None:
        """Recalculate geometric centers for all rings in the model."""
        self._ring_centers.clear()
        ring_to_atoms = {}
        for bond in self.model.bonds.values():
            if bond.ring_id is not None:
                atoms = ring_to_atoms.setdefault(bond.ring_id, set())
                atoms.add(bond.a1_id)
                atoms.add(bond.a2_id)
        
        for ring_id, atom_ids in ring_to_atoms.items():
            sum_x = 0.0
            sum_y = 0.0
            count = 0
            for aid in atom_ids:
                atom = self.model.get_atom(aid)
                if atom:
                    sum_x += atom.x
                    sum_y += atom.y
                    count += 1
            if count > 0:
                self.register_ring_center(ring_id, (sum_x / count, sum_y / count))

    def _compute_ring_bond_pairs(self) -> set[frozenset[int]]:
        """Método auxiliar para  compute ring bond pairs.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        view = MolView(self.model)
        rings = find_rings_simple(view)
        if not rings:
            return set()
        pairs: set[frozenset[int]] = set()
        for ring in rings:
            pairs.update(ring_bonds(view, ring))
        return pairs

    @staticmethod
    def _bond_in_ring_for_pairs(bond: Bond, ring_pairs: set[frozenset[int]]) -> bool:
        """Método auxiliar para  bond in ring for pairs.

        Args:
            bond: Descripción del parámetro.
            ring_pairs: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if bond.ring_id is not None:
            return True
        return frozenset({bond.a1_id, bond.a2_id}) in ring_pairs

    def _refresh_bond_ring_flags(
        self,
        ring_pairs: Optional[set[frozenset[int]]] = None,
    ) -> None:
        """Método auxiliar para  refresh bond ring flags.

        Args:
            ring_pairs: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if ring_pairs is None:
            ring_pairs = self._compute_ring_bond_pairs()
        for bond in self.model.bonds.values():
            item = self.bond_items.get(bond.id)
            if item is None:
                continue
            item.set_bond_in_ring(self._bond_in_ring_for_pairs(bond, ring_pairs))
            if item.style in (BondStyle.WEDGE, BondStyle.HASHED):
                atom1 = self.model.get_atom(bond.a1_id)
                atom2 = self.model.get_atom(bond.a2_id)
                if atom1 and atom2:
                    item.update_positions(atom1, atom2)

    def _insert_molgraph(self, graph: MolGraph, *, select_inserted: bool = False) -> None:
        """Método auxiliar para  insert molgraph.

        Args:
            graph: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if not graph.atoms:
            return
        xs = [atom.x for atom in graph.atoms.values()]
        ys = [atom.y for atom in graph.atoms.values()]
        min_x, max_x = min(xs), max(xs)
        min_y, max_y = min(ys), max(ys)
        center_x = (min_x + max_x) / 2
        center_y = (min_y + max_y) / 2
        target_x = self.paper_width / 2
        target_y = self.paper_height / 2
        dx = target_x - center_x
        dy = target_y - center_y

        self.begin_validation_batch()
        self.undo_stack.beginMacro("Paste molecule")
        inserted_atom_ids: list[int] = []
        try:
            id_map: Dict[int, int] = {}
            inserted_pairs: set[frozenset[int]] = set()
            for atom in graph.atoms.values():
                cmd = AddAtomCommand(
                    self.model,
                    self,
                    atom.element,
                    atom.x + dx,
                    atom.y + dy,
                    is_explicit=atom.is_explicit,
                    charge=atom.charge,
                    isotope=getattr(atom, "isotope", None),
                    radical_electrons=int(getattr(atom, "radical_electrons", 0) or 0),
                    oxidation_state=getattr(atom, "oxidation_state", None),
                    no_implicit=bool(getattr(atom, "no_implicit", False)),
                    auto_hydrogens=False,
                    is_coordination_center=getattr(atom, "is_coordination_center", False),
                    sphere_radius=getattr(atom, "sphere_radius", None),
                    sphere_color=getattr(atom, "sphere_color", None),
                    sphere_filled=bool(getattr(atom, "sphere_filled", True)),
                    sphere_transparent=bool(getattr(atom, "sphere_transparent", False)),
                )
                self.undo_stack.push(cmd)
                if cmd.atom_id is not None:
                    id_map[atom.id] = cmd.atom_id
                    inserted_atom_ids.append(cmd.atom_id)
            for bond in graph.bonds.values():
                a1 = id_map.get(bond.a1_id)
                a2 = id_map.get(bond.a2_id)
                if a1 is None or a2 is None:
                    continue
                pair = frozenset({a1, a2})
                if pair in inserted_pairs:
                    continue
                inserted_pairs.add(pair)
                if self.model.find_bond_between(a1, a2) is not None:
                    continue
                cmd = AddBondCommand(
                    self.model,
                    self,
                    a1,
                    a2,
                    bond.order,
                    bond.style,
                    bond.stereo,
                    is_aromatic=bond.is_aromatic,
                    stroke_px=bond.stroke_px,
                    color=bond.color,
                    donor_atom_id=id_map.get(getattr(bond, "donor_atom_id", -1)),
                    flex_curve_1=getattr(bond, "flex_curve_1", None),
                    flex_curve_2=getattr(bond, "flex_curve_2", None),
                )
                self.undo_stack.push(cmd)
        finally:
            self.undo_stack.endMacro()
            self.end_validation_batch()
        if any(bond.is_aromatic for bond in self.model.bonds.values()):
            self._kekulize_aromatic_bonds()
        if select_inserted:
            self._select_inserted_items(inserted_atom_ids)

    def _pick_template_attachment_atom_id(self, graph: MolGraph) -> Optional[int]:
        """Elige un átomo de la plantilla para conectar al átomo ancla.

        Se priorizan vértices terminales (grado mínimo) y, entre ellos, el más
        alejado del centro para favorecer una orientación externa.
        """
        if not graph.atoms:
            return None
        degrees: Dict[int, int] = {atom_id: 0 for atom_id in graph.atoms.keys()}
        for bond in graph.bonds.values():
            if bond.a1_id in degrees:
                degrees[bond.a1_id] += 1
            if bond.a2_id in degrees:
                degrees[bond.a2_id] += 1
        min_degree = min(degrees.values()) if degrees else 0
        candidate_ids = [atom_id for atom_id, degree in degrees.items() if degree == min_degree]
        if len(candidate_ids) == 1:
            return candidate_ids[0]
        xs = [atom.x for atom in graph.atoms.values()]
        ys = [atom.y for atom in graph.atoms.values()]
        center_x = (min(xs) + max(xs)) / 2.0
        center_y = (min(ys) + max(ys)) / 2.0
        return max(
            candidate_ids,
            key=lambda atom_id: (
                (graph.get_atom(atom_id).x - center_x) ** 2
                + (graph.get_atom(atom_id).y - center_y) ** 2
            ),
        )

    def _insert_molgraph_at(
        self,
        graph: MolGraph,
        target: QPointF,
        attach_to_atom_id: Optional[int] = None,
    ) -> None:
        """Método auxiliar para  insert molgraph at.

        Args:
            graph: Descripción del parámetro.
            target: Descripción del parámetro.
            attach_to_atom_id: Átomo del lienzo al que se conectará la plantilla.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if not graph.atoms:
            return
        xs = [atom.x for atom in graph.atoms.values()]
        ys = [atom.y for atom in graph.atoms.values()]
        min_x, max_x = min(xs), max(xs)
        min_y, max_y = min(ys), max(ys)
        center_x = (min_x + max_x) / 2
        center_y = (min_y + max_y) / 2
        attach_template_id = None
        if attach_to_atom_id is not None and attach_to_atom_id in self.model.atoms:
            attach_template_id = self._pick_template_attachment_atom_id(graph)
        if attach_template_id is not None:
            attach_atom = graph.get_atom(attach_template_id)
            anchor = self.model.get_atom(attach_to_atom_id)
            anchor_pos = QPointF(anchor.x, anchor.y)
            attach_target = self._compute_default_bond_endpoint(anchor_pos, attach_to_atom_id)
            dx = attach_target.x() - attach_atom.x
            dy = attach_target.y() - attach_atom.y
        else:
            dx = target.x() - center_x
            dy = target.y() - center_y

        self.begin_validation_batch()
        self.undo_stack.beginMacro("Paste molecule")
        try:
            id_map: Dict[int, int] = {}
            inserted_pairs: set[frozenset[int]] = set()
            for atom in graph.atoms.values():
                cmd = AddAtomCommand(
                    self.model,
                    self,
                    atom.element,
                    atom.x + dx,
                    atom.y + dy,
                    is_explicit=atom.is_explicit,
                    charge=atom.charge,
                    isotope=getattr(atom, "isotope", None),
                    radical_electrons=int(getattr(atom, "radical_electrons", 0) or 0),
                    oxidation_state=getattr(atom, "oxidation_state", None),
                    no_implicit=bool(getattr(atom, "no_implicit", False)),
                    auto_hydrogens=False,
                    is_coordination_center=getattr(atom, "is_coordination_center", False),
                    sphere_radius=getattr(atom, "sphere_radius", None),
                    sphere_color=getattr(atom, "sphere_color", None),
                    sphere_filled=bool(getattr(atom, "sphere_filled", True)),
                    sphere_transparent=bool(getattr(atom, "sphere_transparent", False)),
                )
                self.undo_stack.push(cmd)
                if cmd.atom_id is not None:
                    id_map[atom.id] = cmd.atom_id
            for bond in graph.bonds.values():
                a1 = id_map.get(bond.a1_id)
                a2 = id_map.get(bond.a2_id)
                if a1 is None or a2 is None:
                    continue
                pair = frozenset({a1, a2})
                if pair in inserted_pairs:
                    continue
                inserted_pairs.add(pair)
                if self.model.find_bond_between(a1, a2) is not None:
                    continue
                cmd = AddBondCommand(
                    self.model,
                    self,
                    a1,
                    a2,
                    bond.order,
                    bond.style,
                    bond.stereo,
                    is_aromatic=bond.is_aromatic,
                    stroke_px=bond.stroke_px,
                    color=bond.color,
                    donor_atom_id=id_map.get(getattr(bond, "donor_atom_id", -1)),
                    flex_curve_1=getattr(bond, "flex_curve_1", None),
                    flex_curve_2=getattr(bond, "flex_curve_2", None),
                )
                self.undo_stack.push(cmd)
            if attach_to_atom_id is not None and attach_template_id is not None:
                template_new_id = id_map.get(attach_template_id)
                if template_new_id is not None:
                    if self.model.find_bond_between(attach_to_atom_id, template_new_id) is None:
                        cmd = AddBondCommand(
                            self.model,
                            self,
                            attach_to_atom_id,
                            template_new_id,
                            1,
                            BondStyle.PLAIN,
                            BondStereo.NONE,
                            is_aromatic=False,
                        )
                        self.undo_stack.push(cmd)
        finally:
            self.undo_stack.endMacro()
            self.end_validation_batch()
        if any(bond.is_aromatic for bond in self.model.bonds.values()):
            self._kekulize_aromatic_bonds()

    def begin_template_insert_mode(self, graph: MolGraph, label: str) -> None:
        """Activa la inserción por clic para una plantilla."""
        self._pending_template_graph = graph
        self._pending_template_label = label
        self.setCursor(Qt.CursorShape.CrossCursor)

    def cancel_template_insert_mode(self) -> None:
        """Cancela la inserción pendiente de plantillas."""
        self._pending_template_graph = None
        self._pending_template_label = None
        if self._space_panning:
            self.setCursor(Qt.CursorShape.OpenHandCursor)
        else:
            self.unsetCursor()

    def _insert_ring_template(self, scene_pos: QPointF) -> None:
        """Método auxiliar para  insert ring template.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        template = self.state.active_ring_template
        if not template:
            return
        anomeric = (self.state.active_ring_anomeric or "beta").lower()
        anomeric_up = anomeric != "alpha"
        graph = None
        if template == "haworth":
            graph = build_haworth_template(
                self.state.bond_length, anomeric_up=anomeric_up, bold_front=True
            )
        elif template == "chair":
            graph = build_chair_template(
                self.state.bond_length, anomeric_up=anomeric_up, bold_front=True
            )
        if graph is None:
            return
        self._insert_molgraph_at(graph, scene_pos)

    def _insert_orbital_item(self, scene_pos: QPointF, kind: str | None = None) -> Optional[OrbitalAnnotationItem]:
        """Inserta un orbital vectorial persistente modelado por dos anchors."""
        orbital_kind = kind or getattr(self.state, "active_orbital_kind", DEFAULT_ORBITAL_KIND)
        extent = orbital_canvas_extent(orbital_kind, self.state.bond_length)
        item = OrbitalAnnotationItem(
            orbital_kind,
            QPointF(scene_pos),
            QPointF(scene_pos.x(), scene_pos.y() - extent),
        )
        self.undo_stack.push(AddOrbitalItemCommand(self, item))
        self._select_inserted_items(items=[item])
        self._refresh_scene_after_image_insert()
        return item

    def _insert_tlc_plate_item(self, scene_pos: QPointF) -> None:
        """Abre el dialogo e inserta una placa TLC."""
        dialog = TLCInsertDialog(self.window())
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return

        lanes = dialog.lanes()
        item = TLCPlateItem(lanes=lanes)
        item.setPos(scene_pos)
        self.scene.addItem(item)
        cmd = AddPlateItemCommand(self, item)
        self.undo_stack.push(cmd)
        item.add_spots_to_lanes(scene=self.scene)
        self.set_current_tool("tool_select")
        self.scene.clearSelection()
        item.setSelected(True)
        self._sync_selection_from_scene()

    def _insert_gel_electrophoresis_item(self, scene_pos: QPointF) -> None:
        """Abre el dialogo e inserta un gel de electroforesis."""
        dialog = GelInsertDialog(self.window())
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return

        lanes = dialog.lanes()
        item = GelElectrophoresisItem(lanes=lanes)
        item.setPos(scene_pos)
        self.scene.addItem(item)
        cmd = AddPlateItemCommand(self, item)
        self.undo_stack.push(cmd)
        item.add_bands_to_lanes(scene=self.scene)
        self.set_current_tool("tool_select")
        self.scene.clearSelection()
        item.setSelected(True)
        self._sync_selection_from_scene()

    def _insert_energy_diagram_item(
        self,
        scene_pos: QPointF,
        kind: str | None = None,
    ) -> Optional[EnergyDiagramItem]:
        """Inserta un diagrama de energia persistente centrado en `scene_pos`."""
        diagram_kind = kind or getattr(
            self.state,
            "active_energy_diagram_kind",
            DEFAULT_ENERGY_DIAGRAM_KIND,
        )
        size = energy_diagram_default_size(diagram_kind, self.state.bond_length)
        item = EnergyDiagramItem(
            diagram_kind,
            label=default_energy_label(diagram_kind),
            label_side=default_energy_label_side(diagram_kind),
            width=float(size.width()),
            height=float(size.height()),
        )
        item.set_display_rect(
            QRectF(
                float(scene_pos.x()) - float(size.width()) * 0.5,
                float(scene_pos.y()) - float(size.height()) * 0.5,
                float(size.width()),
                float(size.height()),
            )
        )
        self.undo_stack.push(AddEnergyDiagramItemCommand(self, item))
        self._select_inserted_items(items=[item])
        self._refresh_scene_after_image_insert()
        return item

    def insert_semantic_diagram(
        self,
        diagram: SemanticDiagram,
        scene_pos: QPointF,
    ) -> list[BaseItem]:
        """Inserta un diagrama electrónico semántico centrado cerca del click."""
        item = CompositeDiagramItem(diagram)
        rect = item.display_rect()
        item.set_display_rect(
            QRectF(
                float(scene_pos.x()) - rect.width() * 0.5,
                float(scene_pos.y()) - rect.height() * 0.5,
                rect.width(),
                rect.height(),
            )
        )
        self.undo_stack.push(AddCompositeDiagramItemCommand(self, item))
        self._select_inserted_items(items=[item])
        self._refresh_scene_after_image_insert()
        return [item]

    def insert_atomic_diagram(self, electron_count: int, scene_pos: QPointF) -> list[BaseItem]:
        """Construye e inserta un diagrama atómico semántico."""
        return self.insert_semantic_diagram(
            build_atomic_subshell_diagram(electron_count),
            scene_pos,
        )

    def insert_semantic_preset(self, preset_name: str, scene_pos: QPointF) -> list[BaseItem]:
        """Construye e inserta un preset químico semántico."""
        return self.insert_semantic_diagram(
            build_semantic_diagram_from_preset(preset_name),
            scene_pos,
        )

    def insert_mo_diatomic_diagram(
        self,
        left_label: str,
        right_label: str,
        total_electrons: int,
        ordering: str = "heavy_2p",
        include_core_1s: bool = False,
        title: str | None = None,
        scene_pos: QPointF | None = None,
    ) -> list[BaseItem]:
        """Construye e inserta un diagrama MO diatómico semántico."""
        target_pos = scene_pos or self.mapToScene(self.viewport().rect().center())
        return self.insert_semantic_diagram(
            build_diatomic_mo_diagram(
                left_label=left_label,
                right_label=right_label,
                total_electrons=total_electrons,
                ordering="light_2p" if ordering == "light_2p" else "heavy_2p",
                include_core_1s=include_core_1s,
                title=title,
            ),
            target_pos,
        )

    def insert_ligand_field_diagram(
        self,
        d_electrons: int,
        geometry: str = "octahedral",
        spin_mode: str = "high",
        title: str | None = None,
        scene_pos: QPointF | None = None,
    ) -> list[BaseItem]:
        """Construye e inserta un diagrama de campo ligando semántico."""
        target_pos = scene_pos or self.mapToScene(self.viewport().rect().center())
        return self.insert_semantic_diagram(
            build_ligand_field_diagram(
                d_electrons=d_electrons,
                geometry=geometry if geometry in {"octahedral", "tetrahedral", "square_planar"} else "octahedral",
                spin_mode="low" if spin_mode == "low" else "high",
                title=title,
            ),
            target_pos,
        )

    def add_atom_item(self, atom) -> None:
        """Añade atom item.

        Args:
            atom: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if atom.id in self.atom_items:
            return
        # Determine visibility based on state preferences
        show_c = self.state.show_implicit_carbons or atom.element != "C"
        show_h = self.state.show_implicit_hydrogens or atom.element != "H" or atom.is_explicit
        item = AtomItem(
            atom,
            show_carbon=show_c,
            show_hydrogen=show_h,
            label_font=self._label_font_for_atom(atom.id),
            style=self.drawing_style,
            use_element_colors=self.state.use_element_colors,
        )
        item.setOpacity(self.effective_atom_opacity(atom))
        self.scene.addItem(item)
        self.atom_items[atom.id] = item
        self._refresh_atom_label(atom.id)
        self.recompute_numbering()

    def update_atom_item(self, atom_id: int, x: float, y: float) -> None:
        """Actualiza atom item.

        Args:
            atom_id: Descripción del parámetro.
            x: Descripción del parámetro.
            y: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = self.atom_items.get(atom_id)
        if item is not None:
            atom_ref = self.model.atoms.get(atom_id)
            if atom_ref is not None and hasattr(item, "atom"):
                item.atom = atom_ref
            item.setPos(x, y)
            item.setOpacity(self.effective_atom_opacity(atom_id))
        self._update_electron_dots_for_atom(atom_id)
        self._update_wavy_anchors_for_atom(atom_id)

    def update_atom_item_element(
        self,
        atom_id: int,
        element: str,
        is_explicit: Optional[bool] = None,
    ) -> None:
        """Actualiza atom item element.

        Args:
            atom_id: Descripción del parámetro.
            element: Descripción del parámetro.
            is_explicit: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = self.atom_items.get(atom_id)
        if item is not None:
            atom_ref = self.model.atoms.get(atom_id)
            if atom_ref is not None and hasattr(item, "atom"):
                item.atom = atom_ref
            if is_explicit is None:
                is_explicit = self.model.get_atom(atom_id).is_explicit
            item.set_element(element, is_explicit=is_explicit)
            item.setOpacity(self.effective_atom_opacity(atom_id))
            self._refresh_atom_label(atom_id)
        self._update_electron_dots_for_atom(atom_id)
        self._update_wavy_anchors_for_atom(atom_id)

    def update_atom_item_charge(self, atom_id: int, charge: int) -> None:
        """Actualiza atom item charge.

        Args:
            atom_id: Descripción del parámetro.
            charge: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = self.atom_items.get(atom_id)
        if item is not None:
            atom_ref = self.model.atoms.get(atom_id)
            if atom_ref is not None and hasattr(item, "atom"):
                item.atom = atom_ref
            item.set_charge(charge)
            item.setOpacity(self.effective_atom_opacity(atom_id))
            self._refresh_atom_label(atom_id)
        self._sync_selection_from_scene()

    def _configure_bond_rendering(self, bond: Bond, item: BondItem) -> None:
        """Método auxiliar para  configure bond rendering.

        Args:
            bond: Descripción del parámetro.
            item: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        prefer_full_length = self.state.show_implicit_carbons and self.state.show_implicit_hydrogens
        is_aromatic_ring = bond.is_aromatic and bond.ring_id is not None
        a1 = self.model.get_atom(bond.a1_id)
        a2 = self.model.get_atom(bond.a2_id)
        has_hetero = False
        if a1 is not None and a2 is not None:
            has_hetero = (a1.element != "C") or (a2.element != "C")
        # Estilo base: vértice tipo C=C (línea sigma central + línea pi desplazada).
        # En heteroátomos la línea pi debe ser completa (sin recorte), pero
        # manteniendo la misma orientación del doble enlace tipo C=C.
        effective_order = bond.display_order if bond.display_order is not None else bond.order
        is_plain_double = (
            bond.style == BondStyle.PLAIN and effective_order == 2 and not is_aromatic_ring
        )
        def _neighbor_bonds(atom_id: int) -> list[Bond]:
            neighbors: list[Bond] = []
            for other in self.model.bonds.values():
                if other.id == bond.id:
                    continue
                if other.style == BondStyle.COORDINATION:
                    continue
                if other.a1_id != atom_id and other.a2_id != atom_id:
                    continue
                other_id = other.a2_id if other.a1_id == atom_id else other.a1_id
                other_atom = self.model.atoms.get(other_id)
                if other_atom is not None and other_atom.element == "H":
                    continue
                neighbors.append(other)
            return neighbors

        def _is_terminal_like(atom_id: int) -> bool:
            total = 0
            for other in self.model.bonds.values():
                if other.style == BondStyle.COORDINATION:
                    continue
                if other.a1_id != atom_id and other.a2_id != atom_id:
                    continue
                other_id = other.a2_id if other.a1_id == atom_id else other.a1_id
                other_atom = self.model.atoms.get(other_id)
                if other_atom is not None and other_atom.element == "H":
                    continue
                total += 1
            return total <= 1

        def _needs_centered_double() -> bool:
            if not is_plain_double or a1 is None or a2 is None:
                return False
            candidates = ((a1.id, a2.id), (a2.id, a1.id))
            for center_id, terminal_id in candidates:
                neighbors = _neighbor_bonds(center_id)
                if len(neighbors) != 2:
                    continue
                if any((nb.display_order if nb.display_order is not None else nb.order) != 1 for nb in neighbors):
                    continue
                if _is_terminal_like(terminal_id):
                    return True
            return False

        if has_hetero and is_plain_double:
            prefer_full_length = True
        symmetric_double = _needs_centered_double()
        item.set_multibond_rendering(prefer_full_length, symmetric_double)

    @staticmethod
    def _ray_ellipse_distance(rect: QRectF, ux: float, uy: float) -> Optional[float]:
        # Ray P = t * D, where D = (ux, uy). Origin (0,0).
        # Ellipse defined by rect: Center C, radii a, b.
        # Equation: ((x - cx)/a)^2 + ((y - cy)/b)^2 = 1
        # Substitute x = t*ux, y = t*uy
        
        """Método auxiliar para  ray ellipse distance.

        Args:
            rect: Descripción del parámetro.
            ux: Descripción del parámetro.
            uy: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        a = rect.width() / 2.0
        b = rect.height() / 2.0
        if a <= 1e-9 or b <= 1e-9:
            return None
            
        cx = rect.center().x()
        cy = rect.center().y()
        
        # Coefficients for At^2 + Bt + K = 0
        # (1/a^2)*(t*ux - cx)^2 + (1/b^2)*(t*uy - cy)^2 - 1 = 0
        
        inv_a2 = 1.0 / (a * a)
        inv_b2 = 1.0 / (b * b)
        
        A = inv_a2 * ux * ux + inv_b2 * uy * uy
        B = -2.0 * (inv_a2 * ux * cx + inv_b2 * uy * cy)
        K = inv_a2 * cx * cx + inv_b2 * cy * cy - 1.0
        
        # Discriminant
        disc = B * B - 4 * A * K
        if disc < 0:
            return None
            
        sqrt_disc = math.sqrt(disc)
        
        # Solutions
        t1 = (-B - sqrt_disc) / (2 * A)
        t2 = (-B + sqrt_disc) / (2 * A)
        
        # We want the smallest positive t that is > 0
        # Since ray starts at 0,0 which is usually inside, one t is neg, one is pos.
        # We want the positive one.
        
        best = None
        if t1 > 1e-6:
            best = t1
        if t2 > 1e-6:
            if best is None or t2 < best:
                best = t2
                
        return best

    @staticmethod
    def _ray_rect_distance(rect: QRectF, ux: float, uy: float) -> Optional[float]:
        """Devuelve la distancia positiva a la primera intersección rayo-rectángulo."""
        t_enter = float("-inf")
        t_exit = float("inf")
        for minimum, maximum, direction in (
            (float(rect.left()), float(rect.right()), float(ux)),
            (float(rect.top()), float(rect.bottom()), float(uy)),
        ):
            if abs(direction) <= 1e-9:
                if minimum <= 0.0 <= maximum:
                    continue
                return None
            t1 = minimum / direction
            t2 = maximum / direction
            if t1 > t2:
                t1, t2 = t2, t1
            t_enter = max(t_enter, t1)
            t_exit = min(t_exit, t2)
            if t_exit < t_enter:
                return None
        if t_exit <= 1e-6:
            return None
        if t_enter > 1e-6:
            return t_enter
        return t_exit

    def _hydrogen_label_side_rect(
        self,
        atom_id: int,
        element_text: str,
        h_text: str,
        *,
        prefix: bool,
    ) -> Optional[QRectF]:
        """Rectángulo aproximado ocupado por el bloque H/H2 a un lado del átomo."""
        if not element_text or not h_text:
            return None
        font = self._label_font_for_atom(atom_id)
        metrics = QFontMetrics(font)
        element_width = float(metrics.horizontalAdvance(element_text))
        hydrogen_width = float(metrics.horizontalAdvance(h_text))
        text_height = float(metrics.height())
        if element_width <= 0.0 or hydrogen_width <= 0.0 or text_height <= 0.0:
            return None

        outer_pad = max(2.0, float(self.drawing_style.stroke_px) * 1.6)
        inner_pad = min(1.0, hydrogen_width * 0.12)
        vertical_pad = max(1.5, float(self.drawing_style.stroke_px) * 1.2)
        half_element = element_width * 0.5

        if prefix:
            right = -half_element + inner_pad
            left = right - hydrogen_width - outer_pad
        else:
            left = half_element - inner_pad
            right = left + hydrogen_width + outer_pad

        top = -text_height * 0.5 - vertical_pad
        bottom = text_height * 0.5 + vertical_pad
        return QRectF(QPointF(left, top), QPointF(right, bottom)).normalized()

    def _hydrogen_side_obstruction_score(
        self,
        atom_id: int,
        element_text: str,
        h_count: int,
        *,
        prefix: bool,
    ) -> float:
        """Penaliza orientaciones donde el bloque H invade enlaces vecinos."""
        atom = self.model.atoms.get(atom_id)
        if atom is None or h_count <= 0:
            return 0.0

        h_text = "H" if h_count == 1 else f"H{int(h_count)}"
        h_rect = self._hydrogen_label_side_rect(atom_id, element_text, h_text, prefix=prefix)
        if h_rect is None:
            return 0.0

        score = 0.0
        side_sign = -1.0 if prefix else 1.0
        for bond in self.model.bonds.values():
            if bond.a1_id == atom_id:
                other = self.model.get_atom(bond.a2_id)
            elif bond.a2_id == atom_id:
                other = self.model.get_atom(bond.a1_id)
            else:
                continue
            if other is None or other.element == "H":
                continue
            dx = float(other.x - atom.x)
            dy = float(other.y - atom.y)
            length = math.hypot(dx, dy)
            if length <= 1e-6:
                continue
            ux = dx / length
            uy = dy / length
            hit_distance = self._ray_rect_distance(h_rect, ux, uy)
            if hit_distance is not None and hit_distance < length - 1e-6:
                proximity = max(0.0, 1.0 - hit_distance / max(length, 1.0))
                score += 100.0 + proximity * 25.0 + abs(ux) * 10.0
                continue

            # Penalización suave si el enlace apunta hacia el mismo lado del H
            # y se mantiene relativamente horizontal.
            same_side = max(0.0, side_sign * ux)
            horizontal = max(0.0, abs(ux) - abs(uy) * 0.35)
            score += same_side * horizontal * 4.0
        return score

    def _prefer_prefix_h(self, atom_id: int, element_text: str, h_count: int) -> bool:
        """Elige H antes/después del átomo minimizando obstrucción visual."""
        prefix_score = self._hydrogen_side_obstruction_score(
            atom_id,
            element_text,
            h_count,
            prefix=True,
        )
        suffix_score = self._hydrogen_side_obstruction_score(
            atom_id,
            element_text,
            h_count,
            prefix=False,
        )
        if abs(prefix_score - suffix_score) > 0.5:
            return prefix_score < suffix_score

        direction = self._label_open_direction(atom_id)
        return direction.x() < -0.2

    def remove_atom_item(self, atom_id: int) -> None:
        """Elimina atom item.

        Args:
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = self.atom_items.pop(atom_id, None)
        if item is not None:
            self.scene.removeItem(item)
        self._clear_implicit_h_overlays([atom_id])
        self._group_anchor_overrides.pop(atom_id, None)
        if self.bond_anchor_id == atom_id:
            self.bond_anchor_id = None
        if self.hovered_atom_id == atom_id:
            self.hovered_atom_id = None
        self.recompute_numbering()

    def add_bond_item(self, bond) -> None:
        """Añade bond item.

        Args:
            bond: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if bond.id in self.bond_items:
            return
        atom1 = self.model.get_atom(bond.a1_id)
        atom2 = self.model.get_atom(bond.a2_id)
        item = BondItem(
            bond,
            atom1,
            atom2,
            render_aromatic_as_circle=self.state.use_aromatic_circles,
            style=self.drawing_style,
        )
        ring_pairs = self._compute_ring_bond_pairs()
        ring_center = self._aromatic_ring_center_for_bond(bond) if bond.is_aromatic else None
        if ring_center is None and bond.ring_id is not None:
            ring_center = self._ring_centers.get(bond.ring_id)
        item.set_ring_context(ring_center)
        item.set_bond_in_ring(self._bond_in_ring_for_pairs(bond, ring_pairs))
        item.set_offset_sign(self._bond_offset_sign(bond))
        self._configure_bond_rendering(bond, item)
        self._set_bond_item_join_context(bond, item)
        trim_start = self._bond_endpoint_trim(bond, bond.a1_id)
        trim_end = self._bond_endpoint_trim(bond, bond.a2_id)
        item.set_endpoint_trim(trim_start, trim_end)
        extend_start = self._bond_endpoint_extend(bond, bond.a1_id)
        extend_end = self._bond_endpoint_extend(bond, bond.a2_id)
        item.set_endpoint_extend(extend_start, extend_end)
        item.update_positions(atom1, atom2)
        item.setOpacity(self.effective_bond_opacity(bond))
        self.scene.addItem(item)
        self.bond_items[bond.id] = item
        self._refresh_atom_label(bond.a1_id)
        self._refresh_atom_label(bond.a2_id)
        self._refresh_bond_ring_flags(ring_pairs)
        self.recompute_numbering()

    @staticmethod
    def _split_bracket_kind(kind: str) -> Optional[tuple[str, str]]:
        """Método auxiliar para  split bracket kind.

        Args:
            kind: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if kind in {"()", "[]", "{}"} and len(kind) == 2:
            return kind[0], kind[1]
        return None

    def _ensure_bracket_preview(self) -> QGraphicsRectItem:
        """Método auxiliar para  ensure bracket preview.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._bracket_preview is None:
            preview = QGraphicsRectItem()
            pen = QPen(QColor("#4A90D9"), 1.0, Qt.PenStyle.DashLine)
            preview.setPen(pen)
            preview.setBrush(QBrush(Qt.BrushStyle.NoBrush))
            preview.setZValue(40)
            self.scene.addItem(preview)
            self._bracket_preview = preview
        self._bracket_preview.setVisible(True)
        return self._bracket_preview

    def _clear_bracket_preview(self) -> None:
        """Método auxiliar para  clear bracket preview.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._bracket_preview is not None:
            self.scene.removeItem(self._bracket_preview)
            self._bracket_preview = None
        self._bracket_drag_start = None

    def update_bond_item(self, bond_id: int) -> None:
        """Actualiza bond item.

        Args:
            bond_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        bond = self.model.get_bond(bond_id)
        atom1 = self.model.get_atom(bond.a1_id)
        atom2 = self.model.get_atom(bond.a2_id)
        item = self.bond_items.get(bond_id)
        if item is not None:
            ring_pairs = self._compute_ring_bond_pairs()
            ring_center = self._aromatic_ring_center_for_bond(bond) if bond.is_aromatic else None
            if ring_center is None and bond.ring_id is not None:
                ring_center = self._ring_centers.get(bond.ring_id)
            item.set_ring_context(ring_center)
            item.set_bond_in_ring(self._bond_in_ring_for_pairs(bond, ring_pairs))
            item.set_offset_sign(self._bond_offset_sign(bond))
            self._configure_bond_rendering(bond, item)
            self._set_bond_item_join_context(bond, item)
            trim_start = self._bond_endpoint_trim(bond, bond.a1_id)
            trim_end = self._bond_endpoint_trim(bond, bond.a2_id)
            item.set_endpoint_trim(trim_start, trim_end)
            extend_start = self._bond_endpoint_extend(bond, bond.a1_id)
            extend_end = self._bond_endpoint_extend(bond, bond.a2_id)
            item.set_endpoint_extend(extend_start, extend_end)
            item.set_bond(bond, atom1, atom2)
            item.set_style(self.drawing_style, atom1, atom2)
            item.setOpacity(self.effective_bond_opacity(bond))
        self._refresh_bond_ring_flags(ring_pairs if item is not None else None)
        self._refresh_atom_label(bond.a1_id)
        self._refresh_atom_label(bond.a2_id)

    def _bond_render_width(self, bond: Bond) -> float:
        """Método auxiliar para  bond render width.

        Args:
            bond: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        base = bond.stroke_px if bond.stroke_px is not None else self.drawing_style.stroke_px
        if bond.style == BondStyle.BOLD:
            return max(base * 2.2, base + 1.0)
        if bond.style == BondStyle.WEDGE:
            scale = base / self.drawing_style.stroke_px if self.drawing_style.stroke_px > 1e-6 else 1.0
            width = self.drawing_style.wedge_width_px * (0.72 + 0.28 * math.sqrt(max(scale, 1e-6)))
            return max(width, base * 2.3)
        if bond.style == BondStyle.HASHED:
            scale = base / self.drawing_style.stroke_px if self.drawing_style.stroke_px > 1e-6 else 1.0
            return max(self.drawing_style.hash_stroke_px * scale, base * 0.85)
        return base

    def _bond_neighbor_vectors(
        self,
        bond: Bond,
        atom_id: int,
        *,
        simple_only: bool = False,
    ) -> list[tuple[float, float, float, float, float]]:
        """Devuelve vecinos con dirección, ancho y centro de extremo renderizado."""
        anchor = self.model.get_atom(atom_id)
        if anchor is None:
            return []
        neighbors: list[tuple[float, float, float, float, float]] = []
        for other in self.model.bonds.values():
            if other.id == bond.id:
                continue
            if simple_only and other.style not in {BondStyle.PLAIN, BondStyle.BOLD}:
                continue
            if other.a1_id == atom_id:
                other_atom = self.model.get_atom(other.a2_id)
            elif other.a2_id == atom_id:
                other_atom = self.model.get_atom(other.a1_id)
            else:
                continue
            if other_atom is None:
                continue
            dx = other_atom.x - anchor.x
            dy = other_atom.y - anchor.y
            length = math.hypot(dx, dy)
            if length <= 1e-6:
                continue
            nux = dx / length
            nuy = dy / length
            nwidth = self._bond_render_width(other)
            shrink = self._label_shrink_for_atom(atom_id, nux, nuy)
            trim = self._bond_endpoint_trim(other, atom_id)
            advance = max(0.0, shrink + trim)
            end_x = anchor.x + nux * advance
            end_y = anchor.y + nuy * advance

            cap_extra = 0.0
            if self.drawing_style.cap_style in (
                Qt.PenCapStyle.RoundCap,
                Qt.PenCapStyle.SquareCap,
            ):
                cap_extra = nwidth * 0.5

            extend = 0.0
            if self.drawing_style.cap_style == Qt.PenCapStyle.FlatCap and advance <= 1e-6:
                extend = self._bond_endpoint_extend(other, atom_id)

            edge_cx = end_x - nux * (cap_extra + extend)
            edge_cy = end_y - nuy * (cap_extra + extend)
            neighbors.append((nux, nuy, nwidth, edge_cx, edge_cy))
        return neighbors

    def _set_bond_item_join_context(self, bond: Bond, item: BondItem) -> None:
        """Inyecta contexto de enlaces vecinos para geometrías dependientes de unión."""
        if bond.style != BondStyle.WEDGE:
            item.set_wedge_join_neighbors([], [])
            return
        start_neighbors = self._bond_neighbor_vectors(bond, bond.a1_id)
        # El extremo ancho de la cuña se integra contra enlaces planos:
        # - Plain (incluye dobles/aromáticos por order/display_order)
        # - Bold
        # Esto permite vértices limpios en uniones cuña->bold.
        end_neighbors = self._bond_neighbor_vectors(bond, bond.a2_id, simple_only=True)
        item.set_wedge_join_neighbors(start_neighbors, end_neighbors)

    def _bond_endpoint_extend(self, bond: Bond, atom_id: int) -> float:
        """Método auxiliar para  bond endpoint extend.

        Args:
            bond: Descripción del parámetro.
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self.drawing_style.cap_style != Qt.PenCapStyle.FlatCap:
            return 0.0
        if self._atom_degree(atom_id) != 2:
            return 0.0
        atom = self.model.get_atom(atom_id)
        if self._atom_label_visible(atom):
            return 0.0
        other_bond = None
        for candidate in self.model.bonds.values():
            if candidate.id == bond.id:
                continue
            if candidate.a1_id == atom_id or candidate.a2_id == atom_id:
                other_bond = candidate
                break
        if other_bond is None:
            return 0.0

        def unit_vector(b: Bond) -> tuple[float, float]:
            """Método auxiliar para unit vector.

            Args:
                b: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
            other_id = b.a2_id if b.a1_id == atom_id else b.a1_id
            other_atom = self.model.get_atom(other_id)
            dx = other_atom.x - atom.x
            dy = other_atom.y - atom.y
            length = math.hypot(dx, dy) or 1.0
            return dx / length, dy / length

        v1x, v1y = unit_vector(bond)
        v2x, v2y = unit_vector(other_bond)
        dot = max(-1.0, min(1.0, v1x * v2x + v1y * v2y))
        theta = math.acos(dot)
        if theta <= 1e-3:
            return 0.0
        width = max(self._bond_render_width(bond), self._bond_render_width(other_bond))
        extension = (width * 0.5) / math.tan(theta / 2.0)
        return max(0.0, min(extension, width * 2.0))

    def _atom_label_visible(self, atom) -> bool:
        """Método auxiliar para  atom label visible.

        Args:
            atom: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if atom.is_explicit:
            return True
        if atom.element == "C":
            return self.state.show_implicit_carbons
        if atom.element == "H":
            return self.state.show_implicit_hydrogens
        return True

    def _bond_endpoint_trim(self, bond: Bond, atom_id: int) -> float:
        """Método auxiliar para  bond endpoint trim.

        Args:
            bond: Descripción del parámetro.
            atom_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        # Keep stereo wedges visually full-length at junctions (ChemDraw-like).
        # Label collision avoidance is handled separately via label_shrink.
        if bond.style in (BondStyle.WEDGE, BondStyle.HASHED, BondStyle.COORDINATION):
            return 0.0
        if self._atom_degree(atom_id) < 2:
            return 0.0
        width = self._bond_render_width(bond)
        max_other = 0.0
        for other in self.model.bonds.values():
            if other.id == bond.id:
                continue
            if other.a1_id == atom_id or other.a2_id == atom_id:
                max_other = max(max_other, self._bond_render_width(other))
        if max_other <= 0.0:
            return 0.0
        if width <= max_other:
            return 0.0
        return max(0.0, (width - max_other) * 0.5)

    def update_bond_items_for_atoms(self, atom_ids: set[int]) -> None:
        """Actualiza bond items for atoms.

        Args:
            atom_ids: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._clear_trackball_reference_if_desynced(atom_ids)
        for bond in self.model.bonds.values():
            if bond.a1_id in atom_ids or bond.a2_id in atom_ids:
                self.update_bond_item(bond.id)
        self._refresh_implicit_h_overlays(atom_ids)
        self._update_aromatic_circles_for_atoms(atom_ids)
        self.recompute_numbering()

    def remove_bond_item(self, bond_id: int) -> None:
        """Elimina bond item.

        Args:
            bond_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = self.bond_items.pop(bond_id, None)
        if item is not None:
            a1_id = item.a1_id
            a2_id = item.a2_id
            self.scene.removeItem(item)
            self._refresh_atom_label(a1_id)
            self._refresh_atom_label(a2_id)
        self._refresh_bond_ring_flags()
        self.recompute_numbering()

    def allocate_ring_id(self) -> int:
        """Método auxiliar para allocate ring id.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        ring_id = self._next_ring_id
        self._next_ring_id += 1
        return ring_id

    def register_ring_center(self, ring_id: int, center: tuple[float, float]) -> None:
        """Método auxiliar para register ring center.

        Args:
            ring_id: Descripción del parámetro.
            center: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._ring_centers[ring_id] = QPointF(center[0], center[1])

    def unregister_ring_center(self, ring_id: int) -> None:
        """Método auxiliar para unregister ring center.

        Args:
            ring_id: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._ring_centers.pop(ring_id, None)

    def _bond_offset_sign(self, bond) -> int:
        """Método auxiliar para  bond offset sign.

        Args:
            bond: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atom1 = self.model.get_atom(bond.a1_id)
        atom2 = self.model.get_atom(bond.a2_id)
        p1 = QPointF(atom1.x, atom1.y)
        p2 = QPointF(atom2.x, atom2.y)
        mid = QPointF((p1.x() + p2.x()) / 2, (p1.y() + p2.y()) / 2)
        dx = p2.x() - p1.x()
        dy = p2.y() - p1.y()
        length = math.hypot(dx, dy) or 1.0
        nx = -dy / length
        ny = dx / length

        if bond.is_aromatic:
            aromatic_center = self._aromatic_ring_center_for_bond(bond)
            if aromatic_center is not None:
                vx = aromatic_center.x() - mid.x()
                vy = aromatic_center.y() - mid.y()
                return 1 if (nx * vx + ny * vy) >= 0 else -1
            neighbor_points = []
            for other_bond in self.model.bonds.values():
                if not other_bond.is_aromatic or other_bond.id == bond.id:
                    continue
                if other_bond.a1_id == atom1.id and other_bond.a2_id != atom2.id:
                    other = self.model.get_atom(other_bond.a2_id)
                elif other_bond.a2_id == atom1.id and other_bond.a1_id != atom2.id:
                    other = self.model.get_atom(other_bond.a1_id)
                elif other_bond.a1_id == atom2.id and other_bond.a2_id != atom1.id:
                    other = self.model.get_atom(other_bond.a2_id)
                elif other_bond.a2_id == atom2.id and other_bond.a1_id != atom1.id:
                    other = self.model.get_atom(other_bond.a1_id)
                else:
                    continue
                neighbor_points.append(QPointF(other.x, other.y))
            if neighbor_points:
                avg_x = sum(p.x() for p in neighbor_points) / len(neighbor_points)
                avg_y = sum(p.y() for p in neighbor_points) / len(neighbor_points)
                vx = avg_x - mid.x()
                vy = avg_y - mid.y()
                return 1 if (nx * vx + ny * vy) >= 0 else -1

        def neighbor_score(atom_id: int) -> float:
            """Método auxiliar para neighbor score.

            Args:
                atom_id: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
            score = 0.0
            for other_bond in self.model.bonds.values():
                if other_bond.id == bond.id:
                    continue
                if other_bond.a1_id == atom_id:
                    other = self.model.get_atom(other_bond.a2_id)
                elif other_bond.a2_id == atom_id:
                    other = self.model.get_atom(other_bond.a1_id)
                else:
                    continue
                vx = other.x - mid.x()
                vy = other.y - mid.y()
                score += nx * vx + ny * vy
            return score

        score = neighbor_score(atom1.id) + neighbor_score(atom2.id)
        return -1 if score > 0 else 1

    def refresh_atom_visibility(self) -> None:
        """Update visibility of all atoms based on current state settings."""
        for atom_id, item in self.atom_items.items():
            atom = self.model.get_atom(atom_id)
            show_c = self.state.show_implicit_carbons or atom.element != "C"
            show_h = self.state.show_implicit_hydrogens or atom.element != "H" or atom.is_explicit
            item.set_visibility_flags(show_c, show_h)
            self._refresh_atom_label(atom_id)
        self._refresh_implicit_h_overlays()

    def clean_2d_fallback(self, atom_ids: Optional[set[int]] = None, iterations: int = 50) -> None:
        """Basic 2D cleanup without RDKit using length + angle relaxation."""
        if atom_ids is None:
            atom_ids = set(self.model.atoms.keys())
        if not atom_ids:
            return

        before = {
            aid: (self.model.get_atom(aid).x, self.model.get_atom(aid).y) for aid in atom_ids
        }
        before_center = self._center_for_atoms(list(atom_ids)) or QPointF(0.0, 0.0)

        positions = {aid: QPointF(x, y) for aid, (x, y) in before.items()}

        target_len = float(self.state.bond_length)
        angle_step_deg = 12.0
        iterations = max(1, int(iterations))

        neighbor_map: dict[int, list[int]] = {}
        bonds_in_selection = []
        bond_lookup: dict[tuple[int, int], Bond] = {}
        for bond in self.model.bonds.values():
            if bond.a1_id in atom_ids and bond.a2_id in atom_ids:
                neighbor_map.setdefault(bond.a1_id, []).append(bond.a2_id)
                neighbor_map.setdefault(bond.a2_id, []).append(bond.a1_id)
                bonds_in_selection.append(bond)
                key = (min(bond.a1_id, bond.a2_id), max(bond.a1_id, bond.a2_id))
                bond_lookup[key] = bond

        components: list[set[int]] = []
        visited: set[int] = set()
        for aid in atom_ids:
            if aid in visited:
                continue
            stack = [aid]
            comp = set()
            while stack:
                node = stack.pop()
                if node in visited:
                    continue
                visited.add(node)
                comp.add(node)
                for nb in neighbor_map.get(node, []):
                    if nb not in visited:
                        stack.append(nb)
            components.append(comp)

        def rotate_point(p: QPointF, center: QPointF, delta_deg: float) -> QPointF:
            """Método auxiliar para rotate point.

            Args:
                p: Descripción del parámetro.
                center: Descripción del parámetro.
                delta_deg: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
            rad = math.radians(delta_deg)
            dx = p.x() - center.x()
            dy = p.y() - center.y()
            cos_t = math.cos(rad)
            sin_t = math.sin(rad)
            return QPointF(center.x() + dx * cos_t - dy * sin_t, center.y() + dx * sin_t + dy * cos_t)

        def atom_geometry(atom_id: int) -> str:
            """Infer geometry at `atom_id` from incident bonds in selection."""
            has_triple = False
            has_double = False
            for bond in bonds_in_selection:
                if bond.a1_id != atom_id and bond.a2_id != atom_id:
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

        def non_h_degree(atom_id: int, component: set[int]) -> int:
            """Degree inside component ignoring hydrogen neighbors."""
            total = 0
            for neigh in neighbor_map.get(atom_id, []):
                if neigh not in component:
                    continue
                atom = self.model.get_atom(neigh)
                if atom is not None and atom.element == "H":
                    continue
                total += 1
            return total

        def choose_tree_root(component: set[int]) -> int:
            """Pick deterministic root favoring terminal heavy atoms."""
            heavy_atoms = [
                aid
                for aid in component
                if (self.model.get_atom(aid) is not None and self.model.get_atom(aid).element != "H")
            ]
            candidates = heavy_atoms if heavy_atoms else list(component)
            leaves = [aid for aid in candidates if non_h_degree(aid, component) <= 1]
            if leaves:
                candidates = leaves
            return min(candidates, key=lambda aid: (before[aid][0], before[aid][1], aid))

        def normalize_tree(component: set[int], roots: list[int]) -> None:
            """Método auxiliar para normalize tree.

            Args:
                component: Descripción del parámetro.
                roots: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
            queue = []
            parent_of: dict[int, Optional[int]] = {}
            placed = set()
            for root in roots:
                if root not in component:
                    continue
                parent_of[root] = None
                placed.add(root)
                queue.append(root)

            while queue:
                parent = queue.pop(0)
                p_parent = positions[parent]

                incoming_angle: Optional[float] = None
                parent_parent = parent_of.get(parent)
                if parent_parent is not None and parent_parent in placed:
                    incoming_angle = angle_deg(p_parent, positions[parent_parent])

                existing_angles: list[float] = []
                for neigh in neighbor_map.get(parent, []):
                    if neigh not in component or neigh not in placed:
                        continue
                    if neigh == parent_parent:
                        continue
                    existing_angles.append(angle_deg(p_parent, positions[neigh]))
                if incoming_angle is not None:
                    existing_angles = [incoming_angle] + existing_angles

                neighbors = [
                    neigh
                    for neigh in neighbor_map.get(parent, [])
                    if neigh in component and neigh not in placed
                ]
                if not neighbors:
                    continue

                def neighbor_sort_key(neigh: int) -> tuple[int, int, int]:
                    atom = self.model.get_atom(neigh)
                    is_h = 1 if (atom is not None and atom.element == "H") else 0
                    return (is_h, -non_h_degree(neigh, component), neigh)

                neighbors.sort(key=neighbor_sort_key)
                parent_atom = self.model.get_atom(parent)
                geometry = atom_geometry(parent)

                def neighbor_bond_length(neigh: int) -> float:
                    key = (min(parent, neigh), max(parent, neigh))
                    bond = bond_lookup.get(key)
                    length = target_len
                    if bond is not None and bond.order >= 2:
                        length = target_len * 0.98
                    return length

                def place_neighbor_with_angle(neigh: int, chosen_angle: float) -> None:
                    length = neighbor_bond_length(neigh)
                    positions[neigh] = endpoint_from_angle_len(p_parent, chosen_angle, length)
                    existing_angles.append(chosen_angle)
                    placed.add(neigh)
                    parent_of[neigh] = parent
                    queue.append(neigh)

                def candidate_collision_metrics(
                    neigh: int,
                    angle_value: float,
                ) -> tuple[int, float, float]:
                    p1 = endpoint_from_angle_len(p_parent, angle_value, neighbor_bond_length(neigh))
                    intersections = 0
                    min_atom_dist = float("inf")
                    min_bond_dist = float("inf")
                    atom_threshold = target_len * MIN_ATOM_DIST_SCALE
                    bond_threshold = target_len * MIN_BOND_DIST_SCALE

                    for atom_id in placed:
                        if atom_id == parent:
                            continue
                        atom_pos = positions.get(atom_id)
                        if atom_pos is None:
                            continue
                        dist = math.hypot(p1.x() - atom_pos.x(), p1.y() - atom_pos.y())
                        min_atom_dist = min(min_atom_dist, dist)
                        if dist < atom_threshold:
                            intersections += 1

                    for bond in bonds_in_selection:
                        if bond.a1_id not in placed or bond.a2_id not in placed:
                            continue
                        if parent in {bond.a1_id, bond.a2_id}:
                            continue
                        p_a = positions.get(bond.a1_id)
                        p_b = positions.get(bond.a2_id)
                        if p_a is None or p_b is None:
                            continue
                        if segments_intersect(p_parent, p1, p_a, p_b):
                            intersections += 1
                        bond_dist = segment_min_distance(p_parent, p1, p_a, p_b)
                        min_bond_dist = min(min_bond_dist, bond_dist)
                        if bond_dist < bond_threshold:
                            intersections += 1

                    return intersections, min_atom_dist, min_bond_dist

                def pick_angle_for_neighbor(neigh: int, mouse_angle: float) -> float:
                    if geometry == "sp3" and len(existing_angles) >= 3:
                        candidates = self._sp3_congested_directions_deg(existing_angles)
                    else:
                        candidates = candidate_directions_deg(geometry, existing_angles, mouse_angle)
                    occupied_tolerance = (
                        SP3_OCCUPIED_TOLERANCE_DEG
                        if geometry == "sp3" and len(existing_angles) >= 3
                        else ANGLE_OCCUPIED_TOLERANCE_DEG
                    )
                    candidates = filter_occupied_angles_deg(
                        candidates,
                        existing_angles,
                        occupied_tolerance,
                    )
                    if not candidates:
                        if geometry == "sp3" and len(existing_angles) >= 3:
                            candidates = self._sp3_congested_directions_deg(existing_angles)
                        if not candidates:
                            candidates = candidate_directions_deg(geometry, [], mouse_angle)
                    if not candidates:
                        candidates = [mouse_angle]

                    preferred: list[float] = []
                    if geometry == "sp3" and incoming_angle is not None:
                        sp3_angle = self._sp3_display_angle_deg()
                        preferred = [
                            (incoming_angle + sp3_angle) % 360.0,
                            (incoming_angle - sp3_angle) % 360.0,
                        ]

                    atom_threshold = target_len * MIN_ATOM_DIST_SCALE
                    bond_threshold = target_len * MIN_BOND_DIST_SCALE

                    def _candidate_score(angle_value: float) -> tuple[int, float]:
                        intersections, min_atom_dist, min_bond_dist = candidate_collision_metrics(
                            neigh,
                            angle_value,
                        )
                        score = angle_distance_deg(angle_value, mouse_angle)
                        if preferred:
                            if (
                                min(
                                    angle_distance_deg(angle_value, pref)
                                    for pref in preferred
                                )
                                <= 15.0
                            ):
                                score -= 15.0
                        score += intersections * 100.0
                        if min_atom_dist < atom_threshold:
                            score += (atom_threshold - min_atom_dist) * 5.0
                        if min_bond_dist < bond_threshold:
                            score += (bond_threshold - min_bond_dist) * 5.0
                        valid = (
                            intersections == 0
                            and min_atom_dist >= atom_threshold
                            and min_bond_dist >= bond_threshold
                        )
                        return (0 if valid else 1, score)

                    return min(candidates, key=_candidate_score)

                non_h_neighbors = []
                h_neighbors = []
                for neigh in neighbors:
                    atom = self.model.get_atom(neigh)
                    if atom is not None and atom.element == "H":
                        h_neighbors.append(neigh)
                    else:
                        non_h_neighbors.append(neigh)

                for neigh in non_h_neighbors:
                    parent_before = before.get(parent)
                    neigh_before = before.get(neigh)
                    if parent_before is not None and neigh_before is not None:
                        mouse_angle = angle_deg(
                            QPointF(parent_before[0], parent_before[1]),
                            QPointF(neigh_before[0], neigh_before[1]),
                        )
                    elif incoming_angle is not None:
                        mouse_angle = incoming_angle
                    else:
                        mouse_angle = 0.0
                    chosen = pick_angle_for_neighbor(neigh, mouse_angle)
                    place_neighbor_with_angle(neigh, chosen)

                if (
                    h_neighbors
                    and parent_atom is not None
                    and parent_atom.element == "C"
                    and geometry == "sp3"
                ):
                    target_h_angles = self._find_best_h_positions(existing_angles, len(h_neighbors))
                    if target_h_angles:
                        h_before_angles: list[float] = []
                        for neigh in h_neighbors:
                            parent_before = before.get(parent)
                            neigh_before = before.get(neigh)
                            if parent_before is not None and neigh_before is not None:
                                h_before_angles.append(
                                    angle_deg(
                                        QPointF(parent_before[0], parent_before[1]),
                                        QPointF(neigh_before[0], neigh_before[1]),
                                    )
                                )
                            elif incoming_angle is not None:
                                h_before_angles.append(incoming_angle)
                            else:
                                h_before_angles.append(0.0)

                        best_perm = tuple(range(len(target_h_angles)))
                        best_score = float("inf")
                        for perm in itertools.permutations(range(len(target_h_angles)), len(h_neighbors)):
                            score = 0.0
                            for idx, target_idx in enumerate(perm):
                                score += angle_distance_deg(h_before_angles[idx], target_h_angles[target_idx])
                            if score < best_score:
                                best_score = score
                                best_perm = perm
                        for idx, neigh in enumerate(h_neighbors):
                            target_idx = best_perm[idx]
                            place_neighbor_with_angle(neigh, target_h_angles[target_idx])
                        h_neighbors = []

                for neigh in h_neighbors:
                    parent_before = before.get(parent)
                    neigh_before = before.get(neigh)
                    if parent_before is not None and neigh_before is not None:
                        mouse_angle = angle_deg(
                            QPointF(parent_before[0], parent_before[1]),
                            QPointF(neigh_before[0], neigh_before[1]),
                        )
                    elif incoming_angle is not None:
                        mouse_angle = incoming_angle
                    else:
                        mouse_angle = 0.0
                    chosen = pick_angle_for_neighbor(neigh, mouse_angle)
                    place_neighbor_with_angle(neigh, chosen)

        def normalize_component_tree(component: set[int]) -> None:
            """Normalize entire acyclic component from a deterministic root."""
            if not component:
                return
            root = choose_tree_root(component)
            normalize_tree(component, [root])

        for component in components:
            if len(component) <= 1:
                continue
            # Find ring core by peeling leaves
            degrees = {aid: len([n for n in neighbor_map.get(aid, []) if n in component]) for aid in component}
            ring_atoms = set(component)
            leaf_queue = [aid for aid, deg in degrees.items() if deg <= 1]
            while leaf_queue:
                node = leaf_queue.pop()
                if node not in ring_atoms:
                    continue
                ring_atoms.remove(node)
                for nb in neighbor_map.get(node, []):
                    if nb in ring_atoms:
                        degrees[nb] -= 1
                        if degrees[nb] <= 1:
                            leaf_queue.append(nb)

            if not ring_atoms:
                # Acyclic: rebuild with ideal local geometry instead of preserving
                # incoming orthogonal/skewed angles from the drawing state.
                normalize_component_tree(component)
                continue

            # Relax ring core
            for _ in range(iterations):
                for bond in bonds_in_selection:
                    if bond.a1_id not in ring_atoms or bond.a2_id not in ring_atoms:
                        continue
                    p1 = positions[bond.a1_id]
                    p2 = positions[bond.a2_id]
                    dx = p2.x() - p1.x()
                    dy = p2.y() - p1.y()
                    dist = math.hypot(dx, dy) or 1e-3
                    delta = (dist - target_len) / dist
                    p1 = QPointF(p1.x() + dx * 0.5 * delta, p1.y() + dy * 0.5 * delta)
                    p2 = QPointF(p2.x() - dx * 0.5 * delta, p2.y() - dy * 0.5 * delta)
                    positions[bond.a1_id] = p1
                    positions[bond.a2_id] = p2

                for aid in ring_atoms:
                    neighbors = [n for n in neighbor_map.get(aid, []) if n in ring_atoms]
                    if len(neighbors) != 2:
                        continue
                    n1, n2 = neighbors
                    center = positions[aid]
                    p1 = positions[n1]
                    p2 = positions[n2]
                    a1 = angle_deg(center, p1)
                    a2 = angle_deg(center, p2)
                    current = angle_distance_deg(a1, a2)

                    has_triple = False
                    has_double = False
                    for bond in bonds_in_selection:
                        if bond.a1_id != aid and bond.a2_id != aid:
                            continue
                        if bond.order >= 3:
                            has_triple = True
                        if bond.order == 2 or bond.is_aromatic:
                            has_double = True
                    if has_triple:
                        target = 180.0
                    elif has_double:
                        target = 120.0
                    else:
                        target = SP3_BOND_ANGLE_DEG

                    delta = target - current
                    if abs(delta) < 1e-2:
                        continue
                    delta = max(-angle_step_deg, min(angle_step_deg, delta))
                    p1_new = rotate_point(p1, center, delta * 0.5)
                    p2_new = rotate_point(p2, center, -delta * 0.5)
                    positions[n1] = p1_new
                    positions[n2] = p2_new

            # Normalize trees attached to ring core
            normalize_tree(component, list(ring_atoms))

        after_center = self._center_for_atoms(list(atom_ids)) or QPointF(0.0, 0.0)
        # Use computed positions for center
        xs = [positions[aid].x() for aid in atom_ids]
        ys = [positions[aid].y() for aid in atom_ids]
        if xs and ys:
            after_center = QPointF(sum(xs) / len(xs), sum(ys) / len(ys))
        shift = QPointF(before_center.x() - after_center.x(), before_center.y() - after_center.y())

        after = {}
        for aid in atom_ids:
            p = positions[aid]
            after[aid] = (p.x() + shift.x(), p.y() + shift.y())

        from chemuson.gui.commands import MoveAtomsCommand
        cmd = MoveAtomsCommand(self.model, self, before, after)
        self.undo_stack.push(cmd)
        self._update_selection_overlay()

    def _find_ring_from(self, start: int, adjacency: dict, max_size: int) -> Optional[list]:
        """Find a ring starting from the given atom using DFS."""
        # Simple BFS to find shortest cycle
        from collections import deque
        
        queue = deque([(start, [start])])
        visited = {start}
        
        while queue:
            current, path = queue.popleft()
            if len(path) > max_size:
                continue
            
            for neighbor in adjacency.get(current, []):
                if neighbor == start and len(path) >= 5:
                    return path
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append((neighbor, path + [neighbor]))
        
        return None

    def _get_atom_at(self, x: float, y: float) -> Optional[int]:
        """Método auxiliar para  get atom at.

        Args:
            x: Descripción del parámetro.
            y: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        for atom in self.model.atoms.values():
            dx = atom.x - x
            dy = atom.y - y
            if (dx * dx + dy * dy) < (ATOM_HIT_RADIUS * ATOM_HIT_RADIUS):
                return atom.id
        return None

    def _get_bond_at(self, scene_pos: QPointF) -> Optional[int]:
        """Método auxiliar para  get bond at.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        for item in self.scene.items(scene_pos):
            if isinstance(item, BondItem):
                return item.bond_id
        return None

    def _structure_bbox(self) -> Optional[QRectF]:
        """Método auxiliar para  structure bbox.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        rect: Optional[QRectF] = None

        def extend(candidate: QRectF) -> None:
            """Método auxiliar para extend.

            Args:
                candidate: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
            nonlocal rect
            if not candidate.isValid() or candidate.isNull():
                return
            rect = candidate if rect is None else rect.united(candidate)

        def extend_atom_bounds(atom_id: int) -> None:
            """Método auxiliar para extend atom bounds.

            Args:
                atom_id: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
            item = self.atom_items.get(atom_id)
            if item is None or item.scene() is not self.scene:
                return
            if item.pen().style() != Qt.PenStyle.NoPen or item.brush().style() != Qt.BrushStyle.NoBrush:
                extend(item.sceneBoundingRect())
            if item.label.isVisible():
                extend(item.label.sceneBoundingRect())
            if item.charge_label.isVisible():
                extend(item.charge_label.sceneBoundingRect())
            overlays = self._implicit_h_overlays.get(atom_id)
            if overlays:
                for line_item, text_item in overlays:
                    if line_item.scene() is self.scene and line_item.isVisible():
                        extend(line_item.sceneBoundingRect())
                    if text_item.scene() is self.scene and text_item.isVisible():
                        extend(text_item.sceneBoundingRect())

        for atom_id in self.model.atoms.keys():
            extend_atom_bounds(atom_id)
        for bond_id in self.model.bonds.keys():
            item = self.bond_items.get(bond_id)
            if item is None:
                continue
            extend(item.sceneBoundingRect())
        return rect

    def _analysis_graph_and_bbox(self) -> tuple[Optional[MolGraph], Optional[QRectF]]:
        """Método auxiliar para  analysis graph and bbox.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atom_ids, bonds = self._selected_structure_ids()
        if atom_ids:
            return self._build_selection_graph(atom_ids, bonds), self._selected_items_bbox()
        if self.model.atoms:
            return self.model, self._structure_bbox()
        return None, None

    def _implicit_h_for_graph(self, graph: MolGraph, atom_id: int, element: str) -> int:
        """Método auxiliar para  implicit h for graph.

        Args:
            graph: Descripción del parámetro.
            atom_id: Descripción del parámetro.
            element: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if element not in IMPLICIT_H_ELEMENTS:
            return 0
        try:
            return int(max(0, graph.implicit_h_count(atom_id)))
        except Exception:
            return 0

    def _analysis_atom_counts(self, graph: MolGraph) -> dict[str, int]:
        """Método auxiliar para  analysis atom counts.

        Args:
            graph: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        counts: dict[str, int] = {}
        for atom in graph.atoms.values():
            counts[atom.element] = counts.get(atom.element, 0) + 1
            assigned_h = (
                graph.assigned_hydrogen_count(atom.id)
                if hasattr(graph, "assigned_hydrogen_count")
                else int(getattr(atom, "explicit_h", 0) or 0)
            )
            if assigned_h:
                counts["H"] = counts.get("H", 0) + int(assigned_h)
        for atom in graph.atoms.values():
            if atom.element == "H":
                continue
            implicit_h = self._implicit_h_for_graph(graph, atom.id, atom.element)
            if implicit_h > 0:
                counts["H"] = counts.get("H", 0) + implicit_h
        return {element: count for element, count in counts.items() if count > 0}

    @staticmethod
    def _analysis_formula(counts: dict[str, int]) -> str:
        """Método auxiliar para  analysis formula.

        Args:
            counts: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if not counts:
            return ""
        order: list[str] = []
        if "C" in counts:
            order.append("C")
        if "H" in counts:
            order.append("H")
        for element in sorted(e for e in counts.keys() if e not in {"C", "H"}):
            order.append(element)
        parts = []
        for element in order:
            count = counts.get(element, 0)
            if count <= 0:
                continue
            parts.append(element if count == 1 else f"{element}{count}")
        return "".join(parts)

    @staticmethod
    def _analysis_exact_mass(counts: dict[str, int]) -> Optional[float]:
        """Método auxiliar para  analysis exact mass.

        Args:
            counts: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        total = 0.0
        for element, count in counts.items():
            mass = MONOISOTOPIC_MASSES.get(element)
            if mass is None:
                return None
            total += mass * count
        return total

    @staticmethod
    def _analysis_molecular_weight(counts: dict[str, int]) -> Optional[float]:
        """Método auxiliar para  analysis molecular weight.

        Args:
            counts: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        total = 0.0
        for element, count in counts.items():
            mass = ATOMIC_WEIGHTS.get(element)
            if mass is None:
                return None
            total += mass * count
        return total

    @staticmethod
    def _analysis_convolve(
        distribution: list[tuple[float, float]],
        isotopes: list[tuple[float, float]],
    ) -> list[tuple[float, float]]:
        """Método auxiliar para  analysis convolve.

        Args:
            distribution: Descripción del parámetro.
            isotopes: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if not distribution:
            return []
        new_dist: dict[float, float] = {}
        for mass, prob in distribution:
            for iso_mass, iso_prob in isotopes:
                combined_mass = round(mass + iso_mass, 4)
                new_dist[combined_mass] = new_dist.get(combined_mass, 0.0) + prob * iso_prob
        if not new_dist:
            return []
        max_prob = max(new_dist.values())
        min_prob = max_prob * 1e-8
        pruned = {m: p for m, p in new_dist.items() if p >= min_prob}
        if len(pruned) > ANALYSIS_DIST_KEEP:
            items = sorted(pruned.items(), key=lambda item: item[1], reverse=True)[:ANALYSIS_DIST_KEEP]
            pruned = dict(items)
        return list(pruned.items())

    def _analysis_isotope_peaks(self, counts: dict[str, int]) -> Optional[list[tuple[float, float]]]:
        """Método auxiliar para  analysis isotope peaks.

        Args:
            counts: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        distribution: list[tuple[float, float]] = [(0.0, 1.0)]
        for element in sorted(counts.keys()):
            isotopes = ISOTOPE_ABUNDANCES.get(element)
            if isotopes is None:
                return None
            for _ in range(counts[element]):
                distribution = self._analysis_convolve(distribution, isotopes)
        if not distribution:
            return None
        binned: dict[float, float] = {}
        for mass, prob in distribution:
            mass_key = round(mass, 2)
            binned[mass_key] = binned.get(mass_key, 0.0) + prob
        max_prob = max(binned.values())
        if max_prob <= 0:
            return None
        peaks = [(mass, (prob / max_prob) * 100.0) for mass, prob in binned.items()]
        peaks = [peak for peak in peaks if peak[1] >= ANALYSIS_MIN_PEAK_PERCENT]
        peaks.sort(key=lambda item: item[1], reverse=True)
        if len(peaks) > ANALYSIS_MAX_PEAKS:
            peaks = peaks[:ANALYSIS_MAX_PEAKS]
        return peaks

    def _analysis_elemental_line(self, counts: dict[str, int], molecular_weight: Optional[float]) -> Optional[str]:
        """Método auxiliar para  analysis elemental line.

        Args:
            counts: Descripción del parámetro.
            molecular_weight: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if molecular_weight is None or molecular_weight <= 0:
            return None
        order: list[str] = []
        if "C" in counts:
            order.append("C")
        if "H" in counts:
            order.append("H")
        for element in sorted(e for e in counts.keys() if e not in {"C", "H"}):
            order.append(element)
        parts = []
        for element in order:
            weight = ATOMIC_WEIGHTS.get(element)
            if weight is None:
                continue
            percent = (weight * counts[element] / molecular_weight) * 100.0
            parts.append(f"{element}, {percent:.2f}")
        if not parts:
            return None
        return "Elemental Analysis: " + "; ".join(parts)

    def _analysis_build_text(self, graph: MolGraph, mode: str) -> str:
        """Método auxiliar para  analysis build text.

        Args:
            graph: Descripción del parámetro.
            mode: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        counts = self._analysis_atom_counts(graph)
        if not counts:
            return ""
        formula = self._analysis_formula(counts)
        exact_mass = self._analysis_exact_mass(counts)
        molecular_weight = self._analysis_molecular_weight(counts)
        peaks = self._analysis_isotope_peaks(counts)

        lines: list[str] = []
        if mode in {"name", "all", "iupac"}:
            iupac = self.current_iupac_name(graph)
            if mode in {"name", "all"}:
                lines.append(f"Nombre IUPAC: {iupac}")
            else:
                lines.append(iupac)
        if mode in {"name", "all"}:
            smiles = ""
            try:
                smiles = molgraph_to_smiles(graph)
            except Exception:
                smiles = ""
            lines.append(f"SMILES: {smiles or 'N/D'}")
        if mode in {"formula", "all"}:
            lines.append(f"Chemical Formula: {formula}")
        if mode in {"exact", "all"}:
            lines.append(
                f"Exact Mass: {exact_mass:.2f}" if exact_mass is not None else "Exact Mass: N/D"
            )
        if mode in {"weight", "all"}:
            lines.append(
                f"Molecular Weight: {molecular_weight:.2f}"
                if molecular_weight is not None
                else "Molecular Weight: N/D"
            )
        if mode in {"mz", "all"}:
            if peaks:
                formatted = ", ".join(f"{mass:.2f} ({percent:.1f}%)" for mass, percent in peaks)
                lines.append(f"m/z: {formatted}")
            else:
                lines.append("m/z: N/D")
        if mode in {"elemental", "all"}:
            elemental_line = self._analysis_elemental_line(counts, molecular_weight)
            lines.append(elemental_line or "Elemental Analysis: N/D")
        return "\n".join(lines)

    def current_name_options(self) -> NameOptions:
        """Construye opciones de nomenclatura según preferencias del documento."""
        advanced = bool(getattr(self, "name_advanced_enabled", True))
        isolated = bool(getattr(self, "name_rdkit_isolated", True))
        return NameOptions(
            enable_experimental=advanced,
            enable_special_templates=advanced,
            enable_advanced_stereo=advanced,
            allow_coordination=advanced,
            rdkit_isolated=isolated,
        )

    def current_iupac_name(self, graph: Optional[MolGraph] = None) -> str:
        """Devuelve nombre IUPAC actual con degradación segura a `N/D`."""
        target = graph or self.model
        try:
            return iupac_name(target, self.current_name_options())
        except Exception:
            return "N/D"

    def _insert_analysis_text(
        self,
        text: str,
        bbox: Optional[QRectF],
        scene_pos: QPointF,
    ) -> None:
        """Método auxiliar para  insert analysis text.

        Args:
            text: Descripción del parámetro.
            bbox: Descripción del parámetro.
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = TextAnnotationItem(text, 0.0, 0.0)
        self._apply_text_settings(item)
        item.setTextInteractionFlags(Qt.TextInteractionFlag.NoTextInteraction)
        try:
            item.document().setDocumentMargin(0)
            item.document().setDefaultStyleSheet("body { background: transparent; }")
        except Exception:
            pass
        if bbox is not None and bbox.isValid():
            x = bbox.left()
            y = bbox.bottom() + ANALYSIS_MARGIN_PX
        else:
            x = scene_pos.x() + 10.0
            y = scene_pos.y() + 10.0
        item.setPos(x, y)
        self.scene.clearSelection()
        self.undo_stack.push(AddTextItemCommand(self, item))
        try:
            item.setSelected(True)
        except RuntimeError:
            pass

    def _run_analysis_action(self, mode: str, scene_pos: QPointF) -> None:
        """Método auxiliar para  run analysis action.

        Args:
            mode: Descripción del parámetro.
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        graph, bbox = self._analysis_graph_and_bbox()
        if graph is None:
            return
        text = self._analysis_build_text(graph, mode)
        if not text:
            return
        self._insert_analysis_text(text, bbox, scene_pos)

    def run_analysis(self, mode: str) -> None:
        """Método auxiliar para run analysis.

        Args:
            mode: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        scene_pos = self._last_scene_pos
        if scene_pos is None:
            try:
                center = self.viewport().rect().center()
                scene_pos = self.mapToScene(center)
            except Exception:
                scene_pos = QPointF(0.0, 0.0)
        self._run_analysis_action(mode, scene_pos)

    @staticmethod
    def _arrow_length(item: ArrowItem) -> float:
        """Devuelve la longitud geométrica de una anotación lineal."""
        start = item.start_point()
        end = item.end_point()
        return math.hypot(end.x() - start.x(), end.y() - start.y())

    @staticmethod
    def _centered_arrow_segment(
        start: QPointF,
        end: QPointF,
        target_length: float,
    ) -> tuple[QPointF, QPointF]:
        """Reconstruye un segmento con la longitud dada preservando centro y ángulo."""
        dx = end.x() - start.x()
        dy = end.y() - start.y()
        length = math.hypot(dx, dy)
        if length <= 1e-6:
            ux, uy = 1.0, 0.0
        else:
            ux, uy = dx / length, dy / length
        center = QPointF((start.x() + end.x()) * 0.5, (start.y() + end.y()) * 0.5)
        half = max(0.5, float(target_length)) * 0.5
        return (
            QPointF(center.x() - ux * half, center.y() - uy * half),
            QPointF(center.x() + ux * half, center.y() + uy * half),
        )

    def _set_line_arrow_length(
        self,
        target_length: float,
        *,
        items: Optional[list[ArrowItem]] = None,
        text: str = "Set line length",
    ) -> bool:
        """Aplica una longitud uniforme a líneas rectas seleccionadas."""
        target = max(0.5, float(target_length))
        line_items = [
            item
            for item in (self._selected_line_arrow_items() if items is None else items)
            if self._is_uniform_line_arrow(item) and item.scene() is self.scene
        ]
        before: dict[ArrowItem, tuple[QPointF, QPointF]] = {}
        after: dict[ArrowItem, tuple[QPointF, QPointF]] = {}
        for item in line_items:
            start = item.start_point()
            end = item.end_point()
            new_start, new_end = self._centered_arrow_segment(start, end, target)
            if self._point_equal(start, new_start) and self._point_equal(end, new_end):
                continue
            before[item] = (start, end)
            after[item] = (new_start, new_end)
        if not after:
            return False
        self.undo_stack.push(MoveArrowItemsCommand(self, before, after, text=text))
        return True

    def _equalize_selected_line_arrow_lengths(
        self,
        reference_item: Optional[ArrowItem] = None,
    ) -> bool:
        """Iguala todas las líneas rectas seleccionadas a una longitud común."""
        items = self._selected_line_arrow_items()
        if len(items) < 2:
            return False
        reference = reference_item if self._is_uniform_line_arrow(reference_item) and reference_item in items else items[0]
        return self._set_line_arrow_length(
            self._arrow_length(reference),
            items=items,
            text="Equalize line lengths",
        )

    def _prompt_selected_line_arrow_length(
        self,
        reference_item: Optional[ArrowItem] = None,
    ) -> None:
        """Solicita una longitud y la aplica a las líneas rectas seleccionadas."""
        items = self._selected_line_arrow_items()
        if not items:
            return
        reference = reference_item if self._is_uniform_line_arrow(reference_item) and reference_item in items else items[0]
        current = self._arrow_length(reference)
        value, ok = QInputDialog.getDouble(
            self,
            "Longitud de líneas",
            "Longitud uniforme:",
            current,
            0.5,
            5000.0,
            1,
        )
        if not ok:
            return
        self._set_line_arrow_length(float(value), items=items, text="Set line length")

    def _sync_scene_with_model(self) -> None:
        """Método auxiliar para  sync scene with model.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._suspend_numbering_refresh = True
        try:
            for item in list(self.scene.items()):
                if isinstance(item, (AtomItem, BondItem)):
                    self.scene.removeItem(item)

            for atom_id in list(self.atom_items.keys()):
                self.remove_atom_item(atom_id)
            for bond_id in list(self.bond_items.keys()):
                self.remove_bond_item(bond_id)

            self._ring_centers.clear()
            ring_atoms: dict[int, set[int]] = {}
            for bond in self.model.bonds.values():
                if bond.ring_id is None:
                    continue
                ring_atoms.setdefault(bond.ring_id, set()).update({bond.a1_id, bond.a2_id})
            for ring_id, atom_ids in ring_atoms.items():
                xs = [self.model.get_atom(aid).x for aid in atom_ids]
                ys = [self.model.get_atom(aid).y for aid in atom_ids]
                center = (sum(xs) / len(xs), sum(ys) / len(ys))
                self.register_ring_center(ring_id, center)

            for atom in self.model.atoms.values():
                self.add_atom_item(atom)
            for bond in self.model.bonds.values():
                self.add_bond_item(bond)

            self.refresh_atom_visibility()
            self.refresh_aromatic_circles()
        finally:
            self._suspend_numbering_refresh = False
        self.recompute_numbering()

    def show_valence_errors(self, errors: Iterable[int]) -> None:
        """Método auxiliar para show valence errors.

        Args:
            errors: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        error_ids = set(errors)
        for atom_id, item in self.atom_items.items():
            item.set_valence_error(atom_id in error_ids)

    def validate_structure(self) -> list[int]:
        """Valida estructura y aplica resaltado de valencias inválidas."""
        errors = list(self.model.validate())
        self.show_valence_errors(errors)
        return errors

    def _kekulize_aromatic_bonds(self, seed_atoms: Optional[Iterable[int]] = None) -> None:
        """Método auxiliar para  kekulize aromatic bonds.

        Args:
            seed_atoms: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        aromatic_bonds = [bond for bond in self.model.bonds.values() if bond.is_aromatic]
        if not aromatic_bonds:
            return

        component_atoms: Optional[set[int]] = None
        if seed_atoms is not None:
            adjacency: dict[int, list[int]] = {}
            for bond in aromatic_bonds:
                adjacency.setdefault(bond.a1_id, []).append(bond.a2_id)
                adjacency.setdefault(bond.a2_id, []).append(bond.a1_id)
            seeds = [atom_id for atom_id in seed_atoms if atom_id in adjacency]
            if not seeds:
                return
            component_atoms = set()
            stack = list(seeds)
            while stack:
                node = stack.pop()
                if node in component_atoms:
                    continue
                component_atoms.add(node)
                for neighbor in adjacency.get(node, []):
                    if neighbor not in component_atoms:
                        stack.append(neighbor)

        display_orders = kekulize_display_orders(self.model, seed_atoms=component_atoms)

        if display_orders is None:
            try:
                import networkx as nx
            except Exception:
                nx = None

            if nx is not None:
                graph = nx.Graph()
                bond_by_pair: dict[frozenset[int], int] = {}
                for bond in aromatic_bonds:
                    if component_atoms is not None and (
                        bond.a1_id not in component_atoms or bond.a2_id not in component_atoms
                    ):
                        continue
                    pair = frozenset({bond.a1_id, bond.a2_id})
                    bond_by_pair[pair] = bond.id
                    weight = 10_000_000 - (min(pair) * 1000 + max(pair))
                    graph.add_edge(bond.a1_id, bond.a2_id, weight=weight)
                matching = nx.max_weight_matching(graph, maxcardinality=True, weight="weight")
                display_orders = {bond_id: 1 for bond_id in bond_by_pair.values()}
                for u, v in matching:
                    bond_id = bond_by_pair.get(frozenset({u, v}))
                    if bond_id is not None:
                        display_orders[bond_id] = 2

        if display_orders is None:
            display_orders = {}
            used_atoms: set[int] = set()
            sorted_bonds = sorted(
                aromatic_bonds,
                key=lambda bond: (min(bond.a1_id, bond.a2_id), max(bond.a1_id, bond.a2_id)),
            )
            for bond in sorted_bonds:
                if component_atoms is not None and (
                    bond.a1_id not in component_atoms or bond.a2_id not in component_atoms
                ):
                    continue
                display_orders[bond.id] = 1
            for bond in sorted_bonds:
                if component_atoms is not None and (
                    bond.a1_id not in component_atoms or bond.a2_id not in component_atoms
                ):
                    continue
                if bond.a1_id in used_atoms or bond.a2_id in used_atoms:
                    continue
                display_orders[bond.id] = 2
                used_atoms.add(bond.a1_id)
                used_atoms.add(bond.a2_id)

        external_double_atoms: set[int] = set()
        for bond in self.model.bonds.values():
            if bond.is_aromatic:
                continue
            if bond.order >= 2:
                external_double_atoms.add(bond.a1_id)
                external_double_atoms.add(bond.a2_id)

        for bond in aromatic_bonds:
            if component_atoms is not None and (
                bond.a1_id not in component_atoms or bond.a2_id not in component_atoms
            ):
                continue
            new_order = display_orders.get(bond.id, 1)
            if bond.a1_id in external_double_atoms or bond.a2_id in external_double_atoms:
                new_order = 1
            if bond.display_order != new_order or bond.order != new_order:
                bond.display_order = new_order
                bond.order = new_order
                self.update_bond_item(bond.id)
        for bond in aromatic_bonds:
            if component_atoms is not None and (
                bond.a1_id not in component_atoms or bond.a2_id not in component_atoms
            ):
                continue
            self.update_bond_item(bond.id)
        self._assign_ring_ids_for_aromatic_cycles()
        self._refresh_aromatic_ring_contexts()

    def _assign_ring_ids_for_aromatic_cycles(self) -> None:
        """Método auxiliar para  assign ring ids for aromatic cycles.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        adjacency: dict[int, list[tuple[int, int]]] = {}
        for bond in self.model.bonds.values():
            if not bond.is_aromatic:
                continue
            adjacency.setdefault(bond.a1_id, []).append((bond.a2_id, bond.id))
            adjacency.setdefault(bond.a2_id, []).append((bond.a1_id, bond.id))

        for bond in self.model.bonds.values():
            if not bond.is_aromatic or bond.ring_id is not None:
                continue
            cycle = self._find_cycle_for_bond(bond, adjacency)
            if not cycle:
                continue
            ring_id = None
            for bond_id in cycle["bond_ids"]:
                existing = self.model.get_bond(bond_id)
                if existing.ring_id is not None:
                    ring_id = existing.ring_id
                    break
            if ring_id is None:
                ring_id = self.allocate_ring_id()
            for bond_id in cycle["bond_ids"]:
                b = self.model.get_bond(bond_id)
                if b.ring_id is None:
                    b.ring_id = ring_id
                    self.update_bond_item(bond_id)
            center = self._center_for_atoms(cycle["atom_ids"])
            if center is not None:
                self.register_ring_center(ring_id, (center.x(), center.y()))

    def _find_cycle_for_bond(
        self, bond: Bond, adjacency: dict[int, list[tuple[int, int]]]
    ) -> Optional[dict[str, list[int]]]:
        """Método auxiliar para  find cycle for bond.

        Args:
            bond: Descripción del parámetro.
            adjacency: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        a1 = bond.a1_id
        a2 = bond.a2_id
        if a1 not in adjacency or a2 not in adjacency:
            return None
        from collections import deque

        queue = deque([(a1, [a1], [])])
        while queue:
            node, path_atoms, path_bonds = queue.popleft()
            if len(path_atoms) > 8:
                continue
            for neighbor, bond_id in adjacency.get(node, []):
                if bond_id == bond.id:
                    continue
                if neighbor == a2:
                    if len(path_atoms) >= 4:
                        atom_ids = path_atoms + [a2]
                        bond_ids = path_bonds + [bond_id, bond.id]
                        return {"atom_ids": atom_ids, "bond_ids": bond_ids}
                    continue
                if neighbor in path_atoms:
                    continue
                queue.append((neighbor, path_atoms + [neighbor], path_bonds + [bond_id]))
        return None

    def _build_aromatic_edges(
        self,
        vertex_defs: List[Tuple[Optional[int], float, float]],
        ring_size: int,
    ) -> List[Tuple[int, int, int, BondStyle, BondStereo, bool]]:
        """Método auxiliar para  build aromatic edges.

        Args:
            vertex_defs: Descripción del parámetro.
            ring_size: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if ring_size % 2 != 0:
            return [
                (i, (i + 1) % ring_size, 1, BondStyle.PLAIN, BondStereo.NONE, True)
                for i in range(ring_size)
            ]

        constraints: dict[int, int] = {}
        for i in range(ring_size):
            j = (i + 1) % ring_size
            a_id = vertex_defs[i][0]
            b_id = vertex_defs[j][0]
            if a_id is None or b_id is None:
                continue
            existing = self.model.find_bond_between(a_id, b_id)
            if existing is not None:
                constraints[i] = 2 if existing.order >= 2 else 1

        def atom_has_double(atom_id: Optional[int], other_id: Optional[int]) -> bool:
            """Método auxiliar para atom has double.

            Args:
                atom_id: Descripción del parámetro.
                other_id: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
            if atom_id is None:
                return False
            for bond in self.model.bonds.values():
                if bond.order < 2:
                    continue
                if bond.a1_id == atom_id and bond.a2_id != other_id:
                    return True
                if bond.a2_id == atom_id and bond.a1_id != other_id:
                    return True
            return False

        def score_phase(phase: int) -> int:
            """Método auxiliar para score phase.

            Args:
                phase: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
            score = 0
            for idx, desired in constraints.items():
                current = 2 if ((idx + phase) % 2 == 0) else 1
                if current != desired:
                    score += 1
            return score

        phase = 0 if score_phase(0) <= score_phase(1) else 1
        edges: List[Tuple[int, int, int, BondStyle, BondStereo, bool]] = []
        for i in range(ring_size):
            j = (i + 1) % ring_size
            order = 2 if ((i + phase) % 2 == 0) else 1
            a_id = vertex_defs[i][0]
            b_id = vertex_defs[j][0]
            if order == 2 and (atom_has_double(a_id, b_id) or atom_has_double(b_id, a_id)):
                order = 1
            edges.append((i, j, order, BondStyle.PLAIN, BondStereo.NONE, True))
        return edges

    def _is_on_paper(self, x: float, y: float) -> bool:
        """Método auxiliar para  is on paper.

        Args:
            x: Descripción del parámetro.
            y: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        return (
            PAPER_MARGIN <= x <= self.paper_width - PAPER_MARGIN
            and PAPER_MARGIN <= y <= self.paper_height - PAPER_MARGIN
        )
