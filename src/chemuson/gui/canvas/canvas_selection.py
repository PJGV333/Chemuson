from __future__ import annotations

from ._shared import *

class CanvasSelectionMixin:
    def _delete_selection(
        self,
        atom_ids: set[int],
        bond_ids: set[int],
        arrow_items: Iterable = (),
        bracket_items: Iterable = (),
        text_items: Iterable = (),
        energy_diagram_items: Iterable = (),
        semantic_diagram_items: Iterable = (),
        orbital_items: Iterable = (),
        wavy_items: Iterable = (),
        image_items: Iterable = (),
        plate_items: Iterable = (),
    ) -> None:
        """Método auxiliar para  delete selection.

        Args:
            atom_ids: Descripción del parámetro.
            bond_ids: Descripción del parámetro.
            arrow_items: Descripción del parámetro.
            bracket_items: Descripción del parámetro.
            text_items: Descripción del parámetro.
            wavy_items: Descripción del parámetro.
            image_items: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        extra_text_items = list(text_items)
        extra_wavy_items = list(wavy_items)
        if atom_ids and self._electron_dots:
            for item in list(self._electron_dots):
                anchor_id = item.data(ELECTRON_ANCHOR_ROLE)
                if anchor_id in atom_ids:
                    extra_text_items.append(item)
        if atom_ids and self._wavy_anchors:
            for item in list(self._wavy_anchors):
                anchor_id = item.data(WAVY_ANCHOR_ROLE)
                if anchor_id in atom_ids:
                    extra_wavy_items.append(item)
        # Prioritize delete logic
        cmd = DeleteSelectionCommand(
            self.model,
            self,
            atom_ids,
            bond_ids,
            arrow_items=arrow_items,
            bracket_items=bracket_items,
            text_items=extra_text_items,
            energy_diagram_items=energy_diagram_items,
            semantic_diagram_items=semantic_diagram_items,
            orbital_items=orbital_items,
            wavy_items=extra_wavy_items,
            image_items=image_items,
            plate_items=plate_items,
        )
        self.undo_stack.push(cmd)
        self.scene.clearSelection()

    def _selected_plate_items(self) -> list:
        return [
            item for item in self.scene.selectedItems() if isinstance(item, (TLCPlateItem, GelElectrophoresisItem))
        ]

    def _selection_snapshot(self) -> dict:
        """Captura la selección activa incluyendo elementos no moleculares."""
        return {
            "atom_ids": set(self.state.selected_atoms),
            "bond_ids": set(self.state.selected_bonds),
            "text_items": list(self._selected_text_items()),
            "arrow_items": list(self._selected_arrow_items()),
            "bracket_items": list(self._selected_bracket_items()),
            "energy_diagram_items": list(self._selected_energy_diagram_items()),
            "semantic_diagram_items": list(self._selected_semantic_diagram_items()),
            "orbital_items": list(self._selected_orbital_items()),
            "image_items": list(self._selected_image_items()),
            "plate_items": list(self._selected_plate_items()),
            "wavy_items": [
                item for item in self.scene.selectedItems() if isinstance(item, WavyAnchorItem)
            ],
        }

    def _restore_selection_snapshot(self, snapshot: dict) -> None:
        """Restaura una selección previamente capturada."""
        try:
            self.scene.blockSignals(True)
            self.scene.clearSelection()
            for atom_id in snapshot.get("atom_ids", ()):
                item = self.atom_items.get(atom_id)
                if item is not None and item.scene() is self.scene:
                    item.setSelected(True)
            for bond_id in snapshot.get("bond_ids", ()):
                item = self.bond_items.get(bond_id)
                if item is not None and item.scene() is self.scene:
                    item.setSelected(True)
            for key in (
                "text_items",
                "arrow_items",
                "bracket_items",
                "energy_diagram_items",
                "semantic_diagram_items",
                "orbital_items",
                "image_items",
                "plate_items",
                "wavy_items",
            ):
                for item in snapshot.get(key, ()):
                    try:
                        if item is not None and item.scene() is self.scene:
                            item.setSelected(True)
                    except RuntimeError:
                        continue
        finally:
            try:
                self.scene.blockSignals(False)
            except RuntimeError:
                pass
        self._sync_selection_from_scene()

    def delete_selection(self) -> None:
        """Método auxiliar para delete selection.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        selected_arrows = [
            item for item in self.scene.selectedItems() if isinstance(item, ArrowItem)
        ]
        selected_brackets = [
            item for item in self.scene.selectedItems() if isinstance(item, BracketItem)
        ]
        selected_text_items = [
            item for item in self.scene.selectedItems() if isinstance(item, TextAnnotationItem)
        ]
        selected_energy_diagrams = [
            item for item in self.scene.selectedItems() if isinstance(item, EnergyDiagramItem)
        ]
        selected_semantic_diagrams = [
            item for item in self.scene.selectedItems() if isinstance(item, CompositeDiagramItem)
        ]
        selected_orbitals = [
            item for item in self.scene.selectedItems() if isinstance(item, OrbitalAnnotationItem)
        ]
        selected_wavy = [
            item for item in self.scene.selectedItems() if isinstance(item, WavyAnchorItem)
        ]
        selected_images = [
            item for item in self.scene.selectedItems() if isinstance(item, ImageAnnotationItem)
        ]
        selected_plates = self._selected_plate_items()
        selected_spots_bands = [
            item for item in self.scene.selectedItems() if isinstance(item, (TLCSpotItem, GelBandItem))
        ]
        plates_from_spots = set()
        for sb in selected_spots_bands:
            lr = getattr(sb, "lane_ref", None)
            if lr is not None:
                parent_plate = lr.parentItem()
                if isinstance(parent_plate, (TLCPlateItem, GelElectrophoresisItem)):
                    plates_from_spots.add(parent_plate)
        all_plates = list(set(selected_plates) | plates_from_spots)
        for plate in all_plates:
            for lane in plate.lane_items:
                if hasattr(lane, "rf_labels"):
                    for spot, lbl in list(lane.rf_labels):
                        try:
                            self.scene.removeItem(spot)
                        except RuntimeError:
                            pass
                        try:
                            self.scene.removeItem(lbl)
                        except RuntimeError:
                            pass
                elif hasattr(lane, "bands"):
                    for band, lbl in list(lane.bands):
                        try:
                            self.scene.removeItem(band)
                        except RuntimeError:
                            pass
                        try:
                            self.scene.removeItem(lbl)
                        except RuntimeError:
                            pass
        if (
            not self.state.selected_atoms
            and not self.state.selected_bonds
            and not selected_arrows
            and not selected_brackets
            and not selected_text_items
            and not selected_energy_diagrams
            and not selected_semantic_diagrams
            and not selected_orbitals
            and not selected_wavy
            and not selected_images
            and not all_plates
        ):
            return
        self._delete_selection(
            set(self.state.selected_atoms),
            set(self.state.selected_bonds),
            arrow_items=selected_arrows,
            bracket_items=selected_brackets,
            text_items=selected_text_items,
            energy_diagram_items=selected_energy_diagrams,
            semantic_diagram_items=selected_semantic_diagrams,
            orbital_items=selected_orbitals,
            wavy_items=selected_wavy,
            image_items=selected_images,
            plate_items=all_plates,
        )

    def _delete_hovered(self) -> bool:
        """Método auxiliar para  delete hovered.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self.hovered_atom_id is not None:
            self._delete_selection({self.hovered_atom_id}, set())
            return True
        if self.hovered_bond_id is not None:
            self._delete_selection(set(), {self.hovered_bond_id})
            return True
        return False

    def _selected_structure_ids(self) -> tuple[set[int], list[Bond]]:
        """Método auxiliar para  selected structure ids.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atom_ids = self._selected_atom_ids_for_transform()
        if not atom_ids:
            return set(), []
        bonds: list[Bond] = []
        for bond in self.model.bonds.values():
            if bond.a1_id in atom_ids and bond.a2_id in atom_ids:
                bonds.append(bond)
        return atom_ids, self._unique_bonds_for_copy(bonds)

    @staticmethod
    def _bond_copy_priority(bond: Bond) -> int:
        """Prioriza la representación más informativa si hay enlaces duplicados."""
        style_bonus = 5 if bond.style == BondStyle.COORDINATION else 0
        aromatic_bonus = 50 if bond.is_aromatic else 0
        display_bonus = int(bond.display_order or 0)
        return int(bond.order or 1) * 10 + style_bonus + aromatic_bonus + display_bonus

    def _unique_bonds_for_copy(self, bonds: Iterable[Bond]) -> list[Bond]:
        """Elimina pares duplicados al copiar/exportar selección."""
        unique: dict[tuple[int, int], Bond] = {}
        for bond in sorted(bonds, key=lambda item: item.id):
            pair = (min(int(bond.a1_id), int(bond.a2_id)), max(int(bond.a1_id), int(bond.a2_id)))
            existing = unique.get(pair)
            if existing is None or self._bond_copy_priority(bond) > self._bond_copy_priority(existing):
                unique[pair] = bond
        return list(unique.values())

    @staticmethod
    def _is_large_clipboard_structure(graph: Optional[MolGraph]) -> bool:
        """Determina si una selección es suficientemente grande para exportación ligera."""
        if graph is None:
            return False
        return (
            len(graph.atoms) >= CLIPBOARD_LARGE_SELECTION_ATOM_THRESHOLD
            or len(graph.bonds) >= CLIPBOARD_LARGE_SELECTION_BOND_THRESHOLD
        )

    def _build_selection_graph(self, atom_ids: set[int], bonds: list[Bond]) -> MolGraph:
        """Método auxiliar para  build selection graph.

        Args:
            atom_ids: Descripción del parámetro.
            bonds: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        graph = MolGraph()
        for atom_id in atom_ids:
            atom = self.model.get_atom(atom_id)
            graph.add_atom(
                atom.element,
                atom.x,
                atom.y,
                atom_id=atom.id,
                charge=atom.charge,
                isotope=atom.isotope,
                radical_electrons=int(getattr(atom, "radical_electrons", 0) or 0),
                oxidation_state=getattr(atom, "oxidation_state", None),
                explicit_h=atom.explicit_h,
                group_h_cap=getattr(atom, "group_h_cap", None),
                mapping=atom.mapping,
                is_query=atom.is_query,
                is_explicit=atom.is_explicit,
                no_implicit=bool(getattr(atom, "no_implicit", False)),
                label_scale=getattr(atom, "label_scale", None),
                is_coordination_center=getattr(atom, "is_coordination_center", False),
                sphere_radius=getattr(atom, "sphere_radius", None),
                sphere_color=getattr(atom, "sphere_color", None),
                sphere_filled=bool(getattr(atom, "sphere_filled", True)),
                sphere_transparent=bool(getattr(atom, "sphere_transparent", False)),
                opacity=self.effective_atom_opacity(atom),
            )
        for bond in self._unique_bonds_for_copy(bonds):
            graph.add_bond(
                bond.a1_id,
                bond.a2_id,
                bond.order,
                bond_id=bond.id,
                style=bond.style,
                stereo=bond.stereo,
                is_aromatic=bond.is_aromatic,
                display_order=bond.display_order,
                is_query=bond.is_query,
                ring_id=bond.ring_id,
                length_px=bond.length_px,
                stroke_px=bond.stroke_px,
                color=bond.color,
                donor_atom_id=getattr(bond, "donor_atom_id", None),
                opacity=self.effective_bond_opacity(bond),
            )
        return graph

    def _build_selection_payload(self) -> Optional[dict]:
        """Método auxiliar para  build selection payload.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atom_ids, bonds = self._selected_structure_ids()
        arrows = self._selected_arrow_items()
        brackets = self._selected_bracket_items()
        texts = self._selected_text_items()
        energy_diagrams = self._selected_energy_diagram_items()
        semantic_diagrams = self._selected_semantic_diagram_items()
        orbitals = self._selected_orbital_items()
        images = self._selected_image_items()
        wavy_items = [
            item for item in self.scene.selectedItems() if isinstance(item, WavyAnchorItem)
        ]
        plates = self._selected_plate_items()

        if (
            not atom_ids
            and not bonds
            and not arrows
            and not brackets
            and not texts
            and not energy_diagrams
            and not semantic_diagrams
            and not orbitals
            and not images
            and not wavy_items
            and not plates
        ):
            return None

        bbox = self._selected_items_bbox()
        if bbox is None:
            return None
        left = bbox.left()
        top = bbox.top()

        atoms_payload = []
        for atom_id in atom_ids:
            atom = self.model.get_atom(atom_id)
            atoms_payload.append(
                {
                    "id": atom_id,
                    "element": atom.element,
                    "x": atom.x - left,
                    "y": atom.y - top,
                    "charge": atom.charge,
                    "formal_charge": int(getattr(atom, "formal_charge", atom.charge)),
                    "isotope": atom.isotope,
                    "radical_electrons": int(getattr(atom, "radical_electrons", 0) or 0),
                    "oxidation_state": getattr(atom, "oxidation_state", None),
                    "explicit_h": atom.explicit_h,
                    "group_h_cap": getattr(atom, "group_h_cap", None),
                    "mapping": atom.mapping,
                    "is_query": atom.is_query,
                    "is_explicit": atom.is_explicit,
                    "no_implicit": bool(getattr(atom, "no_implicit", False)),
                    "label_scale": getattr(atom, "label_scale", None),
                    "is_coordination_center": getattr(atom, "is_coordination_center", False),
                    "sphere_radius": getattr(atom, "sphere_radius", None),
                    "sphere_color": getattr(atom, "sphere_color", None),
                    "sphere_filled": bool(getattr(atom, "sphere_filled", True)),
                    "sphere_transparent": bool(getattr(atom, "sphere_transparent", False)),
                    "opacity": self.effective_atom_opacity(atom),
                    "anchor": self._group_anchor_overrides.get(atom_id),
                }
            )

        bonds_payload = []
        for bond in bonds:
                bonds_payload.append(
                    {
                        "a1": bond.a1_id,
                        "a2": bond.a2_id,
                        "order": bond.order,
                        "style": bond.style.value if bond.style is not None else BondStyle.PLAIN.value,
                        "type": bond.style.value if bond.style is not None else BondStyle.PLAIN.value,
                        "stereo": bond.stereo.value if bond.stereo is not None else BondStereo.NONE.value,
                        "is_aromatic": bond.is_aromatic,
                        "display_order": bond.display_order,
                        "is_query": bond.is_query,
                        "ring_id": bond.ring_id,
                        "length_px": bond.length_px,
                        "stroke_px": bond.stroke_px,
                        "color": bond.color,
                        "donor_atom_id": getattr(bond, "donor_atom_id", None),
                        "flex_curve_1": getattr(bond, "flex_curve_1", None),
                        "flex_curve_2": getattr(bond, "flex_curve_2", None),
                        "opacity": self.effective_bond_opacity(bond),
                    }
                )

        arrows_payload = []
        for item in arrows:
            start = item.start_point()
            end = item.end_point()
            arrows_payload.append(
                {
                    "start": [start.x() - left, start.y() - top],
                    "end": [end.x() - left, end.y() - top],
                    "kind": item.kind(),
                    "curve_factor": item.curve_factor(),
                    "stroke_px": item.stroke_px(),
                    "opacity": self.effective_item_opacity(item),
                }
            )

        brackets_payload = []
        for item in brackets:
            rect = item.base_rect()
            brackets_payload.append(
                {
                    "rect": [rect.x() - left, rect.y() - top, rect.width(), rect.height()],
                    "kind": getattr(item, "_kind", "[]"),
                    "padding": getattr(item, "_padding", None),
                    "stroke_px": item.stroke_px(),
                    "opacity": self.effective_item_opacity(item),
                }
            )

        texts_payload = []
        for item in texts:
            texts_payload.append(
                {
                    "text": item.toPlainText(),
                    "html": item.toHtml(),
                    "x": item.pos().x() - left,
                    "y": item.pos().y() - top,
                    "rotation": item.rotation(),
                    "font": item.font().toString(),
                    "color": item.defaultTextColor().name(),
                    "text_width": item.textWidth(),
                    "opacity": self.effective_item_opacity(item),
                }
            )

        energy_diagrams_payload = []
        for item in energy_diagrams:
            rect = item.display_rect()
            energy_diagrams_payload.append(
                {
                    "x": rect.x() - left,
                    "y": rect.y() - top,
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
                    "opacity": self.effective_item_opacity(item),
                }
            )

        semantic_diagrams_payload = []
        for item in semantic_diagrams:
            payload = item.to_json()
            payload["x"] = float(payload.get("x", 0.0)) - left
            payload["y"] = float(payload.get("y", 0.0)) - top
            payload["opacity"] = self.effective_item_opacity(item)
            semantic_diagrams_payload.append(payload)

        orbitals_payload = []
        for item in orbitals:
            anchor0 = item.anchor0()
            anchor1 = item.anchor1()
            orbitals_payload.append(
                {
                    "kind": item.kind(),
                    "anchor0": [anchor0.x() - left, anchor0.y() - top],
                    "anchor1": [anchor1.x() - left, anchor1.y() - top],
                    "stroke_shaded_lobes": item.stroke_shaded_lobes(),
                    "part_styles": item.part_styles(),
                    "z": item.zValue(),
                    "opacity": self.effective_item_opacity(item),
                }
            )

        images_payload = []
        for item in images:
            rect = item.display_rect()
            images_payload.append(
                {
                    "x": rect.x() - left,
                    "y": rect.y() - top,
                    "width": rect.width(),
                    "height": rect.height(),
                    "rotation": item.rotation(),
                    "z": item.zValue(),
                    "mime_type": item.mime_type(),
                    "data_b64": base64.b64encode(item.data_bytes()).decode("ascii"),
                    "source_name": item.source_name(),
                    "opacity": self.effective_item_opacity(item),
                }
            )

        wavy_payload = []
        for item in wavy_items:
            start = item.start_point()
            end = item.end_point()
            wavy_payload.append(
                {
                    "start": [start.x() - left, start.y() - top],
                    "end": [end.x() - left, end.y() - top],
                    "anchor_id": item.data(WAVY_ANCHOR_ROLE),
                    "angle": item.data(WAVY_ANCHOR_ANGLE_ROLE),
                    "length": item.data(WAVY_ANCHOR_LENGTH_ROLE),
                    "bond_id": item.data(WAVY_ANCHOR_BOND_ROLE),
                    "opacity": self.effective_item_opacity(item),
                }
            )

        plates_payload = []
        for plate in plates:
            plate_data = plate.to_dict()
            plate_data["pos"] = (plate.pos().x() - left, plate.pos().y() - top)
            plates_payload.append(plate_data)

        return {
            "atoms": atoms_payload,
            "bonds": bonds_payload,
            "arrows": arrows_payload,
            "brackets": brackets_payload,
            "texts": texts_payload,
            "energy_diagrams": energy_diagrams_payload,
            "semantic_diagrams": semantic_diagrams_payload,
            "orbitals": orbitals_payload,
            "images": images_payload,
            "wavy_items": wavy_payload,
            "plates": plates_payload,
        }

    def _paste_selection_payload(self, payload: dict) -> None:
        """Método auxiliar para  paste selection payload.

        Args:
            payload: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atoms = payload.get("atoms", [])
        bonds = payload.get("bonds", [])
        arrows = payload.get("arrows", [])
        brackets = payload.get("brackets", [])
        texts = payload.get("texts", [])
        energy_diagrams = payload.get("energy_diagrams", [])
        semantic_diagrams = payload.get("semantic_diagrams", [])
        orbitals = payload.get("orbitals", [])
        images = payload.get("images", [])
        wavy_items = payload.get("wavy_items", [])
        plates = payload.get("plates", [])
        if (
            not atoms
            and not bonds
            and not arrows
            and not brackets
            and not texts
            and not energy_diagrams
            and not semantic_diagrams
            and not orbitals
            and not images
            and not wavy_items
            and not plates
        ):
            return

        target = self._last_scene_pos
        if target is None:
            target = self.mapToScene(self.viewport().rect().center())
        dx = target.x()
        dy = target.y()

        ring_map: dict[int, int] = {}

        def map_ring(ring_id: Optional[int]) -> Optional[int]:
            """Método auxiliar para map ring.

            Args:
                ring_id: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
            if ring_id is None:
                return None
            if ring_id not in ring_map:
                ring_map[ring_id] = self.allocate_ring_id()
            return ring_map[ring_id]

        has_undo_items = bool(
            atoms
            or bonds
            or arrows
            or brackets
            or texts
            or energy_diagrams
            or semantic_diagrams
            or orbitals
            or images
            or wavy_items
            or plates
        )
        self.begin_validation_batch()
        if has_undo_items:
            self.undo_stack.beginMacro("Paste selection")
        id_map: Dict[int, int] = {}
        inserted_atom_ids: list[int] = []
        inserted_items: list[QGraphicsItem] = []
        inserted_pairs: set[frozenset[int]] = set()
        try:
            for atom_d in atoms:
                cmd = AddAtomCommand(
                    self.model,
                    self,
                    atom_d.get("element", "C"),
                    float(atom_d.get("x", 0.0)) + dx,
                    float(atom_d.get("y", 0.0)) + dy,
                    is_explicit=atom_d.get("is_explicit"),
                    charge=atom_d.get("formal_charge", atom_d.get("charge")),
                    isotope=atom_d.get("isotope"),
                    radical_electrons=int(atom_d.get("radical_electrons", 0) or 0),
                    oxidation_state=atom_d.get("oxidation_state"),
                    explicit_h=atom_d.get("explicit_h"),
                    group_h_cap=atom_d.get("group_h_cap"),
                    mapping=atom_d.get("mapping"),
                    is_query=bool(atom_d.get("is_query", False)),
                    anchor_override=atom_d.get("anchor"),
                    auto_hydrogens=False,
                    no_implicit=bool(atom_d.get("no_implicit", False)),
                    label_scale=atom_d.get("label_scale"),
                    is_coordination_center=bool(atom_d.get("is_coordination_center", False)),
                    sphere_radius=atom_d.get("sphere_radius"),
                    sphere_color=atom_d.get("sphere_color"),
                    sphere_filled=bool(atom_d.get("sphere_filled", True)),
                    sphere_transparent=bool(atom_d.get("sphere_transparent", False)),
                    opacity=atom_d.get("opacity"),
                )
                if has_undo_items:
                    self.undo_stack.push(cmd)
                else:
                    cmd.redo()
                if cmd.atom_id is not None:
                    id_map[int(atom_d.get("id"))] = cmd.atom_id
                    inserted_atom_ids.append(cmd.atom_id)

            for bond_d in bonds:
                a1 = id_map.get(int(bond_d.get("a1")))
                a2 = id_map.get(int(bond_d.get("a2")))
                if a1 is None or a2 is None:
                    continue
                pair = frozenset({a1, a2})
                if pair in inserted_pairs:
                    continue
                inserted_pairs.add(pair)
                if self.model.find_bond_between(a1, a2) is not None:
                    continue
                style = self._parse_bond_style_payload(bond_d)
                stereo = self._parse_bond_stereo_payload(bond_d)
                donor_atom_id = None
                raw_donor = bond_d.get("donor_atom_id")
                if raw_donor is not None:
                    try:
                        donor_atom_id = id_map.get(int(raw_donor))
                    except Exception:
                        donor_atom_id = None
                if style == BondStyle.COORDINATION:
                    donor_atom_id = self._infer_coordination_donor_atom(
                        a1,
                        a2,
                        preferred=donor_atom_id,
                    )
                cmd = AddBondCommand(
                    self.model,
                    self,
                    a1,
                    a2,
                    int(bond_d.get("order", 1)),
                    style,
                    stereo,
                    bool(bond_d.get("is_aromatic", False)),
                    display_order=bond_d.get("display_order"),
                    length_px=bond_d.get("length_px"),
                    stroke_px=bond_d.get("stroke_px"),
                    color=bond_d.get("color"),
                    ring_id=map_ring(bond_d.get("ring_id")),
                    donor_atom_id=donor_atom_id,
                    flex_curve_1=bond_d.get("flex_curve_1"),
                    flex_curve_2=bond_d.get("flex_curve_2"),
                    opacity=bond_d.get("opacity"),
                )
                if has_undo_items:
                    self.undo_stack.push(cmd)
                else:
                    cmd.redo()

            for arrow_d in arrows:
                start = arrow_d.get("start", [0.0, 0.0])
                end = arrow_d.get("end", [10.0, 0.0])
                kind = arrow_d.get("kind", "forward")
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
                start_pt = QPointF(float(start[0]) + dx, float(start[1]) + dy)
                end_pt = QPointF(float(end[0]) + dx, float(end[1]) + dy)
                cmd = AddArrowCommand(
                    self,
                    start_pt,
                    end_pt,
                    kind,
                    curve_factor=curve_factor_value,
                    stroke_px=stroke_px_value,
                    opacity=arrow_d.get("opacity"),
                )
                if has_undo_items:
                    self.undo_stack.push(cmd)
                else:
                    cmd.redo()

            for bracket_d in brackets:
                rect_vals = bracket_d.get("rect", [0.0, 0.0, 10.0, 10.0])
                rect = QRectF(
                    float(rect_vals[0]) + dx,
                    float(rect_vals[1]) + dy,
                    float(rect_vals[2]),
                    float(rect_vals[3]),
                )
                kind = bracket_d.get("kind", "[]")
                padding = bracket_d.get("padding")
                stroke_px = bracket_d.get("stroke_px")
                pair = self._split_bracket_kind(kind)
                if pair:
                    for side in pair:
                        cmd = AddBracketCommand(
                            self,
                            rect,
                            side,
                            padding=padding,
                            stroke_px=stroke_px,
                            opacity=bracket_d.get("opacity"),
                        )
                        if has_undo_items:
                            self.undo_stack.push(cmd)
                        else:
                            cmd.redo()
                        if cmd.item is not None:
                            inserted_items.append(cmd.item)
                else:
                    cmd = AddBracketCommand(
                        self,
                        rect,
                        kind,
                        padding=padding,
                        stroke_px=stroke_px,
                        opacity=bracket_d.get("opacity"),
                    )
                    if has_undo_items:
                        self.undo_stack.push(cmd)
                    else:
                        cmd.redo()
                    if cmd.item is not None:
                        inserted_items.append(cmd.item)

            for txt_d in texts:
                text_item = TextAnnotationItem(txt_d.get("text", ""), 0.0, 0.0)
                html = txt_d.get("html")
                if html:
                    text_item.setHtml(html)
                if "rotation" in txt_d:
                    text_item.setRotation(float(txt_d["rotation"]))
                if "font" in txt_d:
                    font = QFont()
                    font.fromString(txt_d["font"])
                    text_item.setFont(font)
                if "color" in txt_d:
                    text_item.setDefaultTextColor(QColor(txt_d["color"]))
                if "text_width" in txt_d:
                    text_item.setTextWidth(float(txt_d["text_width"]))
                self.set_graphics_item_opacity(text_item, txt_d.get("opacity"))
                text_item.setPos(
                    float(txt_d.get("x", 0.0)) + dx, float(txt_d.get("y", 0.0)) + dy
                )
                if has_undo_items:
                    self.undo_stack.push(AddTextItemCommand(self, text_item))
                else:
                    self.scene.addItem(text_item)
                inserted_items.append(text_item)

            for energy_d in energy_diagrams:
                item = EnergyDiagramItem(
                    energy_d.get("kind", DEFAULT_ENERGY_DIAGRAM_KIND),
                    label=energy_d.get("label"),
                    label_side=energy_d.get("label_side"),
                    occupancies=normalize_energy_occupancies(
                        energy_d.get("occupancies", ()),
                        kind=str(energy_d.get("kind", DEFAULT_ENERGY_DIAGRAM_KIND)),
                        box_count=energy_d.get("slot_count"),
                    ),
                    slot_count=energy_d.get("slot_count"),
                    style_payload=energy_d.get("style_payload"),
                    metadata=energy_d.get("metadata"),
                    width=float(energy_d.get("width", 1.0)),
                    height=float(energy_d.get("height", 1.0)),
                )
                item.set_display_rect(
                    QRectF(
                        float(energy_d.get("x", 0.0)) + dx,
                        float(energy_d.get("y", 0.0)) + dy,
                        float(energy_d.get("width", 1.0)),
                        float(energy_d.get("height", 1.0)),
                    )
                )
                item.setRotation(float(energy_d.get("rotation", 0.0)))
                item.setZValue(float(energy_d.get("z", 44.0)))
                self.set_graphics_item_opacity(item, energy_d.get("opacity"))
                if has_undo_items:
                    self.undo_stack.push(AddEnergyDiagramItemCommand(self, item))
                else:
                    self.readd_energy_diagram_item(item)
                inserted_items.append(item)

            for semantic_d in semantic_diagrams:
                try:
                    item = CompositeDiagramItem.from_json(
                        {
                            **dict(semantic_d),
                            "x": float(semantic_d.get("x", 0.0)) + dx,
                            "y": float(semantic_d.get("y", 0.0)) + dy,
                        }
                    )
                except Exception:
                    continue
                self.set_graphics_item_opacity(item, semantic_d.get("opacity"))
                if has_undo_items:
                    self.undo_stack.push(AddCompositeDiagramItemCommand(self, item))
                else:
                    self.readd_semantic_diagram_item(item)
                inserted_items.append(item)

            for orbital_d in orbitals:
                anchor0_vals = orbital_d.get("anchor0", [0.0, 0.0])
                anchor1_vals = orbital_d.get("anchor1", [0.0, -20.0])
                item = OrbitalAnnotationItem(
                    orbital_d.get("kind", DEFAULT_ORBITAL_KIND),
                    QPointF(float(anchor0_vals[0]) + dx, float(anchor0_vals[1]) + dy),
                    QPointF(float(anchor1_vals[0]) + dx, float(anchor1_vals[1]) + dy),
                    stroke_shaded_lobes=orbital_d.get("stroke_shaded_lobes"),
                    part_styles=orbital_d.get("part_styles"),
                )
                item.setZValue(float(orbital_d.get("z", 44.0)))
                self.set_graphics_item_opacity(item, orbital_d.get("opacity"))
                if has_undo_items:
                    self.undo_stack.push(AddOrbitalItemCommand(self, item))
                else:
                    self.readd_orbital_item(item)
                inserted_items.append(item)

            for idx, img_d in enumerate(images):
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
                        float(img_d.get("x", 0.0)) + dx,
                        float(img_d.get("y", 0.0)) + dy,
                        float(img_d.get("width", 1.0)),
                        float(img_d.get("height", 1.0)),
                    )
                )
                item.setRotation(float(img_d.get("rotation", 0.0)))
                item.setZValue(float(img_d.get("z", 8.0)))
                self.set_graphics_item_opacity(item, img_d.get("opacity"))
                if has_undo_items:
                    self.undo_stack.push(AddImageItemCommand(self, item))
                else:
                    self.readd_image_item(item)
                inserted_items.append(item)

            for anchor_d in wavy_items:
                start_vals = anchor_d.get("start", [0.0, 0.0])
                end_vals = anchor_d.get("end", [10.0, 0.0])
                item = WavyAnchorItem(
                    QPointF(float(start_vals[0]) + dx, float(start_vals[1]) + dy),
                    QPointF(float(end_vals[0]) + dx, float(end_vals[1]) + dy),
                    style=self.drawing_style,
                )
                anchor_id = anchor_d.get("anchor_id")
                angle = anchor_d.get("angle")
                length = anchor_d.get("length")
                mapped_anchor_id = None
                if anchor_id is not None:
                    try:
                        mapped_anchor_id = id_map.get(int(anchor_id))
                    except Exception:
                        mapped_anchor_id = None
                    if mapped_anchor_id is None and anchor_id in self.model.atoms:
                        mapped_anchor_id = int(anchor_id)
                if mapped_anchor_id is None:
                    continue
                item.setData(WAVY_ANCHOR_ROLE, mapped_anchor_id)
                if angle is not None:
                    item.setData(WAVY_ANCHOR_ANGLE_ROLE, float(angle))
                if length is not None:
                    item.setData(WAVY_ANCHOR_LENGTH_ROLE, float(length))
                self.set_graphics_item_opacity(item, anchor_d.get("opacity"))
                if has_undo_items:
                    self.undo_stack.push(AddWavyAnchorCommand(self, item))
                else:
                    self.readd_wavy_anchor_item(item)
                inserted_items.append(item)

            for plate_d in plates:
                plate = PlateItem.from_json(plate_d)
                if plate is None:
                    continue
                plate.load_dict(plate_d, scene=self.scene)
                orig_pos = plate_d.get("pos", (0, 0))
                paste_dx = dx if (dx != 0 or dy != 0) else 30.0
                paste_dy = dy if (dx != 0 or dy != 0) else 30.0
                plate.setPos(
                    float(orig_pos[0]) + paste_dx,
                    float(orig_pos[1]) + paste_dy,
                )
                self.scene.addItem(plate)
                if plate not in self.plate_items:
                    self.plate_items.append(plate)
                for lane in plate.lane_items:
                    if hasattr(lane, "rf_labels"):
                        for spot, _ in lane.rf_labels:
                            spot.setPos(spot.pos() + plate.pos())
                            spot.lane_ref = lane
                    elif hasattr(lane, "bands"):
                        for band, _ in lane.bands:
                            band.setPos(band.pos() + plate.pos())
                            band.lane_ref = lane
                if has_undo_items:
                    self.undo_stack.push(AddPlateItemCommand(self, plate))
                inserted_items.append(plate)
        finally:
            if has_undo_items:
                self.undo_stack.endMacro()
            self.end_validation_batch()
        if ring_map:
            self.refresh_ring_centers()
        self._select_inserted_items(inserted_atom_ids, inserted_items)

    def _select_inserted_items(
        self,
        atom_ids: Iterable[int] = (),
        items: Iterable[QGraphicsItem] = (),
    ) -> None:
        """Selecciona por defecto los elementos recién insertados."""
        target_items: list[QGraphicsItem] = []
        seen: set[int] = set()
        for atom_id in atom_ids:
            item = self.atom_items.get(atom_id)
            if item is None:
                continue
            marker = id(item)
            if marker in seen:
                continue
            seen.add(marker)
            target_items.append(item)
        for item in items:
            if item is None:
                continue
            marker = id(item)
            if marker in seen:
                continue
            seen.add(marker)
            target_items.append(item)
        if not target_items:
            return
        self.scene.clearSelection()
        for item in target_items:
            item.setSelected(True)
        self._sync_selection_from_scene()

    def has_copyable_selection(self) -> bool:
        """Indica si hay una selección que puede copiarse/cortarse."""
        return bool(
            self.state.selected_atoms
            or self.state.selected_bonds
            or self.scene.selectedItems()
        )

    def can_paste_from_clipboard(self) -> bool:
        """Indica si el portapapeles contiene un formato pegable por Chemuson."""
        mime = QApplication.clipboard().mimeData()
        if mime is None:
            return False
        return bool(
            mime.hasFormat("application/x-chemuson-selection")
            or mime.hasFormat("application/x-chemuson-text-items")
            or mime.hasFormat("chemical/x-mdl-molfile")
            or mime.hasUrls()
            or mime.hasText()
            or mime.hasFormat("image/png")
            or mime.hasImage()
            or mime.hasFormat("image/svg+xml")
        )

    def selected_semantic_diagram_item(self) -> CompositeDiagramItem | None:
        """Devuelve el diagrama semántico seleccionado si la selección es única."""
        selected_items = list(self.scene.selectedItems())
        if len(selected_items) != 1:
            return None
        item = selected_items[0]
        return item if isinstance(item, CompositeDiagramItem) else None

    def copy_to_clipboard(self) -> None:
        """Método auxiliar para copy to clipboard.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        selected_text_items = self._selected_text_items()
        selected_image_items = self._selected_image_items()
        selected_payload = self._build_selection_payload()
        atom_ids, bonds = self._selected_structure_ids()
        selected_arrows = bool(self._selected_arrow_items())
        selected_brackets = bool(self._selected_bracket_items())
        has_structure_selection = bool(atom_ids or bonds)
        only_text_selection = bool(selected_text_items) and not (
            has_structure_selection or selected_arrows or selected_brackets
        )
        if not self.model.atoms and not selected_text_items and not selected_image_items and selected_payload is None:
            return
        mime = QMimeData()
        if selected_payload is not None:
            mime.setData(
                "application/x-chemuson-selection",
                json.dumps(selected_payload).encode("utf-8"),
            )
        if has_structure_selection:
            graph = self._build_selection_graph(atom_ids, bonds)
        elif selected_payload is None:
            graph = self.model if self.model.atoms and not only_text_selection else None
        else:
            graph = None
        large_structure_selection = bool(
            has_structure_selection
            and selected_payload is not None
            and self._is_large_clipboard_structure(graph)
        )
        smiles = ""
        if graph is not None and graph.atoms and not large_structure_selection:
            try:
                molfile = molgraph_to_molfile(graph)
                smiles = molgraph_to_smiles(graph)
                mime.setData("chemical/x-mdl-molfile", molfile.encode("utf-8"))
            except Exception:
                pass
        if only_text_selection:
            data = {
                "text_items": [
                    {
                        "text": item.toPlainText(),
                        "html": item.toHtml(),
                        "x": item.pos().x(),
                        "y": item.pos().y(),
                        "rotation": item.rotation(),
                        "font": item.font().toString(),
                        "color": item.defaultTextColor().name(),
                        "opacity": self.effective_item_opacity(item),
                    }
                    for item in selected_text_items
                ]
            }
            mime.setData(
                "application/x-chemuson-text-items",
                json.dumps(data).encode("utf-8"),
            )

        has_selection = self.has_copyable_selection()
        image = self._render_scene_image(
            scale=(
                CLIPBOARD_LARGE_SELECTION_RENDER_SCALE
                if large_structure_selection
                else CLIPBOARD_RENDER_SCALE
            ),
            selected_only=has_selection,
            background=Qt.GlobalColor.white,
        )
        if image is not None:
            buffer = QBuffer()
            buffer.open(QBuffer.OpenModeFlag.WriteOnly)
            image.save(buffer, "PNG")
            mime.setData("image/png", buffer.data())
            mime.setImageData(image)
            try:
                png_b64 = base64.b64encode(bytes(buffer.data())).decode("ascii")
                html = f'<img src="data:image/png;base64,{png_b64}" alt="{smiles}">'
                mime.setHtml(html)
            except Exception:
                pass
        elif smiles:
            mime.setText(smiles)
        if not large_structure_selection:
            svg_data = self._render_scene_svg(selected_only=has_selection)
            if svg_data:
                mime.setData("image/svg+xml", svg_data)

        QApplication.clipboard().setMimeData(mime)

    def paste_from_clipboard(self) -> None:
        """Método auxiliar para paste from clipboard.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        clipboard = QApplication.clipboard()
        mime = clipboard.mimeData()
        if mime is None:
            return

        if mime.hasFormat("application/x-chemuson-selection"):
            try:
                payload = json.loads(bytes(mime.data("application/x-chemuson-selection")).decode("utf-8"))
                self._paste_selection_payload(payload)
                return
            except Exception:
                pass

        if mime.hasFormat("application/x-chemuson-text-items"):
            try:
                payload = json.loads(bytes(mime.data("application/x-chemuson-text-items")).decode("utf-8"))
                self._paste_text_items(payload)
                return
            except Exception:
                pass

        if mime.hasFormat("chemical/x-mdl-molfile"):
            molfile = bytes(mime.data("chemical/x-mdl-molfile")).decode("utf-8", errors="ignore")
            try:
                graph = molfile_to_molgraph(molfile)
                self._insert_molgraph(graph, select_inserted=True)
                return
            except Exception:
                pass

        if mime.hasUrls():
            if self._insert_images_from_clipboard(mime):
                return

        if mime.hasText():
            smiles = mime.text().strip()
            if smiles:
                try:
                    graph = smiles_to_molgraph(smiles)
                    self._insert_molgraph(graph, select_inserted=True)
                    return
                except Exception:
                    pass

        if mime.hasFormat("image/png") or mime.hasImage() or mime.hasFormat("image/svg+xml"):
            self._insert_images_from_clipboard(mime)

    def cut_to_clipboard(self) -> None:
        """Copy selected items to clipboard and delete them."""
        self.copy_to_clipboard()
        self.delete_selection()

    def duplicate_selection(self) -> None:
        """Duplicate selected items by copying to clipboard and pasting with offset."""
        if not self.has_copyable_selection():
            return
        self.copy_to_clipboard()
        clipboard = QApplication.clipboard()
        mime = clipboard.mimeData()
        if mime is None:
            return
        if mime.hasFormat("application/x-chemuson-selection"):
            try:
                payload = json.loads(bytes(mime.data("application/x-chemuson-selection")).decode("utf-8"))
                self._paste_selection_payload(payload)
            except Exception:
                pass

    def _paste_text_items(self, payload: dict) -> None:
        """Método auxiliar para  paste text items.

        Args:
            payload: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        items = payload.get("text_items", [])
        if not items:
            return
        created: list[TextAnnotationItem] = []
        for txt_d in items:
            text_item = TextAnnotationItem(txt_d.get("text", ""), 0.0, 0.0)
            html = txt_d.get("html")
            if html:
                text_item.setHtml(html)
            if "rotation" in txt_d:
                text_item.setRotation(float(txt_d["rotation"]))
            if "font" in txt_d:
                font = QFont()
                font.fromString(txt_d["font"])
                text_item.setFont(font)
            if "color" in txt_d:
                text_item.setDefaultTextColor(QColor(txt_d["color"]))
            if "text_width" in txt_d:
                text_item.setTextWidth(float(txt_d["text_width"]))
            self.set_graphics_item_opacity(text_item, txt_d.get("opacity"))
            text_item.setPos(float(txt_d.get("x", 0.0)), float(txt_d.get("y", 0.0)))
            self.scene.addItem(text_item)
            created.append(text_item)

        bbox: Optional[QRectF] = None
        for item in created:
            rect = item.sceneBoundingRect()
            bbox = rect if bbox is None else bbox.united(rect)
        if bbox is None:
            return
        target = self._last_scene_pos
        if target is None:
            target = self.mapToScene(self.viewport().rect().center())
        dx = target.x() - bbox.left()
        dy = target.y() - bbox.top()
        for item in created:
            item.setPos(item.pos() + QPointF(dx, dy))

        self.scene.clearSelection()
        for item in created:
            item.setSelected(True)
        self._sync_selection_from_scene()

    def _sync_selection_from_scene(self) -> None:
        """Método auxiliar para  sync selection from scene.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        try:
            selected_items = list(self.scene.selectedItems())
        except RuntimeError:
            return
        selected_atoms = set()
        selected_bonds = set()
        for item in selected_items:
            if isinstance(item, AtomItem):
                if item.atom_id in self.model.atoms:
                    selected_atoms.add(item.atom_id)
            elif isinstance(item, BondItem):
                if item.bond_id in self.model.bonds:
                    selected_bonds.add(item.bond_id)

        selected_text_items = [item for item in selected_items if isinstance(item, TextAnnotationItem)]
        if not selected_text_items and not selected_atoms and not selected_bonds:
            focus_item = self._active_text_edit_item()
            if focus_item is not None:
                selected_text_items = [focus_item]
        
        self.state.selected_atoms = selected_atoms
        self.state.selected_bonds = selected_bonds

        for atom_id, item in list(self.atom_items.items()):
            try:
                item.set_selected(atom_id in selected_atoms or atom_id == self.bond_anchor_id)
            except RuntimeError:
                self.atom_items.pop(atom_id, None)

        self._update_selection_overlay()

        # Emit selection signal
        details = {}
        if len(selected_atoms) == 1 and not selected_bonds and not selected_text_items:
            atom = self.model.get_atom(next(iter(selected_atoms)))
            details = {
                "type": "atom",
                "id": atom.id,
                "element": atom.element,
                "charge": int(getattr(atom, "formal_charge", atom.charge)),
                "x": atom.x,
                "y": atom.y,
            }
        elif len(selected_bonds) == 1 and not selected_atoms and not selected_text_items:
            bond_id = next(iter(selected_bonds))
            if bond_id in self.model.bonds:
                bond = self.model.get_bond(bond_id)
                details = {"type": "bond", "id": bond.id, "order": bond.order, "style": bond.style, "aromatic": bond.is_aromatic}
        elif (
            len(self._selected_energy_diagram_items()) == 1
            and not selected_atoms
            and not selected_bonds
            and not selected_text_items
        ):
            item = self._selected_energy_diagram_items()[0]
            details = {
                "type": "energy_diagram",
                "kind": energy_diagram_display_name(item.kind()),
                "label": item.label(),
                "boxes": item.box_count(),
                "occupancies": ", ".join(item.occupancies()),
            }
        elif (
            len(self._selected_semantic_diagram_items()) == 1
            and not selected_atoms
            and not selected_bonds
            and not selected_text_items
        ):
            item = self._selected_semantic_diagram_items()[0]
            details = {
                "type": "semantic_diagram",
                "kind": str(item.semantic_diagram.kind),
                "title": str(item.semantic_diagram.title or ""),
                "levels": len(item.semantic_diagram.levels),
                "lanes": len(item.semantic_diagram.lanes),
            }
        elif len(selected_text_items) == 1 and not selected_atoms and not selected_bonds:
            item = selected_text_items[0]
            cursor = item.textCursor()
            fmt = cursor.charFormat()
            font = QFont(item.font())
            families = fmt.fontFamilies()
            if families:
                font.setFamily(families[0])
            point_size = float(fmt.fontPointSize() or 0.0)
            if point_size > 0.0:
                font.setPointSizeF(point_size)
            weight = int(fmt.fontWeight() or font.weight())
            font.setWeight(weight)
            font.setItalic(bool(fmt.fontItalic()))
            font.setUnderline(bool(fmt.fontUnderline()))
            color = fmt.foreground().color()
            if not color.isValid():
                color = item.defaultTextColor()
            details = {
                "type": "text",
                "font": font,
                "color": color,
                "sub": fmt.verticalAlignment() == QTextCharFormat.VerticalAlignment.AlignSubScript,
                "sup": fmt.verticalAlignment() == QTextCharFormat.VerticalAlignment.AlignSuperScript
            }
            
        self.selection_changed.emit(len(selected_atoms), len(selected_bonds), len(selected_text_items), details)

    def _begin_selection_drag(
        self,
        start_pos: QPointF,
        free_select: bool,
        additive: bool,
    ) -> None:
        """Método auxiliar para  begin selection drag.

        Args:
            start_pos: Descripción del parámetro.
            free_select: Descripción del parámetro.
            additive: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._select_drag_mode = "free" if free_select else "rect"
        self._select_start_pos = start_pos
        self._select_additive = additive
        if free_select:
            path = QPainterPath(start_pos)
            self._select_path = path
            if self._select_preview_path is None:
                item = QGraphicsPathItem()
                pen = QPen(QColor("#4A90D9"), 1.2, Qt.PenStyle.DashLine)
                item.setPen(pen)
                item.setBrush(QBrush(Qt.BrushStyle.NoBrush))
                item.setZValue(40)
                self.scene.addItem(item)
                self._select_preview_path = item
            self._select_preview_path.setPath(path)
            self._select_preview_path.setVisible(True)
        else:
            if self._select_preview_rect is None:
                rect_item = QGraphicsRectItem()
                pen = QPen(QColor("#4A90D9"), 1.2, Qt.PenStyle.DashLine)
                rect_item.setPen(pen)
                rect_item.setBrush(QBrush(Qt.BrushStyle.NoBrush))
                rect_item.setZValue(40)
                self.scene.addItem(rect_item)
                self._select_preview_rect = rect_item
            self._select_preview_rect.setRect(QRectF(start_pos, start_pos))
            self._select_preview_rect.setVisible(True)

    def _update_selection_drag(self, scene_pos: QPointF) -> None:
        """Método auxiliar para  update selection drag.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._select_drag_mode == "free" and self._select_path is not None:
            self._select_path.lineTo(scene_pos)
            if self._select_preview_path is not None:
                self._select_preview_path.setPath(self._select_path)
            return
        if self._select_drag_mode == "rect" and self._select_start_pos is not None:
            if self._select_preview_rect is not None:
                self._select_preview_rect.setRect(QRectF(self._select_start_pos, scene_pos).normalized())

    def _finalize_selection_drag(self) -> None:
        """Método auxiliar para  finalize selection drag.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._select_drag_mode is None:
            return
        items: list = []
        if self._select_drag_mode == "free" and self._select_path is not None:
            path = QPainterPath(self._select_path)
            path.closeSubpath()
            candidates = [
                item
                for item in self.scene.items(path)
                if isinstance(
                    item,
                    (AtomItem, BondItem, ArrowItem, BracketItem, TextAnnotationItem, EnergyDiagramItem, OrbitalAnnotationItem, ImageAnnotationItem),
                )
            ]
            items = [
                item
                for item in candidates
                if path.contains(item.sceneBoundingRect().center())
            ]
        elif self._select_drag_mode == "rect" and self._select_start_pos is not None:
            rect = QRectF(self._select_start_pos, self._last_scene_pos).normalized()
            candidates = [
                item
                for item in self.scene.items(rect)
                if isinstance(
                    item,
                    (AtomItem, BondItem, ArrowItem, BracketItem, TextAnnotationItem, EnergyDiagramItem, OrbitalAnnotationItem, ImageAnnotationItem),
                )
            ]
            items = [
                item
                for item in candidates
                if rect.contains(item.sceneBoundingRect().center())
            ]

        if not self._select_additive:
            self.scene.clearSelection()
        for item in items:
            item.setSelected(True)
        self._sync_selection_from_scene()
        self._clear_selection_drag()

    def _clear_selection_drag(self) -> None:
        """Método auxiliar para  clear selection drag.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._select_drag_mode = None
        self._select_start_pos = None
        self._select_path = None
        self._select_additive = False
        if self._select_preview_path is not None:
            self.scene.removeItem(self._select_preview_path)
            self._select_preview_path = None
        if self._select_preview_rect is not None:
            self.scene.removeItem(self._select_preview_rect)
            self._select_preview_rect = None

    def _selected_text_items_for_opacity(self) -> list[TextAnnotationItem]:
        """Devuelve textos manuales seleccionados o en edición para opacidad."""
        items = [
            item
            for item in self.scene.selectedItems()
            if isinstance(item, TextAnnotationItem) and not bool(item.data(NUMBERING_TEXT_ROLE))
        ]
        focus_item = self._active_text_edit_item()
        if (
            isinstance(focus_item, TextAnnotationItem)
            and not bool(focus_item.data(NUMBERING_TEXT_ROLE))
            and focus_item not in items
        ):
            items.append(focus_item)
        return items

    def _selected_wavy_items(self) -> list[WavyAnchorItem]:
        """Devuelve anclas onduladas seleccionadas."""
        return [item for item in self.scene.selectedItems() if isinstance(item, WavyAnchorItem)]

    def _opacity_target_snapshot(self, *, selected_only: bool) -> dict:
        """Recolecta objetivos de opacidad para selección o documento completo."""
        if selected_only:
            return {
                "atom_ids": sorted(
                    atom_id for atom_id in self.state.selected_atoms if atom_id in self.model.atoms
                ),
                "bond_ids": sorted(
                    bond_id for bond_id in self.state.selected_bonds if bond_id in self.model.bonds
                ),
                "text_items": self._selected_text_items_for_opacity(),
                "arrow_items": self._selected_arrow_items(),
                "bracket_items": self._selected_bracket_items(),
                "energy_diagram_items": self._selected_energy_diagram_items(),
                "semantic_diagram_items": self._selected_semantic_diagram_items(),
                "orbital_items": self._selected_orbital_items(),
                "image_items": self._selected_image_items(),
                "wavy_items": self._selected_wavy_items(),
                "numbering_items": [],
            }
        return {
            "atom_ids": sorted(self.model.atoms.keys()),
            "bond_ids": sorted(self.model.bonds.keys()),
            "text_items": [
                item
                for item in self.scene.items()
                if isinstance(item, TextAnnotationItem) and not bool(item.data(NUMBERING_TEXT_ROLE))
            ],
            "arrow_items": [item for item in self.arrow_items if item.scene() is self.scene],
            "bracket_items": [item for item in self.bracket_items if item.scene() is self.scene],
            "energy_diagram_items": [
                item for item in self.energy_diagram_items if item.scene() is self.scene
            ],
            "semantic_diagram_items": [
                item for item in self.semantic_diagram_items if item.scene() is self.scene
            ],
            "orbital_items": [item for item in self.orbital_items if item.scene() is self.scene],
            "image_items": [item for item in self.image_items if item.scene() is self.scene],
            "wavy_items": [item for item in self._wavy_anchors if item.scene() is self.scene],
            "numbering_items": [
                item for item in self._numbering_overlay_items if item is not None and item.scene() is self.scene
            ],
        }

    @staticmethod
    def _snapshot_has_targets(snapshot: dict) -> bool:
        """Indica si el snapshot contiene al menos un objetivo explícito."""
        return any(
            bool(snapshot.get(key))
            for key in (
                "atom_ids",
                "bond_ids",
                "text_items",
                "arrow_items",
                "bracket_items",
                "energy_diagram_items",
                "semantic_diagram_items",
                "orbital_items",
                "image_items",
                "wavy_items",
            )
        )

    def current_opacity_percent(self) -> int:
        """Devuelve la opacidad a mostrar en el control de transparencia."""
        snapshot = self._opacity_target_snapshot(selected_only=True)
        if not self._snapshot_has_targets(snapshot):
            return int(round(self.canvas_default_opacity() * 100.0))

        values: list[float] = []
        values.extend(self.effective_atom_opacity(atom_id) for atom_id in snapshot["atom_ids"])
        values.extend(self.effective_bond_opacity(bond_id) for bond_id in snapshot["bond_ids"])
        for key in (
            "text_items",
            "arrow_items",
            "bracket_items",
            "energy_diagram_items",
            "semantic_diagram_items",
            "orbital_items",
            "image_items",
            "wavy_items",
        ):
            values.extend(self.effective_item_opacity(item) for item in snapshot[key])
        if not values:
            return int(round(self.canvas_default_opacity() * 100.0))
        return int(round((sum(values) / len(values)) * 100.0))

    def apply_opacity_percent(self, percent: float) -> bool:
        """Aplica una transparencia/opacidad a selección o a todo el documento."""
        target_opacity = normalize_opacity(float(percent) / 100.0)
        selection_snapshot = self._opacity_target_snapshot(selected_only=True)
        selection_mode = self._snapshot_has_targets(selection_snapshot)
        targets = selection_snapshot if selection_mode else self._opacity_target_snapshot(selected_only=False)

        atom_values: dict[int, Optional[float]] = {}
        for atom_id in targets["atom_ids"]:
            atom = self.model.atoms.get(atom_id)
            if atom is None:
                continue
            current_raw = getattr(atom, "opacity", None)
            current_effective = self.effective_atom_opacity(atom_id)
            if selection_mode:
                if self._opacity_equal(current_effective, target_opacity):
                    continue
                atom_values[atom_id] = target_opacity
            elif current_raw is not None or not self._opacity_equal(current_effective, target_opacity):
                atom_values[atom_id] = None

        bond_values: dict[int, Optional[float]] = {}
        for bond_id in targets["bond_ids"]:
            bond = self.model.bonds.get(bond_id)
            if bond is None:
                continue
            current_raw = getattr(bond, "opacity", None)
            current_effective = self.effective_bond_opacity(bond_id)
            if selection_mode:
                if self._opacity_equal(current_effective, target_opacity):
                    continue
                bond_values[bond_id] = target_opacity
            elif current_raw is not None or not self._opacity_equal(current_effective, target_opacity):
                bond_values[bond_id] = None

        item_values: dict[object, Optional[float]] = {}
        for key in (
            "text_items",
            "arrow_items",
            "bracket_items",
            "energy_diagram_items",
            "semantic_diagram_items",
            "orbital_items",
            "image_items",
            "wavy_items",
        ):
            for item in targets[key]:
                current_raw = self.item_raw_opacity(item)
                current_effective = self.effective_item_opacity(item)
                if selection_mode:
                    if self._opacity_equal(current_effective, target_opacity):
                        continue
                    item_values[item] = target_opacity
                elif current_raw is not None or not self._opacity_equal(current_effective, target_opacity):
                    item_values[item] = None

        canvas_opacity = target_opacity if not selection_mode else None
        default_changed = not selection_mode and not self._opacity_equal(
            self.canvas_default_opacity(),
            target_opacity,
        )
        if not atom_values and not bond_values and not item_values and not default_changed:
            return False

        command_text = "Change selection opacity" if selection_mode else "Change canvas opacity"
        self.undo_stack.push(
            ChangeCanvasOpacityCommand(
                self.model,
                self,
                atom_values=atom_values,
                bond_values=bond_values,
                item_values=item_values,
                canvas_opacity=canvas_opacity,
                text=command_text,
            )
        )
        return True

    def _ensure_selection_overlay(self) -> None:
        """Método auxiliar para  ensure selection overlay.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._selection_box is None:
            box = QGraphicsRectItem()
            pen = QPen(QColor("#4A90D9"), 0)
            box.setPen(pen)
            box.setBrush(QBrush(Qt.BrushStyle.NoBrush))
            box.setZValue(45)
            self.scene.addItem(box)
            self._selection_box = box
        if self._selection_handle is None:
            radius = SELECTION_HANDLE_RADIUS_PX
            handle = QGraphicsEllipseItem(-radius, -radius, radius * 2.0, radius * 2.0)
            handle.setBrush(QBrush(QColor("#4A90D9")))
            handle.setPen(QPen(QColor("#4A90D9"), 0))
            handle.setZValue(46)
            handle.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIgnoresTransformations, True)
            self.scene.addItem(handle)
            self._selection_handle = handle

        if self._selection_move_handle is None:
            radius = SELECTION_HANDLE_RADIUS_PX
            handle = QGraphicsEllipseItem(-radius, -radius, radius * 2.0, radius * 2.0)
            handle.setBrush(QBrush(QColor("#FFFFFF")))
            handle.setPen(QPen(QColor("#4A90D9"), 1))
            handle.setZValue(46)
            handle.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIgnoresTransformations, True)
            handle.setCursor(Qt.CursorShape.OpenHandCursor)
            self.scene.addItem(handle)
            self._selection_move_handle = handle

        if self._selection_scale_handle is None:
            size = SELECTION_HANDLE_RADIUS_PX * 2.0
            handle = QGraphicsRectItem(-size / 2.0, -size / 2.0, size, size)
            handle.setBrush(QBrush(QColor("#4A90D9")))
            handle.setPen(QPen(QColor("#4A90D9"), 0))
            handle.setZValue(46)
            handle.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIgnoresTransformations, True)
            handle.setCursor(Qt.CursorShape.SizeFDiagCursor)
            self.scene.addItem(handle)
            self._selection_scale_handle = handle

    def _clear_selection_overlay(self) -> None:
        """Método auxiliar para  clear selection overlay.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._selection_box is not None:
            self.scene.removeItem(self._selection_box)
            self._selection_box = None
        if self._selection_handle is not None:
            self.scene.removeItem(self._selection_handle)
            self._selection_handle = None
        if self._selection_move_handle is not None:
            self.scene.removeItem(self._selection_move_handle)
            self._selection_move_handle = None
        if self._selection_scale_handle is not None:
            self.scene.removeItem(self._selection_scale_handle)
            self._selection_scale_handle = None

    def _selected_items_bbox(self) -> Optional[QRectF]:
        """Método auxiliar para  selected items bbox.

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

        for atom_id in self._selected_atom_ids_for_transform():
            extend_atom_bounds(atom_id)
        for bond_id in self.state.selected_bonds:
            item = self.bond_items.get(bond_id)
            if item is None:
                continue
            extend(item.sceneBoundingRect())
        for item in self.scene.selectedItems():
            # Include TextAnnotationItem
            if isinstance(
                item,
                (
                    ArrowItem,
                    BracketItem,
                    TextAnnotationItem,
                    EnergyDiagramItem,
                    CompositeDiagramItem,
                    OrbitalAnnotationItem,
                    ImageAnnotationItem,
                    WavyAnchorItem,
                    TLCPlateItem,
                    GelElectrophoresisItem,
                ),
            ):
                extend(item.sceneBoundingRect())
        return rect

    def _selected_atom_ids_for_transform(self) -> set[int]:
        """Método auxiliar para  selected atom ids for transform.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atom_ids = set(self.state.selected_atoms)
        for bond_id in self.state.selected_bonds:
            if bond_id in self.model.bonds:
                bond = self.model.get_bond(bond_id)
                atom_ids.add(bond.a1_id)
                atom_ids.add(bond.a2_id)
        return atom_ids

    def _apply_selection_overlay_bbox(self, bbox: Optional[QRectF]) -> None:
        """Aplica el overlay de selección usando un bbox ya calculado."""
        if bbox is None:
            for attr in (
                "_selection_box",
                "_selection_handle",
                "_selection_move_handle",
                "_selection_scale_handle",
            ):
                item = getattr(self, attr)
                if item is None:
                    continue
                try:
                    item.setVisible(False)
                except RuntimeError:
                    setattr(self, attr, None)
            return
        self._ensure_selection_overlay()
        padded = QRectF(bbox)
        pad = max(2.0, float(self.drawing_style.stroke_px))
        padded.adjust(-pad, -pad, pad, pad)

        def offset_in_scene(base: QPointF, dx_view: float, dy_view: float) -> QPointF:
            """Método auxiliar para offset in scene.

            Args:
                base: Descripción del parámetro.
                dx_view: Descripción del parámetro.
                dy_view: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
            view_pt = self.mapFromScene(base)
            view_x = float(view_pt.x()) + dx_view
            view_y = float(view_pt.y()) + dy_view
            view_pt = QPoint(int(round(view_x)), int(round(view_y)))
            return self.mapToScene(view_pt)

        if self._selection_box is not None:
            try:
                self._selection_box.setRect(padded)
                self._selection_box.setVisible(True)
            except RuntimeError:
                self._selection_box = None
        if self._selection_handle is not None:
            top_center = QPointF(padded.center().x(), padded.top())
            handle_pos = offset_in_scene(top_center, 0.0, -SELECTION_ROTATE_OFFSET_PX)
            try:
                self._selection_handle.setPos(handle_pos)
                self._selection_handle.setVisible(True)
            except RuntimeError:
                self._selection_handle = None

        if self._selection_move_handle is not None:
            top_center = QPointF(padded.center().x(), padded.top())
            handle_pos = offset_in_scene(top_center, 0.0, SELECTION_MOVE_OFFSET_PX)
            try:
                self._selection_move_handle.setPos(handle_pos)
                self._selection_move_handle.setVisible(True)
            except RuntimeError:
                self._selection_move_handle = None

        if self._selection_scale_handle is not None:
            corner = QPointF(padded.right(), padded.bottom())
            handle_pos = offset_in_scene(
                corner, -SELECTION_HANDLE_RADIUS_PX, -SELECTION_HANDLE_RADIUS_PX
            )
            try:
                self._selection_scale_handle.setPos(handle_pos)
                self._selection_scale_handle.setVisible(True)
            except RuntimeError:
                self._selection_scale_handle = None

    def _update_selection_overlay(self) -> None:
        """Método auxiliar para  update selection overlay.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self._apply_selection_overlay_bbox(self._selected_items_bbox())

    def _update_drag_selection_overlay(self, delta: QPointF) -> None:
        """Actualiza overlay durante drag trasladando el bbox inicial."""
        bbox = self._drag_start_selection_bbox
        if bbox is None:
            self._update_selection_overlay()
            return
        self._apply_selection_overlay_bbox(bbox.translated(delta))

    def _hit_selection_handle(self, scene_pos: QPointF) -> bool:
        """Método auxiliar para  hit selection handle.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._selection_handle is None:
            return False
        try:
            if not self._selection_handle.isVisible():
                return False
        except RuntimeError:
            self._selection_handle = None
            return False
        return self._hit_handle_item(self._selection_handle, scene_pos)

    def _hit_selection_move_handle(self, scene_pos: QPointF) -> bool:
        """Método auxiliar para  hit selection move handle.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._selection_move_handle is None:
            return False
        try:
            if not self._selection_move_handle.isVisible():
                return False
        except RuntimeError:
            self._selection_move_handle = None
            return False
        return self._hit_handle_item(self._selection_move_handle, scene_pos)

    def _hit_selection_scale_handle(self, scene_pos: QPointF) -> bool:
        """Método auxiliar para  hit selection scale handle.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._selection_scale_handle is None:
            return False
        try:
            if not self._selection_scale_handle.isVisible():
                return False
        except RuntimeError:
            self._selection_scale_handle = None
            return False
        return self._hit_handle_item(self._selection_scale_handle, scene_pos)

    def _hit_handle_item(self, handle: QGraphicsItem, scene_pos: QPointF) -> bool:
        # Use a screen-space hit target so handles stay clickable at low zoom.
        """Método auxiliar para  hit handle item.

        Args:
            handle: Descripción del parámetro.
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        distance_sq = self._handle_item_distance_sq(handle, scene_pos)
        if distance_sq is None:
            return False
        radius = self._handle_item_hit_radius(handle)
        return distance_sq <= (radius * radius)

    def _handle_item_distance_sq(self, handle: QGraphicsItem, scene_pos: QPointF) -> Optional[float]:
        """Calcula distancia cuadrática en pantalla entre el puntero y un handle."""
        view_pos = self.mapFromScene(scene_pos)
        try:
            center_scene = handle.mapToScene(handle.boundingRect().center())
        except RuntimeError:
            return None
        center_view = self.mapFromScene(center_scene)
        dx = float(view_pos.x() - center_view.x())
        dy = float(view_pos.y() - center_view.y())
        return dx * dx + dy * dy

    def _handle_item_hit_radius(self, handle: QGraphicsItem) -> float:
        """Devuelve el radio de click efectivo de un handle en píxeles de vista."""
        try:
            handle_rect_scene = handle.mapToScene(handle.boundingRect()).boundingRect()
            top_left_view = self.mapFromScene(handle_rect_scene.topLeft())
            bottom_right_view = self.mapFromScene(handle_rect_scene.bottomRight())
            visual_radius = max(
                abs(bottom_right_view.x() - top_left_view.x()),
                abs(bottom_right_view.y() - top_left_view.y()),
            ) * 0.75
        except Exception:
            visual_radius = 0.0
        return max(float(visual_radius), SELECTION_HANDLE_RADIUS_PX * 3.0, 18.0)

    def _selection_handle_hit_kind(self, scene_pos: QPointF) -> Optional[str]:
        """Resuelve qué handle de selección está más cerca del puntero."""
        candidates: list[tuple[float, str]] = []
        handles = [
            ("scale", self._selection_scale_handle),
            ("rotate", self._selection_handle),
            ("move", self._selection_move_handle),
        ]
        for kind, handle in handles:
            if handle is None:
                continue
            try:
                if not handle.isVisible():
                    continue
            except RuntimeError:
                continue
            distance_sq = self._handle_item_distance_sq(handle, scene_pos)
            if distance_sq is None:
                continue
            radius = self._handle_item_hit_radius(handle)
            if distance_sq <= (radius * radius):
                candidates.append((distance_sq, kind))
        if not candidates:
            return None
        candidates.sort(key=lambda item: item[0])
        return candidates[0][1]

    def _trackball_atom_ids(self) -> tuple[int, ...]:
        """Devuelve IDs objetivo para trackball (solo selección activa)."""
        atom_ids = self._selected_atom_ids_for_transform()
        return tuple(sorted(atom_id for atom_id in atom_ids if atom_id in self.model.atoms))

    def _clear_trackball_reference(self) -> None:
        """Limpia la referencia pseudo-3D cuando deja de representar la geometría actual."""
        self._rotation_3d_ref_atom_ids = tuple()
        self._rotation_3d_ref_positions = {}
        self._rotation_3d_pitch_deg = 0.0
        self._rotation_3d_yaw_deg = 0.0

    def _clear_trackball_reference_if_desynced(self, atom_ids: set[int]) -> None:
        """Invalida trackball si un movimiento 2D rompe la referencia afín."""
        if self._is_rotating_3d:
            return
        if not atom_ids:
            return
        ref_atom_ids = self._rotation_3d_ref_atom_ids
        if not ref_atom_ids or not self._rotation_3d_ref_positions:
            return
        ref_set = set(ref_atom_ids)
        moved_tracked = ref_set.intersection(atom_ids)
        if not moved_tracked:
            return
        if any(atom_id not in self.model.atoms for atom_id in ref_atom_ids):
            self._clear_trackball_reference()
            return
        if any(atom_id not in self._rotation_3d_ref_positions for atom_id in ref_atom_ids):
            self._clear_trackball_reference()
            return

        projected = self._project_trackball_reference(
            ref_atom_ids,
            self._rotation_3d_pitch_deg,
            self._rotation_3d_yaw_deg,
        )
        for atom_id in moved_tracked:
            atom = self.model.get_atom(atom_id)
            expected_x, expected_y = projected[atom_id]
            if math.hypot(expected_x - atom.x, expected_y - atom.y) > 1e-4:
                self._clear_trackball_reference()
                return

    def _connected_component_atom_ids(self, seed_atom_id: int) -> set[int]:
        """Obtiene la componente conectada del átomo semilla."""
        if seed_atom_id not in self.model.atoms:
            return set()
        visited: set[int] = set()
        stack = [seed_atom_id]
        while stack:
            current = stack.pop()
            if current in visited:
                continue
            visited.add(current)
            for bond in self.model.bonds.values():
                if bond.a1_id == current and bond.a2_id not in visited:
                    stack.append(bond.a2_id)
                elif bond.a2_id == current and bond.a1_id not in visited:
                    stack.append(bond.a1_id)
        return visited

    @staticmethod
    def _signed_angle_delta_deg(start_deg: float, end_deg: float) -> float:
        """Devuelve el delta angular firmado mínimo `end - start`."""
        return (float(end_deg) - float(start_deg) + 180.0) % 360.0 - 180.0

    @staticmethod
    def _rotate_scene_point(point: QPointF, center: QPointF, delta_deg: float) -> QPointF:
        """Rota un punto del lienzo alrededor de un centro."""
        rad = math.radians(float(delta_deg))
        dx = point.x() - center.x()
        dy = point.y() - center.y()
        cos_t = math.cos(rad)
        sin_t = math.sin(rad)
        return QPointF(
            center.x() + dx * cos_t - dy * sin_t,
            center.y() + dx * sin_t + dy * cos_t,
        )

    def _selected_single_bond_id(self) -> Optional[int]:
        """Devuelve el único enlace seleccionado si existe."""
        valid = [bond_id for bond_id in self.state.selected_bonds if bond_id in self.model.bonds]
        if len(valid) != 1:
            return None
        return int(valid[0])

    def fragment_pivot_atom_id(self) -> Optional[int]:
        """Devuelve el átomo pivote de fragmento activo, si sigue existiendo."""
        atom_id = self._fragment_pivot_atom_id
        if atom_id is None:
            return None
        if atom_id not in self.model.atoms:
            self._fragment_pivot_atom_id = None
            return None
        return int(atom_id)

    def set_fragment_pivot_atom(self, atom_id: Optional[int]) -> tuple[bool, str]:
        """Define o limpia el átomo pivote usado para rotar fragmentos."""
        if atom_id is None:
            self._fragment_pivot_atom_id = None
            message = "Átomo pivote de fragmento limpiado."
            self._show_status_message(message)
            return True, message

        pivot_atom_id = int(atom_id)
        if pivot_atom_id not in self.model.atoms:
            self._fragment_pivot_atom_id = None
            message = "El átomo pivote seleccionado ya no existe."
            self._show_status_message(message)
            return False, message

        self._fragment_pivot_atom_id = pivot_atom_id
        atom = self.model.get_atom(pivot_atom_id)
        message = (
            f"Pivote de fragmento fijado en el átomo {pivot_atom_id} ({atom.element}). "
            "Selecciona el fragmento y usa Rotar > Fragmento con pivote."
        )
        self._show_status_message(message)
        return True, message

    def _branch_neighbor_angles_deg(
        self,
        anchor_id: int,
        *,
        exclude_bond_id: Optional[int] = None,
    ) -> list[float]:
        """Ángulos de los enlaces vecinos de `anchor_id`, excluyendo el pivote."""
        if anchor_id not in self.model.atoms:
            return []
        anchor = self.model.get_atom(anchor_id)
        origin = QPointF(anchor.x, anchor.y)
        angles: list[float] = []
        for bond in self.model.bonds.values():
            if exclude_bond_id is not None and bond.id == exclude_bond_id:
                continue
            if bond.a1_id == anchor_id:
                other = self.model.get_atom(bond.a2_id)
            elif bond.a2_id == anchor_id:
                other = self.model.get_atom(bond.a1_id)
            else:
                continue
            angles.append(angle_deg(origin, QPointF(other.x, other.y)))
        return angles

    def _bond_components_without_bond(self, bond_id: int) -> Optional[dict[int, set[int]]]:
        """Parte el grafo por `bond_id`; devuelve `None` si el enlace pertenece a un ciclo."""
        if bond_id not in self.model.bonds:
            return None
        bond = self.model.get_bond(bond_id)
        adjacency: dict[int, list[int]] = {}
        for other in self.model.bonds.values():
            if other.id == bond_id:
                continue
            adjacency.setdefault(other.a1_id, []).append(other.a2_id)
            adjacency.setdefault(other.a2_id, []).append(other.a1_id)

        def bfs(seed_atom_id: int) -> set[int]:
            visited: set[int] = set()
            stack = [seed_atom_id]
            while stack:
                current = stack.pop()
                if current in visited:
                    continue
                visited.add(current)
                for neighbor in adjacency.get(current, []):
                    if neighbor not in visited:
                        stack.append(neighbor)
            return visited

        comp_a = bfs(bond.a1_id)
        if bond.a2_id in comp_a:
            return None
        comp_b = bfs(bond.a2_id)
        return {bond.a1_id: comp_a, bond.a2_id: comp_b}

    def _branch_component_sort_key(self, atom_ids: set[int]) -> tuple[int, int, int]:
        """Orden estable para preferir la rama más pequeña."""
        heavy = 0
        for atom_id in atom_ids:
            atom = self.model.atoms.get(atom_id)
            if atom is not None and atom.element != "H":
                heavy += 1
        smallest_id = min(atom_ids) if atom_ids else 0
        return (heavy, len(atom_ids), smallest_id)

    def _selected_fragment_rotation_context(self, pivot_atom_id: int) -> Optional[dict]:
        """Describe un fragmento seleccionado rotatable alrededor de un átomo pivote."""
        if pivot_atom_id not in self.model.atoms:
            return None
        selected_atom_ids = {
            atom_id for atom_id in self._selected_atom_ids_for_transform()
            if atom_id in self.model.atoms
        }
        if not selected_atom_ids:
            return None

        moving_atom_ids = set(selected_atom_ids)
        moving_atom_ids.discard(int(pivot_atom_id))
        if not moving_atom_ids:
            return None

        pivot_component = self._connected_component_atom_ids(int(pivot_atom_id))
        if not moving_atom_ids.issubset(pivot_component):
            return None

        moving_bond_ids: set[int] = set()
        for bond in self.model.bonds.values():
            a_moving = bond.a1_id in moving_atom_ids
            b_moving = bond.a2_id in moving_atom_ids
            if a_moving and b_moving:
                moving_bond_ids.add(bond.id)
                continue
            if not (a_moving or b_moving):
                continue
            static_atom_id = bond.a2_id if a_moving else bond.a1_id
            if static_atom_id != pivot_atom_id:
                return None
            moving_bond_ids.add(bond.id)

        pivot_atom = self.model.get_atom(int(pivot_atom_id))
        fixed_point = QPointF(pivot_atom.x, pivot_atom.y)
        static_atom_ids = set(self.model.atoms.keys()) - moving_atom_ids
        static_bond_ids = {
            bond.id for bond in self.model.bonds.values() if bond.id not in moving_bond_ids
        }
        return {
            "pivot_atom_id": int(pivot_atom_id),
            "selected_atom_ids": selected_atom_ids,
            "moving_atom_ids": moving_atom_ids,
            "moving_bond_ids": moving_bond_ids,
            "static_atom_ids": static_atom_ids,
            "static_bond_ids": static_bond_ids,
            "fixed_point": fixed_point,
        }

    def _selected_fragment_positions_for_delta(
        self,
        context: dict,
        delta_deg: float,
    ) -> dict[int, tuple[float, float]]:
        """Calcula las posiciones finales de un fragmento seleccionado tras rotación rígida."""
        fixed_point = context["fixed_point"]
        after: dict[int, tuple[float, float]] = {}
        for atom_id in context["moving_atom_ids"]:
            if atom_id not in self.model.atoms:
                continue
            atom = self.model.get_atom(atom_id)
            rotated = self._rotate_scene_point(
                QPointF(atom.x, atom.y),
                fixed_point,
                float(delta_deg),
            )
            after[int(atom_id)] = (rotated.x(), rotated.y())
        return after

    def rotate_selected_fragment_degrees(
        self,
        pivot_atom_id: int,
        delta_deg: float,
    ) -> tuple[bool, str]:
        """Rota rígidamente el fragmento seleccionado alrededor de un átomo pivote fijo."""
        context = self._selected_fragment_rotation_context(int(pivot_atom_id))
        if context is None:
            message = (
                "El fragmento seleccionado debe conectarse al resto solo a través del átomo pivote."
            )
            self._show_status_message(message)
            return False, message
        if abs(float(delta_deg)) <= 1e-6:
            message = "No hay rotación que aplicar al fragmento seleccionado."
            self._show_status_message(message)
            return False, message

        before = {
            atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
            for atom_id in context["moving_atom_ids"]
            if atom_id in self.model.atoms
        }
        after = self._selected_fragment_positions_for_delta(context, float(delta_deg))
        changed = {
            atom_id: coords_before
            for atom_id, coords_before in before.items()
            if atom_id in after and (
                abs(after[atom_id][0] - coords_before[0]) > 1e-6
                or abs(after[atom_id][1] - coords_before[1]) > 1e-6
            )
        }
        if not changed:
            message = "No hay rotación que aplicar al fragmento seleccionado."
            self._show_status_message(message)
            return False, message

        self.undo_stack.push(
            MoveAtomsCommand(
                self.model,
                self,
                changed,
                {atom_id: after[atom_id] for atom_id in changed},
            )
        )
        self.validate_structure()
        message = f"Fragmento girado {float(delta_deg):+.0f}° alrededor del átomo pivote."
        self._show_status_message(message)
        return True, message

    def rotate_selected_fragment_around_pivot(self, delta_deg: float) -> tuple[bool, str]:
        """Rota el fragmento seleccionado usando el pivote activo del canvas."""
        pivot_atom_id = self.fragment_pivot_atom_id()
        if pivot_atom_id is None:
            message = "Define primero un átomo pivote para el fragmento."
            self._show_status_message(message)
            return False, message
        return self.rotate_selected_fragment_degrees(pivot_atom_id, float(delta_deg))

    def _branch_rotation_context(
        self,
        bond_id: int,
        *,
        fixed_atom_id: Optional[int] = None,
        moving_atom_id: Optional[int] = None,
    ) -> Optional[dict]:
        """Describe la rama rotatable alrededor de un enlace pivote."""
        if bond_id not in self.model.bonds:
            return None
        bond = self.model.get_bond(bond_id)
        if bond.style == BondStyle.COORDINATION:
            return None
        components = self._bond_components_without_bond(bond_id)
        if components is None:
            return None

        endpoint_ids = (bond.a1_id, bond.a2_id)
        if fixed_atom_id is not None:
            if fixed_atom_id not in endpoint_ids:
                return None
            moving_atom_id = bond.a2_id if fixed_atom_id == bond.a1_id else bond.a1_id
        elif moving_atom_id is not None:
            if moving_atom_id not in endpoint_ids:
                return None
            fixed_atom_id = bond.a2_id if moving_atom_id == bond.a1_id else bond.a1_id
        else:
            comp_a = components[bond.a1_id]
            comp_b = components[bond.a2_id]
            key_a = self._branch_component_sort_key(comp_a)
            key_b = self._branch_component_sort_key(comp_b)
            if key_a == key_b:
                if self._last_scene_pos is not None:
                    atom_a = self.model.get_atom(bond.a1_id)
                    atom_b = self.model.get_atom(bond.a2_id)
                    dist_a = math.hypot(
                        self._last_scene_pos.x() - atom_a.x,
                        self._last_scene_pos.y() - atom_a.y,
                    )
                    dist_b = math.hypot(
                        self._last_scene_pos.x() - atom_b.x,
                        self._last_scene_pos.y() - atom_b.y,
                    )
                    moving_atom_id = bond.a1_id if dist_a <= dist_b else bond.a2_id
                else:
                    moving_atom_id = min(endpoint_ids)
            else:
                moving_atom_id = bond.a1_id if key_a <= key_b else bond.a2_id
            fixed_atom_id = bond.a2_id if moving_atom_id == bond.a1_id else bond.a1_id

        if fixed_atom_id is None or moving_atom_id is None:
            return None

        fixed_atom = self.model.get_atom(fixed_atom_id)
        moving_atom = self.model.get_atom(moving_atom_id)
        fixed_point = QPointF(fixed_atom.x, fixed_atom.y)
        moving_point = QPointF(moving_atom.x, moving_atom.y)
        moving_atom_ids = set(components[moving_atom_id])
        fixed_component_atom_ids = set(components[fixed_atom_id])

        moving_bond_ids: set[int] = set()
        for other in self.model.bonds.values():
            endpoints = {other.a1_id, other.a2_id}
            if endpoints == {fixed_atom_id, moving_atom_id}:
                moving_bond_ids.add(other.id)
                continue
            if other.a1_id in moving_atom_ids and other.a2_id in moving_atom_ids:
                moving_bond_ids.add(other.id)
        static_atom_ids = set(self.model.atoms.keys()) - moving_atom_ids
        static_bond_ids = {
            other.id for other in self.model.bonds.values() if other.id not in moving_bond_ids
        }

        return {
            "bond_id": int(bond_id),
            "bond": bond,
            "fixed_atom_id": int(fixed_atom_id),
            "moving_atom_id": int(moving_atom_id),
            "moving_atom_ids": moving_atom_ids,
            "fixed_component_atom_ids": fixed_component_atom_ids,
            "static_atom_ids": static_atom_ids,
            "moving_bond_ids": moving_bond_ids,
            "static_bond_ids": static_bond_ids,
            "fixed_point": fixed_point,
            "current_angle_deg": angle_deg(fixed_point, moving_point),
            "pivot_length": math.hypot(moving_point.x() - fixed_point.x(), moving_point.y() - fixed_point.y()),
        }

    def _branch_context_title(self, context: dict) -> str:
        """Etiqueta legible para el lado móvil de una rama."""
        atom = self.model.get_atom(int(context["moving_atom_id"]))
        label = self._editable_atom_label(atom)
        count = len(set(context["moving_atom_ids"]))
        return f"Rama de {label} ({count} át.)"

    def _branch_candidate_angles_deg(self, context: dict, reference_angle_deg: float) -> list[float]:
        """Devuelve candidatos químicos para reorientar una rama."""
        fixed_atom_id = int(context["fixed_atom_id"])
        bond = context["bond"]
        existing_angles = self._branch_neighbor_angles_deg(
            fixed_atom_id,
            exclude_bond_id=int(context["bond_id"]),
        )
        geometry = self._bond_geometry(fixed_atom_id, bond.order, bond.is_aromatic)
        sp3_exact_mode = geometry == "sp3" and len(existing_angles) >= 3
        if sp3_exact_mode:
            candidates = self._sp3_congested_directions_deg(existing_angles)
        else:
            candidates = candidate_directions_deg(geometry, existing_angles, reference_angle_deg)
            if self.state.fixed_angles:
                candidates = self._snap_angles_to_grid(candidates)
        occupied_tolerance = (
            SP3_OCCUPIED_TOLERANCE_DEG if sp3_exact_mode else ANGLE_OCCUPIED_TOLERANCE_DEG
        )
        filtered = filter_occupied_angles_deg(candidates, existing_angles, occupied_tolerance)
        if filtered:
            candidates = filtered
        elif not candidates:
            candidates = [reference_angle_deg]
        candidates.append(float(context["current_angle_deg"]))
        seen: set[float] = set()
        deduped: list[float] = []
        for candidate in candidates:
            normalized = candidate % 360.0
            key = round(normalized, 6)
            if key in seen:
                continue
            seen.add(key)
            deduped.append(normalized)
        return deduped

    def _pick_branch_step_angle_deg(self, context: dict, delta_deg: float) -> float:
        """Avanza la rama al siguiente ángulo químicamente válido en la dirección dada."""
        current = float(context["current_angle_deg"]) % 360.0
        candidates = [
            angle
            for angle in self._branch_candidate_angles_deg(context, current + float(delta_deg))
            if angle_distance_deg(angle, current) > BRANCH_ROTATION_NOOP_TOLERANCE_DEG
        ]
        if not candidates:
            return current
        step_mag = abs(float(delta_deg))
        if delta_deg >= 0:
            def sort_key(angle_value: float) -> tuple[float, float, float]:
                travel = (angle_value - current) % 360.0
                return (
                    travel,
                    abs(travel - step_mag),
                    angle_distance_deg(angle_value, current + delta_deg),
                )
        else:
            def sort_key(angle_value: float) -> tuple[float, float, float]:
                travel = (current - angle_value) % 360.0
                return (
                    travel,
                    abs(travel - step_mag),
                    angle_distance_deg(angle_value, current + delta_deg),
                )
        return min(candidates, key=sort_key)

    def _branch_positions_for_angle(
        self,
        context: dict,
        target_angle_deg: float,
    ) -> tuple[dict[int, tuple[float, float]], float]:
        """Calcula las nuevas coordenadas de la rama móvil para un ángulo objetivo."""
        current = float(context["current_angle_deg"])
        delta_deg = self._signed_angle_delta_deg(current, target_angle_deg)
        center = context["fixed_point"]
        after: dict[int, tuple[float, float]] = {}
        for atom_id in context["moving_atom_ids"]:
            atom = self.model.get_atom(atom_id)
            rotated = self._rotate_scene_point(QPointF(atom.x, atom.y), center, -delta_deg)
            after[int(atom_id)] = (rotated.x(), rotated.y())
        return after, delta_deg

    def _branch_collision_score(
        self,
        context: dict,
        after_positions: dict[int, tuple[float, float]],
        target_angle_deg: float,
    ) -> float:
        """Evalúa qué tan limpia queda una rama respecto al resto del lienzo."""
        moving_atom_ids = set(context["moving_atom_ids"])
        static_atom_ids = set(context["static_atom_ids"])
        moving_bond_ids = set(context["moving_bond_ids"])
        static_bond_ids = set(context["static_bond_ids"])
        fixed_atom_id = int(context["fixed_atom_id"])
        atom_threshold = self.state.bond_length * MIN_ATOM_DIST_SCALE
        bond_threshold = self.state.bond_length * MIN_BOND_DIST_SCALE

        def point_for_atom(atom_id: int) -> QPointF:
            coords = after_positions.get(atom_id)
            if coords is not None:
                return QPointF(coords[0], coords[1])
            atom = self.model.get_atom(atom_id)
            return QPointF(atom.x, atom.y)

        intersections = 0
        atom_penalty = 0.0
        atom_bond_penalty = 0.0
        bond_penalty = 0.0

        for moving_atom_id in moving_atom_ids:
            p_moving = point_for_atom(moving_atom_id)
            for static_atom_id in static_atom_ids:
                if static_atom_id == fixed_atom_id:
                    continue
                p_static = point_for_atom(static_atom_id)
                dist = math.hypot(p_moving.x() - p_static.x(), p_moving.y() - p_static.y())
                if dist < atom_threshold:
                    intersections += 1
                    atom_penalty += atom_threshold - dist

        for moving_bond_id in moving_bond_ids:
            moving_bond = self.model.get_bond(moving_bond_id)
            moved_endpoints = {moving_bond.a1_id, moving_bond.a2_id}
            p0 = point_for_atom(moving_bond.a1_id)
            p1 = point_for_atom(moving_bond.a2_id)
            for static_bond_id in static_bond_ids:
                static_bond = self.model.get_bond(static_bond_id)
                if moved_endpoints.intersection({static_bond.a1_id, static_bond.a2_id}):
                    continue
                p_a = point_for_atom(static_bond.a1_id)
                p_b = point_for_atom(static_bond.a2_id)
                if segments_intersect(p0, p1, p_a, p_b):
                    intersections += 1
                    bond_penalty += atom_threshold
                bond_dist = segment_min_distance(p0, p1, p_a, p_b)
                if bond_dist < bond_threshold:
                    bond_penalty += bond_threshold - bond_dist
                point_dist_a = segment_min_distance(p_a, p_a, p0, p1)
                if point_dist_a < bond_threshold:
                    atom_bond_penalty += bond_threshold - point_dist_a
                point_dist_b = segment_min_distance(p_b, p_b, p0, p1)
                if point_dist_b < bond_threshold:
                    atom_bond_penalty += bond_threshold - point_dist_b

        angular_penalty = angle_distance_deg(target_angle_deg, float(context["current_angle_deg"])) * 0.05
        return (
            intersections * 1000.0
            + atom_penalty * 15.0
            + atom_bond_penalty * 10.0
            + bond_penalty * 8.0
            + angular_penalty
        )

    def _apply_branch_rotation_context(
        self,
        context: Optional[dict],
        *,
        target_angle_deg: float,
        status_message: str,
    ) -> tuple[bool, str]:
        """Aplica una rotación rígida de rama usando undo/redo."""
        if context is None:
            return False, "Selecciona un solo enlace acíclico para reordenar una rama."
        after, delta_deg = self._branch_positions_for_angle(context, target_angle_deg)
        if abs(delta_deg) <= BRANCH_ROTATION_NOOP_TOLERANCE_DEG:
            return False, "No hay una orientación química alternativa para esa rama."
        before = {
            atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
            for atom_id in context["moving_atom_ids"]
            if atom_id in self.model.atoms
        }
        self.undo_stack.push(MoveAtomsCommand(self.model, self, before, after))
        self.validate_structure()
        self._show_status_message(status_message)
        return True, status_message

    def rotate_branch_side_degrees(
        self,
        bond_id: int,
        fixed_atom_id: int,
        delta_deg: float,
    ) -> tuple[bool, str]:
        """Rota una rama alrededor de un enlace usando snap químico."""
        context = self._branch_rotation_context(bond_id, fixed_atom_id=fixed_atom_id)
        if context is None:
            return False, "No se puede rotar una rama alrededor de ese enlace."
        target_angle = self._pick_branch_step_angle_deg(context, delta_deg)
        direction = "+" if float(delta_deg) > 0 else "-"
        return self._apply_branch_rotation_context(
            context,
            target_angle_deg=target_angle,
            status_message=f"Rama girada {direction}{abs(float(delta_deg)):.0f}° con snap químico.",
        )

    def invert_branch_side(self, bond_id: int, fixed_atom_id: int) -> tuple[bool, str]:
        """Invierte una rama 180° exactos alrededor del enlace pivote."""
        context = self._branch_rotation_context(bond_id, fixed_atom_id=fixed_atom_id)
        if context is None:
            return False, "No se puede invertir una rama alrededor de ese enlace."
        target_angle = (float(context["current_angle_deg"]) + 180.0) % 360.0
        return self._apply_branch_rotation_context(
            context,
            target_angle_deg=target_angle,
            status_message="Rama invertida 180°.",
        )

    def auto_arrange_branch_side(self, bond_id: int, fixed_atom_id: int) -> tuple[bool, str]:
        """Elige la mejor orientación local de una rama sin deformar su geometría interna."""
        context = self._branch_rotation_context(bond_id, fixed_atom_id=fixed_atom_id)
        if context is None:
            return False, "No se puede autoacomodar una rama alrededor de ese enlace."
        current_angle = float(context["current_angle_deg"])
        candidates = self._branch_candidate_angles_deg(context, current_angle)
        if not candidates:
            return False, "No hay orientaciones químicas disponibles para esa rama."

        best_angle = current_angle
        best_score = None
        current_score = None
        for angle_value in candidates:
            after_positions, _delta = self._branch_positions_for_angle(context, angle_value)
            score = self._branch_collision_score(context, after_positions, angle_value)
            if angle_distance_deg(angle_value, current_angle) <= BRANCH_ROTATION_NOOP_TOLERANCE_DEG:
                current_score = score
            if best_score is None or score < best_score:
                best_score = score
                best_angle = angle_value

        if angle_distance_deg(best_angle, current_angle) <= BRANCH_ROTATION_NOOP_TOLERANCE_DEG:
            if current_score is None:
                current_score = best_score
            if current_score is None or best_score is None or best_score >= current_score - 1e-6:
                return False, "La rama ya está en la mejor orientación local."

        return self._apply_branch_rotation_context(
            context,
            target_angle_deg=best_angle,
            status_message="Rama autoacomodada sin alterar su geometría interna.",
        )

    def _default_branch_rotation_context(self) -> Optional[dict]:
        """Contexto por defecto para acciones globales de rama: un único enlace seleccionado."""
        bond_id = self._selected_single_bond_id()
        if bond_id is None:
            return None
        return self._branch_rotation_context(bond_id)

    def rotate_selected_branch_degrees(self, delta_deg: float) -> tuple[bool, str]:
        """Rota la rama menor del enlace seleccionado usando snap químico."""
        context = self._default_branch_rotation_context()
        if context is None:
            message = "Selecciona un solo enlace acíclico para reordenar una rama."
            self._show_status_message(message)
            return False, message
        return self.rotate_branch_side_degrees(
            int(context["bond_id"]),
            int(context["fixed_atom_id"]),
            delta_deg,
        )

    def invert_selected_branch(self) -> tuple[bool, str]:
        """Invierte 180° la rama menor del enlace seleccionado."""
        context = self._default_branch_rotation_context()
        if context is None:
            message = "Selecciona un solo enlace acíclico para invertir una rama."
            self._show_status_message(message)
            return False, message
        return self.invert_branch_side(
            int(context["bond_id"]),
            int(context["fixed_atom_id"]),
        )

    def auto_arrange_selected_branch(self) -> tuple[bool, str]:
        """Autoacomoda la rama menor del enlace seleccionado."""
        context = self._default_branch_rotation_context()
        if context is None:
            message = "Selecciona un solo enlace acíclico para autoacomodar una rama."
            self._show_status_message(message)
            return False, message
        return self.auto_arrange_branch_side(
            int(context["bond_id"]),
            int(context["fixed_atom_id"]),
        )

    def _precise_3d_target_atom_ids(
        self,
        clicked_atom_id: Optional[int] = None,
        clicked_bond_id: Optional[int] = None,
    ) -> tuple[int, ...]:
        """Resuelve el conjunto de átomos objetivo (solo selección activa)."""
        selected_atoms = self._selected_atom_ids_for_transform()
        if selected_atoms:
            return tuple(sorted(atom_id for atom_id in selected_atoms if atom_id in self.model.atoms))

        selected_bonds = {
            bond_id for bond_id in self.state.selected_bonds if bond_id in self.model.bonds
        }
        if selected_bonds:
            atom_ids: set[int] = set()
            for bond_id in selected_bonds:
                bond = self.model.get_bond(bond_id)
                if bond.a1_id in self.model.atoms:
                    atom_ids.add(bond.a1_id)
                if bond.a2_id in self.model.atoms:
                    atom_ids.add(bond.a2_id)
            if atom_ids:
                return tuple(sorted(atom_ids))

        return tuple()

    def _ensure_trackball_reference(
        self,
        atom_ids: tuple[int, ...],
        before_positions: Dict[int, Tuple[float, float]],
    ) -> None:
        """Garantiza referencia estable para proyecciones trackball."""
        reset_reference = atom_ids != self._rotation_3d_ref_atom_ids
        if not reset_reference:
            for atom_id in atom_ids:
                if atom_id not in self._rotation_3d_ref_positions:
                    reset_reference = True
                    break
        if not reset_reference:
            projected = self._project_trackball_reference(
                atom_ids,
                self._rotation_3d_pitch_deg,
                self._rotation_3d_yaw_deg,
            )
            max_delta = 0.0
            for atom_id in atom_ids:
                px, py = projected[atom_id]
                cx, cy = before_positions[atom_id]
                max_delta = max(max_delta, math.hypot(px - cx, py - cy))
            if max_delta > TRACKBALL_REFERENCE_MATCH_TOLERANCE_PX:
                # La geometría cambió por otro flujo (mover, editar, undo, etc.).
                reset_reference = True
        if reset_reference:
            self._rotation_3d_ref_atom_ids = atom_ids
            self._rotation_3d_ref_positions = dict(before_positions)
            self._rotation_3d_pitch_deg = 0.0
            self._rotation_3d_yaw_deg = 0.0

    def _project_trackball_reference(
        self,
        atom_ids: tuple[int, ...],
        pitch_deg: float,
        yaw_deg: float,
    ) -> Dict[int, Tuple[float, float]]:
        """Proyecta la referencia 2D con ángulos pseudo-3D absolutos."""
        if not atom_ids:
            return {}
        points = [self._rotation_3d_ref_positions[atom_id] for atom_id in atom_ids]
        count = len(points)
        cx = sum(x for x, _ in points) / count
        cy = sum(y for _, y in points) / count
        rotated = project_3d_rotation(
            points,
            (cx, cy),
            math.radians(pitch_deg),
            math.radians(yaw_deg),
        )
        return {atom_id: point for atom_id, point in zip(atom_ids, rotated)}

    def _begin_3d_rotation_drag(
        self,
        scene_pos: QPointF,
        click_atom_id: Optional[int] = None,
    ) -> bool:
        """Inicia rotación pseudo-3D usando una referencia estable para evitar colapso."""
        atom_ids = self._trackball_atom_ids()
        if not atom_ids:
            return False

        before_positions = {
            atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
            for atom_id in atom_ids
        }
        self._ensure_trackball_reference(atom_ids, before_positions)

        self._is_rotating_3d = True
        self._rotation_3d_start_view_pos = self.mapFromScene(scene_pos)
        self._rotation_3d_before_positions = before_positions
        self._rotation_3d_click_atom_id = click_atom_id
        self._rotation_3d_has_moved = False
        self._rotation_3d_drag_start_pitch_deg = self._rotation_3d_pitch_deg
        self._rotation_3d_drag_start_yaw_deg = self._rotation_3d_yaw_deg
        self.setCursor(Qt.CursorShape.ClosedHandCursor)
        return True

    def _update_3d_rotation_drag(self, scene_pos: QPointF) -> None:
        """Actualiza rotación pseudo-3D con ángulos absolutos sobre referencia fija."""
        if not self._is_rotating_3d or self._rotation_3d_start_view_pos is None:
            return

        atom_ids = self._rotation_3d_ref_atom_ids
        if not atom_ids:
            return

        view_pos = self.mapFromScene(scene_pos)
        dx = float(view_pos.x() - self._rotation_3d_start_view_pos.x())
        dy = float(view_pos.y() - self._rotation_3d_start_view_pos.y())

        pitch_deg = self._rotation_3d_drag_start_pitch_deg + dy * TRACKBALL_ROTATION_DEG_PER_PIXEL
        yaw_deg = self._rotation_3d_drag_start_yaw_deg + dx * TRACKBALL_ROTATION_DEG_PER_PIXEL
        pitch_deg = max(-TRACKBALL_MAX_TILT_DEG, min(TRACKBALL_MAX_TILT_DEG, pitch_deg))
        yaw_deg = max(-TRACKBALL_MAX_TILT_DEG, min(TRACKBALL_MAX_TILT_DEG, yaw_deg))

        projected = self._project_trackball_reference(atom_ids, pitch_deg, yaw_deg)
        changed = False
        moved_atom_ids: set[int] = set()
        for atom_id in atom_ids:
            if atom_id not in self.model.atoms:
                continue
            new_x, new_y = projected[atom_id]
            atom = self.model.get_atom(atom_id)
            if not changed and (abs(new_x - atom.x) > 1e-9 or abs(new_y - atom.y) > 1e-9):
                changed = True
            self.model.update_atom_position(atom_id, new_x, new_y)
            self.update_atom_item(atom_id, new_x, new_y)
            moved_atom_ids.add(atom_id)

        self._rotation_3d_pitch_deg = pitch_deg
        self._rotation_3d_yaw_deg = yaw_deg

        if changed and moved_atom_ids:
            self.update_bond_items_for_atoms(moved_atom_ids)
            self._rotation_3d_has_moved = True

    def _reset_3d_rotation_perspective(self) -> bool:
        """Restaura la perspectiva base del trackball para selección o toda la molécula."""
        selected_atom_ids = self._selected_atom_ids_for_transform()
        if selected_atom_ids:
            atom_ids = tuple(
                sorted(atom_id for atom_id in selected_atom_ids if atom_id in self.model.atoms)
            )
        elif self._rotation_3d_ref_atom_ids:
            atom_ids = self._rotation_3d_ref_atom_ids
        else:
            atom_ids = self._trackball_atom_ids()
        if not atom_ids:
            return False
        if atom_ids != self._rotation_3d_ref_atom_ids:
            return False
        if any(atom_id not in self.model.atoms for atom_id in atom_ids):
            return False
        if any(atom_id not in self._rotation_3d_ref_positions for atom_id in atom_ids):
            return False

        before = {
            atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
            for atom_id in atom_ids
        }
        after = {
            atom_id: self._rotation_3d_ref_positions[atom_id]
            for atom_id in atom_ids
        }
        changed = any(
            abs(before[atom_id][0] - after[atom_id][0]) > 1e-9
            or abs(before[atom_id][1] - after[atom_id][1]) > 1e-9
            for atom_id in atom_ids
        )
        self._rotation_3d_pitch_deg = 0.0
        self._rotation_3d_yaw_deg = 0.0
        self._rotation_3d_drag_start_pitch_deg = 0.0
        self._rotation_3d_drag_start_yaw_deg = 0.0
        if changed:
            self.undo_stack.push(MoveAtomsCommand(self.model, self, before, after))
        else:
            self._update_selection_overlay()
        return True

    def _finalize_3d_rotation_drag(self) -> None:
        """Finaliza la rotación pseudo-3D y registra un único comando de undo."""
        if not self._is_rotating_3d:
            return

        if self._rotation_3d_has_moved and self._rotation_3d_before_positions:
            before: Dict[int, Tuple[float, float]] = {}
            after: Dict[int, Tuple[float, float]] = {}
            changed = False
            for atom_id, (x0, y0) in self._rotation_3d_before_positions.items():
                if atom_id not in self.model.atoms:
                    continue
                atom = self.model.get_atom(atom_id)
                before[atom_id] = (x0, y0)
                after[atom_id] = (atom.x, atom.y)
                if not changed and (abs(atom.x - x0) > 1e-9 or abs(atom.y - y0) > 1e-9):
                    changed = True
            if changed and before:
                self.undo_stack.push(
                    MoveAtomsCommand(
                        self.model,
                        self,
                        before,
                        after,
                        skip_first_redo=True,
                    )
                )
        elif not self._rotation_3d_has_moved and self._rotation_3d_click_atom_id is not None:
            self._cycle_anchor_override(self._rotation_3d_click_atom_id)

        self._is_rotating_3d = False
        self._rotation_3d_start_view_pos = None
        self._rotation_3d_before_positions = {}
        self._rotation_3d_click_atom_id = None
        self._rotation_3d_has_moved = False
        self._rotation_3d_drag_start_pitch_deg = 0.0
        self._rotation_3d_drag_start_yaw_deg = 0.0
        if self.current_tool in {"tool_select", "tool_select_lasso"}:
            self.setCursor(Qt.CursorShape.ArrowCursor)
        self._update_selection_overlay()

    def _prompt_precise_3d_rotation(
        self,
        clicked_atom_id: Optional[int] = None,
        clicked_bond_id: Optional[int] = None,
    ) -> None:
        """Solicita ángulos precisos y aplica rotación pseudo-3D reproducible."""
        atom_ids = self._precise_3d_target_atom_ids(
            clicked_atom_id=clicked_atom_id,
            clicked_bond_id=clicked_bond_id,
        )
        if not atom_ids:
            return
        before = {
            atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
            for atom_id in atom_ids
            if atom_id in self.model.atoms
        }
        if not before:
            return
        atom_ids = tuple(sorted(before.keys()))
        self._ensure_trackball_reference(atom_ids, before)

        default_pitch = self._rotation_3d_pitch_deg
        default_yaw = self._rotation_3d_yaw_deg

        def apply_preview(pitch_deg: float, yaw_deg: float) -> None:
            """Aplica proyección en vivo para previsualización."""
            pitch_deg = max(-TRACKBALL_MAX_TILT_DEG, min(TRACKBALL_MAX_TILT_DEG, float(pitch_deg)))
            yaw_deg = max(-TRACKBALL_MAX_TILT_DEG, min(TRACKBALL_MAX_TILT_DEG, float(yaw_deg)))
            # Mantener estado trackball activo durante preview para que los
            # círculos aromáticos usen el mismo contexto afín que Alt+Drag.
            self._rotation_3d_pitch_deg = pitch_deg
            self._rotation_3d_yaw_deg = yaw_deg
            projected = self._project_trackball_reference(atom_ids, pitch_deg, yaw_deg)
            moved_atom_ids: set[int] = set()
            for atom_id in atom_ids:
                if atom_id not in self.model.atoms:
                    continue
                new_x, new_y = projected[atom_id]
                self.model.update_atom_position(atom_id, new_x, new_y)
                self.update_atom_item(atom_id, new_x, new_y)
                moved_atom_ids.add(atom_id)
            if moved_atom_ids:
                self.update_bond_items_for_atoms(moved_atom_ids)
            else:
                self._update_selection_overlay()

        dialog = TrackballRotationDialog(
            float(default_pitch),
            float(default_yaw),
            TRACKBALL_MAX_TILT_DEG,
            parent=self,
        )
        dialog.preview_changed.connect(apply_preview)
        result = dialog.exec()
        if result == QDialog.DialogCode.Accepted:
            pitch_deg, yaw_deg = dialog.angles()
            self._rotation_3d_pitch_deg = float(pitch_deg)
            self._rotation_3d_yaw_deg = float(yaw_deg)
            after = {
                atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
                for atom_id in atom_ids
                if atom_id in self.model.atoms
            }
            changed_ids = {
                atom_id
                for atom_id in after.keys()
                if (
                    abs(after[atom_id][0] - before[atom_id][0]) > 1e-9
                    or abs(after[atom_id][1] - before[atom_id][1]) > 1e-9
                )
            }
            if not changed_ids:
                self._update_selection_overlay()
                return
            self.undo_stack.push(
                MoveAtomsCommand(
                    self.model,
                    self,
                    {atom_id: before[atom_id] for atom_id in changed_ids},
                    {atom_id: after[atom_id] for atom_id in changed_ids},
                    skip_first_redo=True,
                )
            )
            return

        restored_ids: set[int] = set()
        for atom_id, (x0, y0) in before.items():
            if atom_id not in self.model.atoms:
                continue
            self.model.update_atom_position(atom_id, x0, y0)
            self.update_atom_item(atom_id, x0, y0)
            restored_ids.add(atom_id)
        self._rotation_3d_pitch_deg = default_pitch
        self._rotation_3d_yaw_deg = default_yaw
        if restored_ids:
            self.update_bond_items_for_atoms(restored_ids)
        else:
            self._update_selection_overlay()

    def _begin_rotation_drag(self, scene_pos: QPointF) -> None:
        """Método auxiliar para  begin rotation drag.

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
            and not self._selected_energy_diagram_items()
            and not self._selected_semantic_diagram_items()
            and not self._selected_orbital_items()
            and not self._selected_image_items()
        ):
            return
        bbox = self._selected_items_bbox()
        if bbox is None:
            return
        self._rotation_dragging = True
        self._grab_interaction_mouse()
        self._rotation_center = bbox.center()
        dx = scene_pos.x() - self._rotation_center.x()
        dy = scene_pos.y() - self._rotation_center.y()
        self._rotation_start_angle = math.atan2(dy, dx)
        
        # Capture atoms
        atom_ids = self._selected_atom_ids_for_transform()
        self._rotation_start_positions = {
            atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
            for atom_id in atom_ids
            if atom_id in self.model.atoms
        }
        
        # Capture text items
        self._rotation_start_text_transforms = {} # item -> (pos, rotation)
        for item in self._selected_text_items():
            self._rotation_start_text_transforms[item] = (item.pos(), item.rotation())
        self._rotation_start_arrow_positions = {
            item: (item.start_point(), item.end_point())
            for item in self._selected_arrow_items()
        }
        self._rotation_start_energy_diagram_snapshots = {
            item: self._energy_diagram_transform_snapshot(item)
            for item in self._selected_energy_diagram_items()
        }
        self._rotation_start_semantic_diagram_snapshots = {
            item: self._semantic_diagram_transform_snapshot(item)
            for item in self._selected_semantic_diagram_items()
        }
        self._rotation_start_orbital_snapshots = {
            item: self._orbital_transform_snapshot(item)
            for item in self._selected_orbital_items()
        }
        self._rotation_start_image_snapshots = {
            item: self._image_transform_snapshot(item)
            for item in self._selected_image_items()
        }

    def _update_rotation_drag(self, scene_pos: QPointF) -> None:
        """Método auxiliar para  update rotation drag.

        Args:
            scene_pos: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if self._rotation_center is None or self._rotation_start_angle is None:
            return
        dx = scene_pos.x() - self._rotation_center.x()
        dy = scene_pos.y() - self._rotation_center.y()
        current_angle = math.atan2(dy, dx)
        delta_angle = current_angle - self._rotation_start_angle
        degrees = math.degrees(delta_angle)

        self.rotate_selection_degrees(degrees, use_start_positions=True)

    def _finalize_rotation_drag(self) -> None:
        """Método auxiliar para  finalize rotation drag.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        atom_before = dict(self._rotation_start_positions)
        atom_after = {
            atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
            for atom_id in atom_before
            if atom_id in self.model.atoms
        }
        moved_atoms = {
            atom_id: atom_before[atom_id]
            for atom_id, (after_x, after_y) in atom_after.items()
            if (
                abs(after_x - atom_before[atom_id][0]) > 1e-4
                or abs(after_y - atom_before[atom_id][1]) > 1e-4
            )
        }

        text_before = dict(getattr(self, "_rotation_start_text_transforms", {}))
        text_after = {item: (item.pos(), item.rotation()) for item in text_before}
        moved_text = any(
            not (
                self._point_equal(text_before[item][0], text_after[item][0])
                and abs(text_before[item][1] - text_after[item][1]) <= 1e-4
            )
            for item in text_before
        )

        arrow_before = dict(self._rotation_start_arrow_positions)
        arrow_after = {item: (item.start_point(), item.end_point()) for item in arrow_before}
        moved_arrows = any(
            not (
                self._point_equal(arrow_before[item][0], arrow_after[item][0])
                and self._point_equal(arrow_before[item][1], arrow_after[item][1])
            )
            for item in arrow_before
        )

        energy_before = dict(self._rotation_start_energy_diagram_snapshots)
        energy_after = {item: self._energy_diagram_transform_snapshot(item) for item in energy_before}
        moved_energy_diagrams = any(
            not (
                self._point_equal(energy_before[item][0], energy_after[item][0])
                and abs(energy_before[item][1] - energy_after[item][1]) <= 0.05
                and abs(energy_before[item][2] - energy_after[item][2]) <= 0.05
                and abs(energy_before[item][3] - energy_after[item][3]) <= 0.05
            )
            for item in energy_before
        )

        orbital_before = dict(self._rotation_start_orbital_snapshots)
        orbital_after = {item: self._orbital_transform_snapshot(item) for item in orbital_before}
        moved_orbitals = any(
            not (
                self._point_equal(orbital_before[item][0], orbital_after[item][0])
                and self._point_equal(orbital_before[item][1], orbital_after[item][1])
            )
            for item in orbital_before
        )

        image_before = dict(self._rotation_start_image_snapshots)
        image_after = {item: self._image_transform_snapshot(item) for item in image_before}
        moved_images = any(
            not (
                self._point_equal(image_before[item][0], image_after[item][0])
                and abs(image_before[item][1] - image_after[item][1]) <= 0.05
                and abs(image_before[item][2] - image_after[item][2]) <= 0.05
                and abs(image_before[item][3] - image_after[item][3]) <= 0.05
            )
            for item in image_before
        )

        semantic_before = dict(self._rotation_start_semantic_diagram_snapshots)
        semantic_after = {
            item: self._semantic_diagram_transform_snapshot(item)
            for item in semantic_before
        }
        moved_semantic_diagrams = any(
            not (
                self._point_equal(semantic_before[item][0], semantic_after[item][0])
                and abs(semantic_before[item][1] - semantic_after[item][1]) <= 0.05
                and abs(semantic_before[item][2] - semantic_after[item][2]) <= 0.05
                and abs(semantic_before[item][3] - semantic_after[item][3]) <= 0.05
            )
            for item in semantic_before
        )

        move_count = sum(
            [
                bool(moved_atoms),
                bool(moved_text),
                bool(moved_arrows),
                bool(moved_energy_diagrams),
                bool(moved_semantic_diagrams),
                bool(moved_orbitals),
                bool(moved_images),
            ]
        )
        if move_count > 1:
            self.undo_stack.beginMacro("Rotate selection")

        if moved_atoms:
            self.undo_stack.push(
                MoveAtomsCommand(
                    self.model,
                    self,
                    moved_atoms,
                    {atom_id: atom_after[atom_id] for atom_id in moved_atoms},
                    skip_first_redo=True,
                )
            )

        if moved_text:
            self.undo_stack.push(MoveTextItemsCommand(self, text_before, text_after))

        if moved_arrows:
            self.undo_stack.push(MoveArrowItemsCommand(self, arrow_before, arrow_after))

        if moved_energy_diagrams:
            self.undo_stack.push(
                TransformEnergyDiagramItemsCommand(
                    self,
                    energy_before,
                    energy_after,
                    "Rotate energy diagrams",
                )
            )

        if moved_orbitals:
            self.undo_stack.push(
                TransformOrbitalItemsCommand(
                    self,
                    orbital_before,
                    orbital_after,
                    "Rotate orbitals",
                )
            )

        if moved_semantic_diagrams:
            self.undo_stack.push(
                TransformImageItemsCommand(
                    self,
                    semantic_before,
                    semantic_after,
                    "Rotate semantic diagrams",
                )
            )

        if moved_images:
            self.undo_stack.push(
                TransformImageItemsCommand(
                    self,
                    image_before,
                    image_after,
                    "Rotate images",
                )
            )

        if move_count > 1:
            self.undo_stack.endMacro()

        self._rotation_dragging = False
        self._rotation_center = None
        self._rotation_start_angle = None
        self._rotation_start_positions = {}
        self._rotation_start_text_transforms = {}
        self._rotation_start_arrow_positions = {}
        self._rotation_start_energy_diagram_snapshots = {}
        self._rotation_start_semantic_diagram_snapshots = {}
        self._rotation_start_orbital_snapshots = {}
        self._rotation_start_image_snapshots = {}
        self._release_interaction_mouse()
        self._update_selection_overlay()

    @staticmethod
    def _optional_float_equal(a: Optional[float], b: Optional[float], tol: float = 0.05) -> bool:
        """Compara dos flotantes opcionales con tolerancia visual."""
        if a is None and b is None:
            return True
        if a is None or b is None:
            return False
        return abs(float(a) - float(b)) <= tol

    @staticmethod
    def _point_equal(a: QPointF, b: QPointF, tol: float = 1e-4) -> bool:
        """Compara dos puntos con tolerancia."""
        return abs(a.x() - b.x()) <= tol and abs(a.y() - b.y()) <= tol

    @staticmethod
    def _scale_point_from_anchor(anchor: QPointF, point: QPointF, scale: float) -> QPointF:
        """Escala un punto alrededor de un ancla."""
        return QPointF(
            anchor.x() + (point.x() - anchor.x()) * scale,
            anchor.y() + (point.y() - anchor.y()) * scale,
        )

    def _normalize_label_scale(self, value: float) -> Optional[float]:
        """Normaliza una escala local de etiqueta para herencia/global."""
        scale = max(0.2, float(value))
        return None if abs(scale - 1.0) < 0.02 else scale

    def _normalize_custom_stroke(self, value: float) -> Optional[float]:
        """Convierte un grosor efectivo en override local o herencia."""
        stroke = max(0.6, float(value))
        return None if abs(stroke - float(self.drawing_style.stroke_px)) < 0.05 else stroke

    @staticmethod
    def _text_effective_width(item: TextAnnotationItem) -> float:
        """Devuelve el ancho visual útil actual de un cuadro de texto."""
        width = float(item.textWidth())
        if math.isfinite(width) and width > 0.0:
            return max(1.0, width)
        rect_width = float(item.boundingRect().width())
        if not math.isfinite(rect_width) or rect_width <= 1e-3:
            rect_width = float(item.document().idealWidth())
        return max(1.0, rect_width)

    def _selection_is_text_reflow_only(self) -> bool:
        """Indica si la selección actual solo contiene cuadros de texto libres."""
        return (
            bool(self._selected_text_items())
            and not self._selected_atom_ids_for_transform()
            and not self._selected_arrow_items()
            and not self._selected_bracket_items()
            and not self._selected_energy_diagram_items()
            and not self._selected_semantic_diagram_items()
            and not self._selected_orbital_items()
            and not self._selected_image_items()
        )

    def _selected_bonds_for_scale(self, atom_ids: set[int]) -> list[Bond]:
        """Lista enlaces completamente contenidos en una selección de átomos."""
        if not atom_ids:
            return []
        return [
            bond
            for bond in self.model.bonds.values()
            if bond.a1_id in atom_ids and bond.a2_id in atom_ids
        ]

    def _effective_coordination_sphere_radius(self, atom_id: int) -> float:
        """Obtiene el radio visual efectivo actual de una esfera de coordinación."""
        item = self.atom_items.get(atom_id)
        if item is not None and hasattr(item, "_coordination_draw_radius"):
            try:
                return max(4.0, float(item._coordination_draw_radius()))
            except Exception:
                pass
        atom = self.model.atoms.get(atom_id)
        configured = getattr(atom, "sphere_radius", None) if atom is not None else None
        if configured is not None:
            try:
                return max(4.0, float(configured))
            except Exception:
                return 16.0
        return 16.0

    def _text_scale_snapshot(self, item: TextAnnotationItem) -> tuple[QPointF, float, str, float]:
        """Captura estado geométrico y tipográfico de un texto."""
        width = float(item.textWidth())
        if not math.isfinite(width):
            width = -1.0
        return QPointF(item.pos()), float(item.rotation()), item.font().toString(), width

    def _image_transform_snapshot(self, item: ImageAnnotationItem) -> tuple[QPointF, float, float, float]:
        """Captura posición, tamaño y rotación de una imagen anotada."""
        rect = item.display_rect()
        return QPointF(rect.topLeft()), float(rect.width()), float(rect.height()), float(item.rotation())

    def _energy_diagram_transform_snapshot(self, item: EnergyDiagramItem) -> tuple[QPointF, float, float, float]:
        """Captura posición, tamaño y rotación de un diagrama de energia."""
        rect = item.display_rect()
        return QPointF(rect.topLeft()), float(rect.width()), float(rect.height()), float(item.rotation())

    def _semantic_diagram_transform_snapshot(self, item: CompositeDiagramItem) -> tuple[QPointF, float, float, float]:
        """Captura posición, tamaño y rotación de un diagrama semántico compuesto."""
        rect = item.display_rect()
        return QPointF(rect.topLeft()), float(rect.width()), float(rect.height()), float(item.rotation())

    def _orbital_transform_snapshot(self, item: OrbitalAnnotationItem) -> tuple[QPointF, QPointF]:
        """Captura anchors de un orbital persistente."""
        return item.anchor0(), item.anchor1()

    def _arrow_scale_snapshot(self, item: ArrowItem) -> tuple[QPointF, QPointF, Optional[float]]:
        """Captura geometría y grosor de una flecha."""
        return item.start_point(), item.end_point(), item.stroke_px()

    def _bracket_scale_snapshot(self, item: BracketItem) -> tuple[QRectF, float, Optional[float]]:
        """Captura rectángulo, padding y grosor de un corchete."""
        return item.base_rect(), float(getattr(item, "_padding", 8.0)), item.stroke_px()

    def _capture_scale_state(
        self,
        *,
        atom_ids: Optional[set[int]] = None,
        text_items: Optional[list[TextAnnotationItem]] = None,
        arrow_items: Optional[list[ArrowItem]] = None,
        bracket_items: Optional[list[BracketItem]] = None,
        energy_diagram_items: Optional[list[EnergyDiagramItem]] = None,
        semantic_diagram_items: Optional[list[CompositeDiagramItem]] = None,
        orbital_items: Optional[list[OrbitalAnnotationItem]] = None,
        image_items: Optional[list[ImageAnnotationItem]] = None,
        plate_items: Optional[list[TLCPlateItem | GelElectrophoresisItem]] = None,
    ) -> dict:
        """Captura el estado base usado por el motor de escalado."""
        if atom_ids is None:
            atom_ids = self._selected_atom_ids_for_transform()
        if text_items is None:
            text_items = self._selected_text_items()
        if arrow_items is None:
            arrow_items = self._selected_arrow_items()
        if bracket_items is None:
            bracket_items = self._selected_bracket_items()
        if energy_diagram_items is None:
            energy_diagram_items = self._selected_energy_diagram_items()
        if semantic_diagram_items is None:
            semantic_diagram_items = self._selected_semantic_diagram_items()
        if orbital_items is None:
            orbital_items = self._selected_orbital_items()
        if image_items is None:
            image_items = self._selected_image_items()
        if plate_items is None:
            plate_items = self._selected_plate_items()

        atom_label_scales: Dict[int, Optional[float]] = {}
        atom_sphere_radii: Dict[int, Tuple[Optional[float], float]] = {}
        for atom_id in atom_ids:
            atom = self.model.atoms.get(atom_id)
            if atom is None:
                continue
            atom_label_scales[atom_id] = getattr(atom, "label_scale", None)
            if bool(getattr(atom, "is_coordination_center", False)):
                atom_sphere_radii[atom_id] = (
                    getattr(atom, "sphere_radius", None),
                    self._effective_coordination_sphere_radius(atom_id),
                )

        bond_strokes: Dict[int, Tuple[Optional[float], float]] = {}
        bond_lengths: Dict[int, Optional[float]] = {}
        for bond in self._selected_bonds_for_scale(atom_ids):
            configured = bond.stroke_px
            effective = configured if configured is not None else float(self.drawing_style.stroke_px)
            bond_strokes[bond.id] = (configured, float(effective))
            bond_lengths[bond.id] = bond.length_px

        return {
            "atom_positions": {
                atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
                for atom_id in atom_ids
                if atom_id in self.model.atoms
            },
            "atom_label_scales": atom_label_scales,
            "atom_sphere_radii": atom_sphere_radii,
            "bond_strokes": bond_strokes,
            "bond_lengths": bond_lengths,
            "text_snapshots": {item: self._text_scale_snapshot(item) for item in text_items},
            "text_effective_widths": {
                item: self._text_effective_width(item) for item in text_items
            },
            "arrow_snapshots": {
                item: (
                    item.start_point(),
                    item.end_point(),
                    item.stroke_px(),
                    float(item.stroke_px() if item.stroke_px() is not None else self.drawing_style.stroke_px),
                )
                for item in arrow_items
            },
            "bracket_snapshots": {
                item: (
                    item.base_rect(),
                    float(getattr(item, "_padding", 8.0)),
                    item.stroke_px(),
                    float(item.stroke_px() if item.stroke_px() is not None else self.drawing_style.stroke_px),
                )
                for item in bracket_items
            },
            "energy_diagram_snapshots": {
                item: self._energy_diagram_transform_snapshot(item)
                for item in energy_diagram_items
            },
            "semantic_diagram_snapshots": {
                item: self._semantic_diagram_transform_snapshot(item)
                for item in semantic_diagram_items
            },
            "orbital_snapshots": {
                item: self._orbital_transform_snapshot(item)
                for item in orbital_items
            },
            "image_snapshots": {
                item: self._image_transform_snapshot(item)
                for item in image_items
            },
            "plate_snapshots": {
                item: (item.display_rect(), item.rotation())
                for item in plate_items
            },
        }

    def _apply_scale_state(
        self,
        state: dict,
        anchor: QPointF,
        scale: float,
        *,
        include_style: bool = True,
        text_reflow_only: bool = False,
    ) -> None:
        """Aplica en vivo una escala partiendo de un snapshot inicial."""
        atom_positions = state.get("atom_positions", {})
        atom_ids = set(atom_positions.keys())
        for atom_id, (x0, y0) in atom_positions.items():
            if atom_id not in self.model.atoms:
                continue
            new_pos = self._scale_point_from_anchor(anchor, QPointF(x0, y0), scale)
            self.model.update_atom_position(atom_id, new_pos.x(), new_pos.y())
            self.update_atom_item(atom_id, new_pos.x(), new_pos.y())

        if include_style:
            atom_label_scales = state.get("atom_label_scales", {})
            if atom_label_scales:
                for atom_id, base_scale in atom_label_scales.items():
                    factor = float(base_scale) if base_scale is not None else 1.0
                    self.model.update_atom_label_scale(
                        atom_id,
                        self._normalize_label_scale(factor * scale),
                    )
                self.refresh_label_fonts(atom_label_scales.keys())

            for atom_id, (_configured, effective_radius) in state.get("atom_sphere_radii", {}).items():
                atom = self.model.atoms.get(atom_id)
                if atom is None:
                    continue
                atom.sphere_radius = max(4.0, float(effective_radius) * scale)
                item = self.atom_items.get(atom_id)
                if item is not None and hasattr(item, "refresh_coordination_visual"):
                    item.refresh_coordination_visual()

            for bond_id, (_configured, effective) in state.get("bond_strokes", {}).items():
                if bond_id not in self.model.bonds:
                    continue
                self.model.update_bond(
                    bond_id,
                    stroke_px=self._normalize_custom_stroke(float(effective) * scale),
                )
            for bond_id, base_length in state.get("bond_lengths", {}).items():
                if bond_id not in self.model.bonds or base_length is None:
                    continue
                self.model.update_bond_length(bond_id, max(1.0, float(base_length) * scale))

        if atom_ids:
            self.update_bond_items_for_atoms(atom_ids)
            self.recompute_numbering()

        text_effective_widths = state.get("text_effective_widths", {})
        for item, (pos, rot, font_str, text_width) in state.get("text_snapshots", {}).items():
            item.setPos(self._scale_point_from_anchor(anchor, pos, scale))
            item.setRotation(rot)
            if text_reflow_only:
                base_width = float(text_effective_widths.get(item, 0.0))
                if base_width <= 0.0:
                    base_width = text_width if text_width > 0.0 else self._text_effective_width(item)
                item.setTextWidth(max(1.0, float(base_width) * scale))
            elif include_style:
                font = QFont()
                if font_str:
                    font.fromString(font_str)
                size = font.pointSizeF()
                if size <= 0.0:
                    size = float(font.pointSize()) if font.pointSize() > 0 else 12.0
                font.setPointSizeF(max(1.0, size * scale))
                item.setFont(font)
                item.setTextWidth(text_width * scale if text_width > 0.0 else -1.0)

        for item, (start, end, _configured, effective_stroke) in state.get("arrow_snapshots", {}).items():
            item.update_positions(
                self._scale_point_from_anchor(anchor, start, scale),
                self._scale_point_from_anchor(anchor, end, scale),
            )
            if include_style:
                item.set_stroke_px(self._normalize_custom_stroke(float(effective_stroke) * scale))

        for item, (rect, padding, _configured, effective_stroke) in state.get("bracket_snapshots", {}).items():
            top_left = self._scale_point_from_anchor(anchor, rect.topLeft(), scale)
            bottom_right = self._scale_point_from_anchor(anchor, rect.bottomRight(), scale)
            item.set_rect(QRectF(top_left, bottom_right).normalized(), padding=max(0.0, padding * scale))
            if include_style:
                item.set_stroke_px(self._normalize_custom_stroke(float(effective_stroke) * scale))

        for item, (pos, width, height, rotation) in state.get("energy_diagram_snapshots", {}).items():
            top_left = self._scale_point_from_anchor(anchor, pos, scale)
            item.set_display_rect(
                QRectF(
                    top_left.x(),
                    top_left.y(),
                    max(1.0, float(width) * scale),
                    max(1.0, float(height) * scale),
                )
            )
            item.setRotation(rotation)

        for item, (pos, width, height, rotation) in state.get("semantic_diagram_snapshots", {}).items():
            top_left = self._scale_point_from_anchor(anchor, pos, scale)
            item.set_display_rect(
                QRectF(
                    top_left.x(),
                    top_left.y(),
                    max(1.0, float(width) * scale),
                    max(1.0, float(height) * scale),
                )
            )
            item.setRotation(rotation)

        for item, (anchor0, anchor1) in state.get("orbital_snapshots", {}).items():
            item.set_anchors(
                self._scale_point_from_anchor(anchor, anchor0, scale),
                self._scale_point_from_anchor(anchor, anchor1, scale),
            )

        for item, (pos, width, height, rotation) in state.get("image_snapshots", {}).items():
            top_left = self._scale_point_from_anchor(anchor, pos, scale)
            item.set_display_rect(
                QRectF(
                    top_left.x(),
                    top_left.y(),
                    max(1.0, float(width) * scale),
                    max(1.0, float(height) * scale),
                )
            )
            item.setRotation(rotation)

        for item, (rect0, rotation) in state.get("plate_snapshots", {}).items():
            item.set_display_rect(
                self._scale_rect_from_anchor(anchor, rect0, scale)
            )
            item.setRotation(rotation)

        self._update_selection_overlay()

    def _push_scale_commands(
        self,
        state: dict,
        *,
        macro_name: str,
        skip_first_redo_atoms: bool,
    ) -> bool:
        """Convierte el estado final de una escala en comandos de undo/redo."""
        command_count = 0
        atom_positions_before = state.get("atom_positions", {})
        atom_positions_after = {
            atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
            for atom_id in atom_positions_before
            if atom_id in self.model.atoms
        }
        moved_atoms = {
            atom_id: atom_positions_before[atom_id]
            for atom_id, after in atom_positions_after.items()
            if (
                abs(after[0] - atom_positions_before[atom_id][0]) > 1e-4
                or abs(after[1] - atom_positions_before[atom_id][1]) > 1e-4
            )
        }
        if moved_atoms:
            command_count += 1

        changed_label_scales = []
        for atom_id, before_scale in state.get("atom_label_scales", {}).items():
            atom = self.model.atoms.get(atom_id)
            if atom is None:
                continue
            after_scale = getattr(atom, "label_scale", None)
            if not self._optional_float_equal(before_scale, after_scale, tol=0.02):
                changed_label_scales.append((atom_id, before_scale, after_scale))
        command_count += len(changed_label_scales)

        changed_spheres = []
        for atom_id, (before_radius, _effective) in state.get("atom_sphere_radii", {}).items():
            atom = self.model.atoms.get(atom_id)
            if atom is None:
                continue
            after_radius = getattr(atom, "sphere_radius", None)
            if not self._optional_float_equal(before_radius, after_radius, tol=0.05):
                changed_spheres.append((atom_id, before_radius, after_radius))
        command_count += len(changed_spheres)

        changed_bond_strokes = []
        for bond_id, (before_stroke, _effective) in state.get("bond_strokes", {}).items():
            bond = self.model.bonds.get(bond_id)
            if bond is None:
                continue
            after_stroke = bond.stroke_px
            if not self._optional_float_equal(before_stroke, after_stroke, tol=0.05):
                changed_bond_strokes.append((bond_id, before_stroke, after_stroke))
        command_count += len(changed_bond_strokes)

        changed_bond_lengths = []
        for bond_id, before_length in state.get("bond_lengths", {}).items():
            bond = self.model.bonds.get(bond_id)
            if bond is None:
                continue
            after_length = bond.length_px
            if not self._optional_float_equal(before_length, after_length, tol=0.05):
                changed_bond_lengths.append((bond_id, before_length, after_length))
        command_count += len(changed_bond_lengths)

        text_before = dict(state.get("text_snapshots", {}))
        text_after = {item: self._text_scale_snapshot(item) for item in text_before}
        changed_text = any(
            not (
                self._point_equal(text_before[item][0], text_after[item][0])
                and abs(text_before[item][1] - text_after[item][1]) <= 1e-4
                and text_before[item][2] == text_after[item][2]
                and abs(text_before[item][3] - text_after[item][3]) <= 0.05
            )
            for item in text_before
        )
        if changed_text:
            command_count += 1

        arrow_before = {
            item: (start, end, configured)
            for item, (start, end, configured, _effective) in state.get("arrow_snapshots", {}).items()
        }
        arrow_after = {item: self._arrow_scale_snapshot(item) for item in arrow_before}
        changed_arrows = any(
            not (
                self._point_equal(arrow_before[item][0], arrow_after[item][0])
                and self._point_equal(arrow_before[item][1], arrow_after[item][1])
                and self._optional_float_equal(arrow_before[item][2], arrow_after[item][2], tol=0.05)
            )
            for item in arrow_before
        )
        if changed_arrows:
            command_count += 1

        bracket_before = {
            item: (rect, padding, configured)
            for item, (rect, padding, configured, _effective) in state.get("bracket_snapshots", {}).items()
        }
        bracket_after = {item: self._bracket_scale_snapshot(item) for item in bracket_before}
        changed_brackets = any(
            not (
                self._point_equal(bracket_before[item][0].topLeft(), bracket_after[item][0].topLeft())
                and self._point_equal(bracket_before[item][0].bottomRight(), bracket_after[item][0].bottomRight())
                and abs(bracket_before[item][1] - bracket_after[item][1]) <= 0.05
                and self._optional_float_equal(bracket_before[item][2], bracket_after[item][2], tol=0.05)
            )
            for item in bracket_before
        )
        if changed_brackets:
            command_count += 1

        energy_before = dict(state.get("energy_diagram_snapshots", {}))
        energy_after = {item: self._energy_diagram_transform_snapshot(item) for item in energy_before}
        changed_energy_diagrams = any(
            not (
                self._point_equal(energy_before[item][0], energy_after[item][0])
                and abs(energy_before[item][1] - energy_after[item][1]) <= 0.05
                and abs(energy_before[item][2] - energy_after[item][2]) <= 0.05
                and abs(energy_before[item][3] - energy_after[item][3]) <= 0.05
            )
            for item in energy_before
        )
        if changed_energy_diagrams:
            command_count += 1

        semantic_before = dict(state.get("semantic_diagram_snapshots", {}))
        semantic_after = {
            item: self._semantic_diagram_transform_snapshot(item)
            for item in semantic_before
        }
        changed_semantic_diagrams = any(
            not (
                self._point_equal(semantic_before[item][0], semantic_after[item][0])
                and abs(semantic_before[item][1] - semantic_after[item][1]) <= 0.05
                and abs(semantic_before[item][2] - semantic_after[item][2]) <= 0.05
                and abs(semantic_before[item][3] - semantic_after[item][3]) <= 0.05
            )
            for item in semantic_before
        )
        if changed_semantic_diagrams:
            command_count += 1

        orbital_before = dict(state.get("orbital_snapshots", {}))
        orbital_after = {item: self._orbital_transform_snapshot(item) for item in orbital_before}
        changed_orbitals = any(
            not (
                self._point_equal(orbital_before[item][0], orbital_after[item][0])
                and self._point_equal(orbital_before[item][1], orbital_after[item][1])
            )
            for item in orbital_before
        )
        if changed_orbitals:
            command_count += 1

        image_before = dict(state.get("image_snapshots", {}))
        image_after = {item: self._image_transform_snapshot(item) for item in image_before}
        changed_images = any(
            not (
                self._point_equal(image_before[item][0], image_after[item][0])
                and abs(image_before[item][1] - image_after[item][1]) <= 0.05
                and abs(image_before[item][2] - image_after[item][2]) <= 0.05
                and abs(image_before[item][3] - image_after[item][3]) <= 0.05
            )
            for item in image_before
        )
        if changed_images:
            command_count += 1
            
        plate_before = dict(state.get("plate_snapshots", {}))
        plate_after = {
            item: (item.display_rect(), item.rotation())
            for item in plate_before
        }
        changed_plates = any(
            not (
                self._rect_equal(plate_before[item][0], plate_after[item][0])
                and abs(plate_before[item][1] - plate_after[item][1]) <= 0.05
            )
            for item in plate_before
        )
        if changed_plates:
            command_count += 1

        if command_count <= 0:
            return False

        if command_count > 1:
            self.undo_stack.beginMacro(macro_name)

        if moved_atoms:
            self.undo_stack.push(
                MoveAtomsCommand(
                    self.model,
                    self,
                    moved_atoms,
                    {atom_id: atom_positions_after[atom_id] for atom_id in moved_atoms},
                    skip_first_redo=skip_first_redo_atoms,
                )
            )

        for atom_id, before_scale, after_scale in changed_label_scales:
            self.undo_stack.push(
                ChangeAtomLabelScaleCommand(
                    self.model,
                    self,
                    atom_id,
                    after_scale,
                    old_scale=before_scale,
                )
            )

        for atom_id, before_radius, after_radius in changed_spheres:
            self.undo_stack.push(
                ChangeCoordinationSphereStyleCommand(
                    self.model,
                    self,
                    atom_id,
                    new_radius=after_radius,
                    old_radius=before_radius,
                )
            )

        for bond_id, before_stroke, after_stroke in changed_bond_strokes:
            self.undo_stack.push(
                ChangeBondStrokeCommand(
                    self.model,
                    self,
                    bond_id,
                    after_stroke,
                    old_stroke_px=before_stroke,
                )
            )

        for bond_id, before_length, after_length in changed_bond_lengths:
            self.undo_stack.push(
                ChangeBondLengthCommand(
                    self.model,
                    self,
                    bond_id,
                    after_length,
                    old_length=before_length,
                )
            )

        if changed_text:
            self.undo_stack.push(ScaleTextItemsCommand(self, text_before, text_after))

        if changed_arrows:
            self.undo_stack.push(ScaleArrowItemsCommand(self, arrow_before, arrow_after))

        if changed_brackets:
            self.undo_stack.push(ScaleBracketItemsCommand(self, bracket_before, bracket_after))

        if changed_energy_diagrams:
            self.undo_stack.push(
                TransformEnergyDiagramItemsCommand(
                    self,
                    energy_before,
                    energy_after,
                    "Scale energy diagrams",
                )
            )

        if changed_semantic_diagrams:
            self.undo_stack.push(
                TransformImageItemsCommand(
                    self,
                    semantic_before,
                    semantic_after,
                    "Scale semantic diagrams",
                )
            )

        if changed_orbitals:
            self.undo_stack.push(
                TransformOrbitalItemsCommand(self, orbital_before, orbital_after, "Scale orbitals")
            )

        if changed_images:
            self.undo_stack.push(
                TransformImageItemsCommand(self, image_before, image_after, "Scale images")
            )
            
        if changed_plates:
            self.undo_stack.push(
                TransformPlateItemsCommand(self, plate_before, plate_after, "Scale plates")
            )

        if command_count > 1:
            self.undo_stack.endMacro()
        return True

    def _begin_scale_drag(self, scene_pos: QPointF) -> None:
        """Inicia el escalado interactivo de la selección."""
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
        bbox = self._selected_items_bbox()
        if bbox is None or bbox.width() <= 1e-6 or bbox.height() <= 1e-6:
            return
        anchor = QPointF(bbox.left(), bbox.top())
        handle = QPointF(scene_pos.x(), scene_pos.y())
        start_len = math.hypot(handle.x() - anchor.x(), handle.y() - anchor.y())
        if start_len <= 1e-6:
            return

        self._scale_dragging = True
        self._grab_interaction_mouse()
        self._scale_anchor = anchor
        self._scale_start_handle = handle
        self._scale_start_length = start_len

        state = self._capture_scale_state()
        self._scale_text_reflow_only = self._selection_is_text_reflow_only()
        self._scale_start_positions = state["atom_positions"]
        self._scale_start_atom_label_scales = state["atom_label_scales"]
        self._scale_start_atom_sphere_radii = state["atom_sphere_radii"]
        self._scale_start_bond_strokes = state["bond_strokes"]
        self._scale_start_bond_lengths = state["bond_lengths"]
        self._scale_start_text_positions = {
            item: (snapshot[0], snapshot[1]) for item, snapshot in state["text_snapshots"].items()
        }
        self._scale_start_text_styles = {
            item: (snapshot[2], snapshot[3]) for item, snapshot in state["text_snapshots"].items()
        }
        self._scale_start_text_effective_widths = dict(state.get("text_effective_widths", {}))
        self._scale_start_arrow_positions = {
            item: (snapshot[0], snapshot[1]) for item, snapshot in state["arrow_snapshots"].items()
        }
        self._scale_start_arrow_strokes = {
            item: (snapshot[2], snapshot[3]) for item, snapshot in state["arrow_snapshots"].items()
        }
        self._scale_start_bracket_rects = {
            item: snapshot[0] for item, snapshot in state["bracket_snapshots"].items()
        }
        self._scale_start_bracket_styles = {
            item: (snapshot[1], snapshot[2], snapshot[3])
            for item, snapshot in state["bracket_snapshots"].items()
        }
        self._scale_start_energy_diagram_snapshots = dict(state["energy_diagram_snapshots"])
        self._scale_start_semantic_diagram_snapshots = dict(state["semantic_diagram_snapshots"])
        self._scale_start_orbital_snapshots = dict(state["orbital_snapshots"])
        self._scale_start_image_snapshots = dict(state["image_snapshots"])
        self._scale_has_moved = False

    def _update_scale_drag(self, scene_pos: QPointF) -> None:
        """Actualiza la previsualización de escala interactiva."""
        if not self._scale_dragging or self._scale_anchor is None:
            return
        if self._scale_start_length <= 1e-6:
            return
        dx = scene_pos.x() - self._scale_anchor.x()
        dy = scene_pos.y() - self._scale_anchor.y()
        current_len = math.hypot(dx, dy)
        scale = max(0.05, current_len / self._scale_start_length)
        if abs(scale - 1.0) > 1e-4:
            self._scale_has_moved = True
        self._apply_scale_state(
            {
                "atom_positions": self._scale_start_positions,
                "atom_label_scales": self._scale_start_atom_label_scales,
                "atom_sphere_radii": self._scale_start_atom_sphere_radii,
                "bond_strokes": self._scale_start_bond_strokes,
                "bond_lengths": self._scale_start_bond_lengths,
                "text_snapshots": {
                    item: (
                        pos,
                        rot,
                        self._scale_start_text_styles[item][0],
                        self._scale_start_text_styles[item][1],
                    )
                    for item, (pos, rot) in self._scale_start_text_positions.items()
                },
                "text_effective_widths": dict(self._scale_start_text_effective_widths),
                "arrow_snapshots": {
                    item: (
                        start,
                        end,
                        self._scale_start_arrow_strokes[item][0],
                        self._scale_start_arrow_strokes[item][1],
                    )
                    for item, (start, end) in self._scale_start_arrow_positions.items()
                },
                "bracket_snapshots": {
                    item: (
                        rect,
                        self._scale_start_bracket_styles[item][0],
                        self._scale_start_bracket_styles[item][1],
                        self._scale_start_bracket_styles[item][2],
                    )
                    for item, rect in self._scale_start_bracket_rects.items()
                },
                "energy_diagram_snapshots": dict(self._scale_start_energy_diagram_snapshots),
                "semantic_diagram_snapshots": dict(self._scale_start_semantic_diagram_snapshots),
                "orbital_snapshots": dict(self._scale_start_orbital_snapshots),
                "image_snapshots": dict(self._scale_start_image_snapshots),
            },
            self._scale_anchor,
            scale,
            include_style=True,
            text_reflow_only=self._scale_text_reflow_only,
        )

    def _finalize_scale_drag(self) -> None:
        """Finaliza el escalado interactivo y lo registra en undo/redo."""
        if not self._scale_dragging:
            return

        if self._scale_has_moved:
            self._push_scale_commands(
                {
                    "atom_positions": self._scale_start_positions,
                    "atom_label_scales": self._scale_start_atom_label_scales,
                    "atom_sphere_radii": self._scale_start_atom_sphere_radii,
                    "bond_strokes": self._scale_start_bond_strokes,
                    "bond_lengths": self._scale_start_bond_lengths,
                    "text_snapshots": {
                        item: (
                            pos,
                            rot,
                            self._scale_start_text_styles[item][0],
                            self._scale_start_text_styles[item][1],
                        )
                        for item, (pos, rot) in self._scale_start_text_positions.items()
                    },
                    "text_effective_widths": dict(self._scale_start_text_effective_widths),
                    "arrow_snapshots": {
                        item: (
                            start,
                            end,
                            self._scale_start_arrow_strokes[item][0],
                            self._scale_start_arrow_strokes[item][1],
                        )
                        for item, (start, end) in self._scale_start_arrow_positions.items()
                    },
                    "bracket_snapshots": {
                        item: (
                            rect,
                            self._scale_start_bracket_styles[item][0],
                            self._scale_start_bracket_styles[item][1],
                            self._scale_start_bracket_styles[item][2],
                        )
                        for item, rect in self._scale_start_bracket_rects.items()
                    },
                    "energy_diagram_snapshots": dict(self._scale_start_energy_diagram_snapshots),
                    "semantic_diagram_snapshots": dict(self._scale_start_semantic_diagram_snapshots),
                    "orbital_snapshots": dict(self._scale_start_orbital_snapshots),
                    "image_snapshots": dict(self._scale_start_image_snapshots),
                },
                macro_name="Scale selection",
                skip_first_redo_atoms=True,
            )

        self._scale_dragging = False
        self._scale_anchor = None
        self._scale_start_handle = None
        self._scale_start_length = 0.0
        self._scale_start_positions = {}
        self._scale_start_text_positions = {}
        self._scale_start_text_styles = {}
        self._scale_start_text_effective_widths = {}
        self._scale_start_arrow_positions = {}
        self._scale_start_arrow_strokes = {}
        self._scale_start_bracket_rects = {}
        self._scale_start_bracket_styles = {}
        self._scale_start_energy_diagram_snapshots = {}
        self._scale_start_semantic_diagram_snapshots = {}
        self._scale_start_orbital_snapshots = {}
        self._scale_start_image_snapshots = {}
        self._scale_start_atom_label_scales = {}
        self._scale_start_atom_sphere_radii = {}
        self._scale_start_bond_strokes = {}
        self._scale_start_bond_lengths = {}
        self._scale_text_reflow_only = False
        self._scale_has_moved = False
        self._release_interaction_mouse()
        self._update_selection_overlay()

    def _selected_text_items(self) -> list[TextAnnotationItem]:
        """Método auxiliar para  selected text items.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        return [item for item in self.scene.selectedItems() if isinstance(item, TextAnnotationItem)]

    def _selected_arrow_items(self) -> list[ArrowItem]:
        """Método auxiliar para  selected arrow items.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        return [item for item in self.scene.selectedItems() if isinstance(item, ArrowItem)]

    @staticmethod
    def _is_uniform_line_arrow(item: object) -> bool:
        """Indica si un ArrowItem es una línea recta homogeneizable."""
        return isinstance(item, ArrowItem) and item.kind() in {"line", "line_dashed"}

    def _selected_line_arrow_items(self) -> list[ArrowItem]:
        """Devuelve solo las líneas rectas seleccionadas."""
        return [item for item in self._selected_arrow_items() if self._is_uniform_line_arrow(item)]

    def _selected_bracket_items(self) -> list[BracketItem]:
        """Método auxiliar para  selected bracket items.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        return [item for item in self.scene.selectedItems() if isinstance(item, BracketItem)]

    def _selected_energy_diagram_items(self) -> list[EnergyDiagramItem]:
        """Devuelve diagramas de energia seleccionados."""
        return [
            item
            for item in self.scene.selectedItems()
            if isinstance(item, EnergyDiagramItem)
            and not isinstance(item.parentItem(), CompositeDiagramItem)
        ]

    def _selected_semantic_diagram_items(self) -> list[CompositeDiagramItem]:
        """Devuelve diagramas semánticos compuestos seleccionados."""
        return [item for item in self.scene.selectedItems() if isinstance(item, CompositeDiagramItem)]

    def _selected_orbital_items(self) -> list[OrbitalAnnotationItem]:
        """Devuelve los orbitales persistentes actualmente seleccionados."""
        return [item for item in self.scene.selectedItems() if isinstance(item, OrbitalAnnotationItem)]

    def _selected_image_items(self) -> list[ImageAnnotationItem]:
        """Devuelve las imágenes anotadas actualmente seleccionadas."""
        return [item for item in self.scene.selectedItems() if isinstance(item, ImageAnnotationItem)]

    def _all_scalable_text_items(self) -> list[TextAnnotationItem]:
        """Devuelve textos libres persistentes, excluyendo overlays automáticos."""
        return [
            item
            for item in self.scene.items()
            if isinstance(item, TextAnnotationItem)
            and not bool(item.data(NUMBERING_TEXT_ROLE))
            and item not in self._electron_dots
        ]

    def _targets_bbox(
        self,
        *,
        atom_ids: Iterable[int] = (),
        bond_ids: Iterable[int] = (),
        text_items: Iterable[TextAnnotationItem] = (),
        arrow_items: Iterable[ArrowItem] = (),
        bracket_items: Iterable[BracketItem] = (),
        energy_diagram_items: Iterable[EnergyDiagramItem] = (),
        semantic_diagram_items: Iterable[CompositeDiagramItem] = (),
        orbital_items: Iterable[OrbitalAnnotationItem] = (),
        image_items: Iterable[ImageAnnotationItem] = (),
    ) -> Optional[QRectF]:
        """Calcula un bounding box para un conjunto explícito de elementos."""
        rect: Optional[QRectF] = None

        def extend(candidate: QRectF) -> None:
            nonlocal rect
            if not candidate.isValid() or candidate.isNull():
                return
            rect = candidate if rect is None else rect.united(candidate)

        def extend_atom_bounds(atom_id: int) -> None:
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

        for atom_id in atom_ids:
            extend_atom_bounds(atom_id)
        for bond_id in bond_ids:
            item = self.bond_items.get(bond_id)
            if item is not None:
                extend(item.sceneBoundingRect())
        for item in text_items:
            if item.scene() is self.scene:
                extend(item.sceneBoundingRect())
        for item in arrow_items:
            if item.scene() is self.scene:
                extend(item.sceneBoundingRect())
        for item in bracket_items:
            if item.scene() is self.scene:
                extend(item.sceneBoundingRect())
        for item in energy_diagram_items:
            if item.scene() is self.scene:
                extend(item.sceneBoundingRect())
        for item in semantic_diagram_items:
            if item.scene() is self.scene:
                extend(item.sceneBoundingRect())
        for item in orbital_items:
            if item.scene() is self.scene:
                extend(item.sceneBoundingRect())
        for item in image_items:
            if item.scene() is self.scene:
                extend(item.sceneBoundingRect())
        return rect

    def scale_current_selection(self, scale: float, *, include_style: bool = True) -> bool:
        """Escala la selección actual alrededor de su centro visual."""
        bbox = self._selected_items_bbox()
        if bbox is None:
            return False
        scale_factor = max(0.05, float(scale))
        state = self._capture_scale_state()
        text_reflow_only = self._selection_is_text_reflow_only()
        if not include_style:
            state["atom_label_scales"] = {}
            state["atom_sphere_radii"] = {}
            state["bond_strokes"] = {}
            state["bond_lengths"] = {}
        self._apply_scale_state(
            state,
            bbox.center(),
            scale_factor,
            include_style=include_style,
            text_reflow_only=text_reflow_only,
        )
        return self._push_scale_commands(
            state,
            macro_name="Scale selection",
            skip_first_redo_atoms=True,
        )

    def apply_document_dimensions(
        self,
        *,
        style: DrawingStyle,
        label_font_size: float,
        numbering_font_size: float,
        scale_existing: bool,
        scale_factor: float,
    ) -> None:
        """Aplica dimensiones globales del documento y opcionalmente reescala su contenido."""
        state = None
        anchor = None
        if scale_existing and abs(float(scale_factor) - 1.0) > 1e-4:
            atom_ids = set(self.model.atoms.keys())
            text_items = self._all_scalable_text_items()
            arrow_items = list(self.arrow_items)
            bracket_items = list(self.bracket_items)
            energy_diagram_items = list(self.energy_diagram_items)
            semantic_diagram_items = list(self.semantic_diagram_items)
            orbital_items = list(self.orbital_items)
            image_items = list(self.image_items)
            state = self._capture_scale_state(
                atom_ids=atom_ids,
                text_items=text_items,
                arrow_items=arrow_items,
                bracket_items=bracket_items,
                energy_diagram_items=energy_diagram_items,
                semantic_diagram_items=semantic_diagram_items,
                orbital_items=orbital_items,
                image_items=image_items,
            )
            state["atom_label_scales"] = {}
            anchor_bbox = self._targets_bbox(
                atom_ids=atom_ids,
                bond_ids=self.bond_items.keys(),
                text_items=text_items,
                arrow_items=arrow_items,
                bracket_items=bracket_items,
                energy_diagram_items=energy_diagram_items,
                semantic_diagram_items=semantic_diagram_items,
                orbital_items=orbital_items,
                image_items=image_items,
            )
            anchor = anchor_bbox.center() if anchor_bbox is not None else None

        self.drawing_style = style
        self.state.bond_length = style.bond_length_px
        self.state.label_font_size = max(1.0, float(label_font_size))
        self.state.numbering_font_size = max(1.0, float(numbering_font_size))
        for item in self.atom_items.values():
            item.set_style(style)
        self.refresh_label_fonts()
        self.update_bond_items_for_atoms(set(self.model.atoms.keys()))
        for arrow in self.arrow_items:
            arrow.set_style(style)
        for bracket in self.bracket_items:
            bracket.set_style(style)
        self.refresh_atom_visibility()

        if state is not None and anchor is not None:
            self._apply_scale_state(
                state,
                anchor,
                max(0.05, float(scale_factor)),
                include_style=True,
            )
        self._update_selection_overlay()

    def open_selection_scale_dialog(self) -> None:
        """Abre el diálogo de redimensionado de selección."""
        if self._selected_items_bbox() is None:
            return
        dialog = SelectionScaleDialog(self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        scale_factor, include_style = dialog.values()
        self.scale_current_selection(scale_factor, include_style=include_style)

    def rotate_selection_degrees(self, degrees: float, use_start_positions: bool = False) -> None:
        """Rotate selected items by degrees around their collective center."""
        selected_atom_ids = self._selected_atom_ids_for_transform()
        selected_text_items = self._selected_text_items()
        selected_arrow_items = self._selected_arrow_items()
        selected_energy_diagram_items = self._selected_energy_diagram_items()
        selected_semantic_diagram_items = self._selected_semantic_diagram_items()
        selected_orbital_items = self._selected_orbital_items()
        selected_image_items = self._selected_image_items()
        
        if (
            not selected_atom_ids
            and not selected_text_items
            and not selected_arrow_items
            and not selected_energy_diagram_items
            and not selected_semantic_diagram_items
            and not selected_orbital_items
            and not selected_image_items
        ):
            return

        cx, cy = 0.0, 0.0
        
        # If dragging, usage rotation center
        if use_start_positions and self._rotation_center is not None:
             cx = self._rotation_center.x()
             cy = self._rotation_center.y()
             # We use the fixed center from drag start
        else:
             # Calculate collective center
            bbox = self._selected_items_bbox()
            if bbox is None:
                return
            center = bbox.center()
            cx, cy = center.x(), center.y()
            
        rad = math.radians(degrees)
        cos_a = math.cos(rad)
        sin_a = math.sin(rad)

        def rotate_point(pt: QPointF) -> QPointF:
            dx = pt.x() - cx
            dy = pt.y() - cy
            return QPointF(dx * cos_a - dy * sin_a + cx, dx * sin_a + dy * cos_a + cy)

        # Rotate Atoms
        if selected_atom_ids:
            new_positions = {}
            for aid in selected_atom_ids:
                if use_start_positions and aid in self._rotation_start_positions:
                    ox, oy = self._rotation_start_positions[aid]
                else:
                    atom = self.model.get_atom(aid)
                    ox, oy = atom.x, atom.y
                
                dx = ox - cx
                dy = oy - cy
                nx = dx * cos_a - dy * sin_a + cx
                ny = dx * sin_a + dy * cos_a + cy
                new_positions[aid] = (nx, ny)

            # If this is live drag, update directly without undo stack (handled at end)
            # Actually, canvas usually updates model directly during drag and pushes single command at end?
            # Existing _update_rotation_drag called rotate_selection_degrees.
            # But the original implementation of rotate_selection_degrees PUSHED COMMANDS.
            # We should probably separate "preview" from "commit".
            # For this refactor, let's assume this method updates positions efficiently.
            
            # Use MoveAtomsCommand for undo stack if valid drop?
            # Or just update positions if use_start_positions is True (dragging).
            if use_start_positions:
                for aid, (nx, ny) in new_positions.items():
                    self.model.update_atom_position(aid, nx, ny)
                    self.update_atom_item(aid, nx, ny)
                self.update_bond_items_for_atoms(selected_atom_ids)
            else:
                before = {
                    atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
                    for atom_id in selected_atom_ids
                }
                after = new_positions
                self.undo_stack.push(
                    MoveAtomsCommand(self.model, self, before, after)
                )

        # Rotate Text Items
        if selected_text_items:
             for item in selected_text_items:
                if use_start_positions and hasattr(self, "_rotation_start_text_transforms") and item in self._rotation_start_text_transforms:
                    start_pos, start_rot = self._rotation_start_text_transforms[item]
                    
                    # Logic: We are rotating around (cx, cy) by 'degrees' relative to START pose.
                    # QGraphicsItem rotation is around its transform origin.
                    
                    # 1. Calculate new center position (orbit)
                    # Item's scene pos is top-left usually? TextItem might be different.
                    # Let's use mapToScene(center) if we stored it?
                    # Getting intricate. Simplified:
                    # Treat item.pos() as the anchor.
                    
                    # Better: define item center relative to scene.
                    # Start Logic:
                    # We stored start_pos (item.pos) and start_rot.
                    
                    # We need the item's center at start state.
                    # But calculating that is hard if we don't store it.
                    # Let's assume start_pos is good enough for now, or improve storage.
                    
                    # Let's do delta rotation on CURRENT state? 
                    # use_start_positions implies absolute rotation from start.
                    
                    # Re-implement correctly for visual drag:
                    
                    # Calculate the item's center in scene coordinates based on its START position and rotation
                    # This is complex because QGraphicsTextItem's bounding rect changes with rotation.
                    # For simplicity during drag, let's apply the rotation to the item's current state
                    # relative to the rotation center, and then reset its rotation.
                    
                    # This is a tricky part. The current implementation of _update_rotation_drag
                    # calls this function with the *total* delta_angle from the start.
                    # So, we need to apply the rotation from the *original* state.
                    
                    # Reset to start state
                    item.setPos(start_pos)
                    item.setRotation(start_rot)
                    
                    # Now apply the full rotation from the start state
                    # Set origin to center for rotation
                    br = item.boundingRect()
                    center = br.center()
                    item.setTransformOriginPoint(center)
                    
                    # 1. Rotate in place
                    item.setRotation(start_rot + degrees) # Apply total rotation
                    
                    # 2. Orbit if we are rotating around a collective center that isn't just this item's center
                    item_scene_center = item.mapToScene(center)
                    # Only orbit if the item's center is not the rotation center (i.e., it's not a single item rotation around its own center)
                    if math.hypot(item_scene_center.x() - cx, item_scene_center.y() - cy) > 1.0:
                        # Calculate new center position after orbiting
                        dx_orbit = item_scene_center.x() - cx
                        dy_orbit = item_scene_center.y() - cy
                        nx_center_orbit = dx_orbit * cos_a - dy_orbit * sin_a + cx
                        ny_center_orbit = dx_orbit * sin_a + dy_orbit * cos_a + cy
                        
                        # Move item so its center is at nx_center_orbit, ny_center_orbit
                        offset = QPointF(nx_center_orbit, ny_center_orbit) - item_scene_center
                        item.moveBy(offset.x(), offset.y())
                else:
                    # Increment rotation (toolbar)
                    current_rot = item.rotation()
                    # Rotate around center of bbox? Or item center?
                    # Users usually expect in-place rotation if handled via toolbar unless group.
                    # For consistency with atoms, rotate around group center.
                    
                    item_center = item.sceneBoundingRect().center()
                    dx = item_center.x() - cx
                    dy = item_center.y() - cy
                    
                    nx_center = dx * cos_a - dy * sin_a + cx
                    ny_center = dx * sin_a + dy * cos_a + cy
                    
                    # Move item center to new position
                    # Diff
                    offset = QPointF(nx_center, ny_center) - item_center
                    item.moveBy(offset.x(), offset.y())
                    
                    # Add rotation
                    item.setTransformOriginPoint(item.boundingRect().center())
                    item.setRotation(current_rot + degrees)

        # Rotate Arrow Items
        if selected_arrow_items:
            if use_start_positions:
                for item in selected_arrow_items:
                    if item in self._rotation_start_arrow_positions:
                        start_pt, end_pt = self._rotation_start_arrow_positions[item]
                    else:
                        start_pt, end_pt = item.start_point(), item.end_point()
                    item.update_positions(rotate_point(start_pt), rotate_point(end_pt))
            else:
                before = {
                    item: (item.start_point(), item.end_point())
                    for item in selected_arrow_items
                }
                after = {
                    item: (rotate_point(start), rotate_point(end))
                    for item, (start, end) in before.items()
                }
                if before:
                    self.undo_stack.push(MoveArrowItemsCommand(self, before, after))

        if selected_energy_diagram_items:
            if use_start_positions:
                for item in selected_energy_diagram_items:
                    if item in self._rotation_start_energy_diagram_snapshots:
                        start_pos, width, height, start_rotation = self._rotation_start_energy_diagram_snapshots[item]
                    else:
                        start_pos, width, height, start_rotation = self._energy_diagram_transform_snapshot(item)
                    start_center = QPointF(start_pos.x() + width / 2.0, start_pos.y() + height / 2.0)
                    rotated_center = rotate_point(start_center)
                    item.set_display_rect(
                        QRectF(
                            rotated_center.x() - width / 2.0,
                            rotated_center.y() - height / 2.0,
                            width,
                            height,
                        )
                    )
                    item.setRotation(start_rotation + degrees)
            else:
                before = {
                    item: self._energy_diagram_transform_snapshot(item)
                    for item in selected_energy_diagram_items
                }
                after = {}
                for item, (pos, width, height, rotation) in before.items():
                    center_point = QPointF(pos.x() + width / 2.0, pos.y() + height / 2.0)
                    rotated_center = rotate_point(center_point)
                    after[item] = (
                        QPointF(rotated_center.x() - width / 2.0, rotated_center.y() - height / 2.0),
                        width,
                        height,
                        rotation + degrees,
                    )
                if before:
                    self.undo_stack.push(
                        TransformEnergyDiagramItemsCommand(
                            self,
                            before,
                            after,
                            "Rotate energy diagrams",
                        )
                    )

        if selected_semantic_diagram_items:
            if use_start_positions:
                for item in selected_semantic_diagram_items:
                    if item in self._rotation_start_semantic_diagram_snapshots:
                        start_pos, width, height, start_rotation = self._rotation_start_semantic_diagram_snapshots[item]
                    else:
                        start_pos, width, height, start_rotation = self._semantic_diagram_transform_snapshot(item)
                    start_center = QPointF(start_pos.x() + width / 2.0, start_pos.y() + height / 2.0)
                    rotated_center = rotate_point(start_center)
                    item.set_display_rect(
                        QRectF(
                            rotated_center.x() - width / 2.0,
                            rotated_center.y() - height / 2.0,
                            width,
                            height,
                        )
                    )
                    item.setRotation(start_rotation + degrees)
            else:
                before = {
                    item: self._semantic_diagram_transform_snapshot(item)
                    for item in selected_semantic_diagram_items
                }
                after = {}
                for item, (pos, width, height, rotation) in before.items():
                    center_point = QPointF(pos.x() + width / 2.0, pos.y() + height / 2.0)
                    rotated_center = rotate_point(center_point)
                    after[item] = (
                        QPointF(rotated_center.x() - width / 2.0, rotated_center.y() - height / 2.0),
                        width,
                        height,
                        rotation + degrees,
                    )
                if before:
                    self.undo_stack.push(
                        TransformImageItemsCommand(self, before, after, "Rotate semantic diagrams")
                    )

        if selected_orbital_items:
            if use_start_positions:
                for item in selected_orbital_items:
                    if item in self._rotation_start_orbital_snapshots:
                        start_anchor0, start_anchor1 = self._rotation_start_orbital_snapshots[item]
                    else:
                        start_anchor0, start_anchor1 = self._orbital_transform_snapshot(item)
                    item.set_anchors(rotate_point(start_anchor0), rotate_point(start_anchor1))
            else:
                before = {item: self._orbital_transform_snapshot(item) for item in selected_orbital_items}
                after = {
                    item: (rotate_point(anchor0), rotate_point(anchor1))
                    for item, (anchor0, anchor1) in before.items()
                }
                if before:
                    self.undo_stack.push(
                        TransformOrbitalItemsCommand(self, before, after, "Rotate orbitals")
                    )

        if selected_image_items:
            if use_start_positions:
                for item in selected_image_items:
                    if item in self._rotation_start_image_snapshots:
                        start_pos, width, height, start_rotation = self._rotation_start_image_snapshots[item]
                    else:
                        start_pos, width, height, start_rotation = self._image_transform_snapshot(item)
                    start_center = QPointF(start_pos.x() + width / 2.0, start_pos.y() + height / 2.0)
                    rotated_center = rotate_point(start_center)
                    item.set_display_rect(
                        QRectF(
                            rotated_center.x() - width / 2.0,
                            rotated_center.y() - height / 2.0,
                            width,
                            height,
                        )
                    )
                    item.setRotation(start_rotation + degrees)
            else:
                before = {
                    item: self._image_transform_snapshot(item)
                    for item in selected_image_items
                }
                after = {}
                for item, (pos, width, height, rotation) in before.items():
                    center_point = QPointF(pos.x() + width / 2.0, pos.y() + height / 2.0)
                    rotated_center = rotate_point(center_point)
                    after[item] = (
                        QPointF(rotated_center.x() - width / 2.0, rotated_center.y() - height / 2.0),
                        width,
                        height,
                        rotation + degrees,
                    )
                if before:
                    self.undo_stack.push(
                        TransformImageItemsCommand(self, before, after, "Rotate images")
                    )

        # Update overlay at the end
        if use_start_positions:
             self._update_selection_overlay()

    def flip_selection_horizontal(self) -> None:
        """Método auxiliar para flip selection horizontal.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if not self._selected_atom_ids_for_transform() and not self._selected_text_items():
            return
        bbox = self._selected_items_bbox()
        if bbox is None:
            return
        cx = bbox.center().x()
        
        # Atoms
        atom_ids = self._selected_atom_ids_for_transform()
        # ... existing atom flip logic ...
        if atom_ids:
            before = {
                atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
                for atom_id in atom_ids
                if atom_id in self.model.atoms
            }
            after = {
                atom_id: (cx - (x0 - cx), y0)
                for atom_id, (x0, y0) in before.items()
            }
            cmd = MoveAtomsCommand(self.model, self, before, after)
            self.undo_stack.push(cmd)
            
        # Text Items
        # Flip position relative to cx
        for item in self._selected_text_items():
            # Flip center X
            br = item.sceneBoundingRect()
            center = br.center()
            new_cx = cx - (center.x() - cx)
            # Move
            item.moveBy(new_cx - center.x(), 0)
            # Potentially flip content? Usually text isn't mirrored.
            
        self._update_selection_overlay()

    def flip_selection_vertical(self) -> None:
        """Método auxiliar para flip selection vertical.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if not self._selected_atom_ids_for_transform() and not self._selected_text_items():
            return
        bbox = self._selected_items_bbox()
        if bbox is None:
            return
        cy = bbox.center().y()
        atom_ids = self._selected_atom_ids_for_transform()
        before = {
            atom_id: (self.model.get_atom(atom_id).x, self.model.get_atom(atom_id).y)
            for atom_id in atom_ids
            if atom_id in self.model.atoms
        }
        after = {
            atom_id: (x0, cy - (y0 - cy))
            for atom_id, (x0, y0) in before.items()
        }
        cmd = MoveAtomsCommand(self.model, self, before, after)
        self.undo_stack.push(cmd)
        self._update_selection_overlay()
