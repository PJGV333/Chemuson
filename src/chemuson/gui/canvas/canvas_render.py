from __future__ import annotations

from ._shared import (
    AROMATIC_CIRCLE_ATOMS_ROLE,
    AddImageItemCommand,
    AromaticCircleItem,
    ArrowItem,
    AtomItem,
    Bond,
    BondItem,
    BondStyle,
    BracketItem,
    CompositeDiagramItem,
    DrawingStyle,
    EnergyDiagramItem,
    GelBandItem,
    GelElectrophoresisItem,
    ITEM_OPACITY_ROLE,
    ITEM_OPACITY_UNSET,
    ImageAnnotationItem,
    Iterable,
    MolView,
    NUMBERING_TEXT_ROLE,
    Optional,
    OrbitalAnnotationItem,
    PAPER_MARGIN,
    PASTE_IMAGE_MAX_PAPER_FRACTION,
    PASTE_IMAGE_OFFSET_PX,
    QBuffer,
    QColor,
    QGraphicsItem,
    QImage,
    QMimeData,
    QPainter,
    QPixmap,
    QPointF,
    QRect,
    QRectF,
    QSize,
    Qt,
    SUPPORTED_IMAGE_FILE_MIME_TYPES,
    TLCPlateItem,
    TLCSpotItem,
    TRACKBALL_REFERENCE_MATCH_TOLERANCE_PX,
    TextAnnotationItem,
    WavyAnchorItem,
    angle_deg,
    find_rings_simple,
    math,
    normalize_opacity,
    normalize_optional_opacity,
    os,
    ring_bonds,
)

class CanvasRenderMixin:
    def set_canvas_opacity_default(self, opacity: float) -> None:
        """Define la opacidad por defecto del documento."""
        self.state.canvas_opacity = normalize_opacity(opacity)

    def canvas_default_opacity(self) -> float:
        """Devuelve la opacidad por defecto efectiva del documento."""
        return normalize_opacity(getattr(self.state, "canvas_opacity", 1.0))

    @staticmethod
    def _opacity_equal(left: object, right: object, tol: float = 0.004) -> bool:
        """Compara opacidades normalizadas con una tolerancia pequeña."""
        try:
            return abs(normalize_opacity(left) - normalize_opacity(right)) <= float(tol)
        except Exception:
            return False

    def item_raw_opacity(self, item: QGraphicsItem) -> Optional[float]:
        """Lee la opacidad local persistida de un item o `None` si hereda."""
        if item is None:
            return None
        try:
            return normalize_optional_opacity(item.data(ITEM_OPACITY_ROLE))
        except Exception:
            return None

    def effective_item_opacity(self, item: QGraphicsItem) -> float:
        """Resuelve la opacidad efectiva de un item gráfico."""
        raw = self.item_raw_opacity(item)
        if raw is None:
            return self.canvas_default_opacity()
        return normalize_opacity(raw)

    def effective_atom_opacity(self, atom_id: int | AtomItem | object) -> float:
        """Resuelve la opacidad efectiva de un átomo."""
        atom = atom_id
        if isinstance(atom_id, AtomItem):
            atom = getattr(atom_id, "atom", None)
        elif not hasattr(atom_id, "opacity"):
            atom = self.model.atoms.get(int(atom_id))
        raw = getattr(atom, "opacity", None) if atom is not None else None
        if raw is None:
            return self.canvas_default_opacity()
        return normalize_opacity(raw)

    def effective_bond_opacity(self, bond_id: int | BondItem | object) -> float:
        """Resuelve la opacidad efectiva de un enlace."""
        bond = bond_id
        if isinstance(bond_id, BondItem):
            bond = self.model.bonds.get(bond_id.bond_id)
        elif not hasattr(bond_id, "opacity"):
            bond = self.model.bonds.get(int(bond_id))
        raw = getattr(bond, "opacity", None) if bond is not None else None
        if raw is None:
            return self.canvas_default_opacity()
        return normalize_opacity(raw)

    def set_graphics_item_opacity(
        self,
        item: QGraphicsItem,
        opacity: Optional[float] | object,
    ) -> None:
        """Fija una opacidad local persistente para un item gráfico."""
        if item is None:
            return
        raw_opacity = normalize_optional_opacity(opacity)
        try:
            item.setData(ITEM_OPACITY_ROLE, raw_opacity)
            item.setOpacity(
                self.canvas_default_opacity()
                if raw_opacity is None
                else normalize_opacity(raw_opacity)
            )
        except RuntimeError:
            return

    def ensure_graphics_item_opacity(
        self,
        item: QGraphicsItem,
        opacity: Optional[float] | object = ITEM_OPACITY_UNSET,
    ) -> None:
        """Aplica opacidad efectiva preservando la local existente si ya existe."""
        if opacity is ITEM_OPACITY_UNSET:
            opacity = self.item_raw_opacity(item)
        self.set_graphics_item_opacity(item, opacity)

    def update_atom_item_opacity(self, atom_id: int) -> None:
        """Sincroniza la opacidad efectiva de un átomo ya dibujado."""
        item = self.atom_items.get(atom_id)
        if item is None:
            return
        try:
            item.setOpacity(self.effective_atom_opacity(atom_id))
        except RuntimeError:
            return

    def refresh_numbering_opacity(self) -> None:
        """Reaplica la opacidad global a los overlays de numeración."""
        opacity = self.canvas_default_opacity()
        for item in list(self._numbering_overlay_items):
            if item is None or item.scene() is not self.scene:
                continue
            try:
                item.setOpacity(opacity)
            except RuntimeError:
                continue

    def _aromatic_circle_effective_opacity(self, atom_ids: Iterable[int]) -> float:
        """Calcula una opacidad efectiva para el círculo aromático de un anillo."""
        ring_set = {int(atom_id) for atom_id in atom_ids}
        if not ring_set:
            return self.canvas_default_opacity()
        opacities: list[float] = []
        for bond in self.model.bonds.values():
            if not bond.is_aromatic:
                continue
            if bond.a1_id in ring_set and bond.a2_id in ring_set:
                opacities.append(self.effective_bond_opacity(bond))
        if not opacities:
            return self.canvas_default_opacity()
        return sum(opacities) / len(opacities)

    def refresh_aromatic_circle_opacities(self) -> None:
        """Reaplica la opacidad efectiva de los círculos aromáticos visibles."""
        for circle in list(self.aromatic_circles):
            if circle is None or circle.scene() is not self.scene:
                continue
            atom_ids = circle.data(AROMATIC_CIRCLE_ATOMS_ROLE)
            if not atom_ids:
                continue
            try:
                circle.setOpacity(self._aromatic_circle_effective_opacity(atom_ids))
            except RuntimeError:
                continue

    def _refresh_scene_after_image_insert(self) -> None:
        """Sincroniza visuales químicos tras pegar imágenes sobre un canvas ocupado."""
        if not self.model.atoms and not self.model.bonds:
            self._update_hover(self._last_scene_pos)
            self._update_selection_overlay()
            self.scene.invalidate()
            self.viewport().update()
            return
        snapshot = self._selection_snapshot()
        self._sync_scene_with_model()
        self._restore_selection_snapshot(snapshot)
        self.validate_structure()
        self._update_hover(self._last_scene_pos)
        self._update_selection_overlay()
        self.scene.invalidate()
        self.viewport().update()

    def add_image_item(self, item: ImageAnnotationItem) -> None:
        """Añade una imagen anotada persistente al lienzo."""
        self.ensure_graphics_item_opacity(item)
        if item.scene() is not self.scene:
            self.scene.addItem(item)
        if item not in self.image_items:
            self.image_items.append(item)

    def add_energy_diagram_item(self, item: EnergyDiagramItem) -> None:
        """Añade un diagrama de energia persistente al lienzo."""
        self.ensure_graphics_item_opacity(item)
        if item.scene() is not self.scene:
            self.scene.addItem(item)
        if item not in self.energy_diagram_items:
            self.energy_diagram_items.append(item)

    def add_semantic_diagram_item(self, item: CompositeDiagramItem) -> None:
        """Añade un diagrama semántico compuesto al lienzo."""
        self.ensure_graphics_item_opacity(item)
        if item.scene() is not self.scene:
            self.scene.addItem(item)
        if item not in self.semantic_diagram_items:
            self.semantic_diagram_items.append(item)

    def add_orbital_item(self, item: OrbitalAnnotationItem) -> None:
        """Añade un orbital vectorial persistente al lienzo."""
        self.ensure_graphics_item_opacity(item)
        if item.scene() is not self.scene:
            self.scene.addItem(item)
        if item not in self.orbital_items:
            self.orbital_items.append(item)

    def add_wavy_anchor_item(self, item: WavyAnchorItem) -> None:
        """Añade wavy anchor item."""
        self.ensure_graphics_item_opacity(item)
        if item.scene() is not self.scene:
            self.scene.addItem(item)
        self._wavy_anchors.add(item)
        self._position_wavy_anchor(item)

    def add_plate_item(self, item: TLCPlateItem | GelElectrophoresisItem) -> None:
        """Añade una placa de análisis al lienzo."""
        if item.scene() is not self.scene:
            self.scene.addItem(item)
        if item not in self.plate_items:
            self.plate_items.append(item)

    def remove_image_item(self, item: ImageAnnotationItem) -> None:
        """Elimina una imagen anotada del lienzo."""
        if item in self.image_items:
            self.image_items.remove(item)
        if item.scene() is self.scene:
            self.scene.removeItem(item)

    def remove_energy_diagram_item(self, item: EnergyDiagramItem) -> None:
        """Elimina un diagrama de energia persistente."""
        if item in self.energy_diagram_items:
            self.energy_diagram_items.remove(item)
        if item.scene() is self.scene:
            self.scene.removeItem(item)

    def remove_semantic_diagram_item(self, item: CompositeDiagramItem) -> None:
        """Elimina un diagrama semántico compuesto."""
        if item in self.semantic_diagram_items:
            self.semantic_diagram_items.remove(item)
        if item.scene() is self.scene:
            self.scene.removeItem(item)

    def remove_orbital_item(self, item: OrbitalAnnotationItem) -> None:
        """Elimina un orbital persistente del lienzo."""
        if item in self.orbital_items:
            self.orbital_items.remove(item)
        if item.scene() is self.scene:
            self.scene.removeItem(item)

    def remove_wavy_anchor_item(self, item: WavyAnchorItem) -> None:
        """Elimina wavy anchor item."""
        if item.scene() is self.scene:
            self.scene.removeItem(item)
        if item in self._wavy_anchors:
            self._wavy_anchors.discard(item)

    def remove_plate_item(self, item: TLCPlateItem | GelElectrophoresisItem) -> None:
        """Elimina una placa de analisis del lienzo.

        Removing the plate removes all its children (lanes, spots, bands, labels)
        from the scene automatically. We don't need to do it manually.
        """
        if item in self.plate_items:
            self.plate_items.remove(item)
        try:
            if item.scene() is self.scene:
                self.scene.removeItem(item)
        except RuntimeError:
            pass

    def readd_image_item(self, item: ImageAnnotationItem) -> None:
        """Reintroduce una imagen anotada en el lienzo."""
        self.ensure_graphics_item_opacity(item)
        if item.scene() is not self.scene:
            self.scene.addItem(item)
        if item not in self.image_items:
            self.image_items.append(item)

    def readd_energy_diagram_item(self, item: EnergyDiagramItem) -> None:
        """Reintroduce un diagrama de energia persistente en el lienzo."""
        self.ensure_graphics_item_opacity(item)
        if item.scene() is not self.scene:
            self.scene.addItem(item)
        if item not in self.energy_diagram_items:
            self.energy_diagram_items.append(item)

    def readd_semantic_diagram_item(self, item: CompositeDiagramItem) -> None:
        """Reintroduce un diagrama semántico compuesto en el lienzo."""
        self.ensure_graphics_item_opacity(item)
        if item.scene() is not self.scene:
            self.scene.addItem(item)
        if item not in self.semantic_diagram_items:
            self.semantic_diagram_items.append(item)

    def readd_orbital_item(self, item: OrbitalAnnotationItem) -> None:
        """Reintroduce un orbital persistente en el lienzo."""
        self.ensure_graphics_item_opacity(item)
        if item.scene() is not self.scene:
            self.scene.addItem(item)
        if item not in self.orbital_items:
            self.orbital_items.append(item)

    def readd_plate_item(self, item) -> None:
        """Reintroduce una placa en el lienzo."""
        self.ensure_graphics_item_opacity(item)
        if item.scene() is not self.scene:
            self.scene.addItem(item)
        if item not in self.plate_items:
            self.plate_items.append(item)

    def readd_wavy_anchor_item(self, item: WavyAnchorItem) -> None:
        """Método auxiliar para readd wavy anchor item.

        Args:
            item: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        self.ensure_graphics_item_opacity(item)
        if item.scene() is not self.scene:
            self.scene.addItem(item)
        self._wavy_anchors.add(item)
        self._position_wavy_anchor(item)

    def _render_scene_bounds(self, selected_only: bool = False) -> Optional[QRectF]:
        """Método auxiliar para  render scene bounds.

        Args:
            selected_only: Descripción del parámetro.

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
            if rect is None:
                rect = candidate
            else:
                rect = rect.united(candidate)

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
            atom = self.model.atoms.get(atom_id)
            draws_coordination_sphere = bool(
                atom is not None
                and getattr(atom, "is_coordination_center", False)
                and not getattr(atom, "sphere_transparent", False)
            )
            if (
                item.pen().style() != Qt.PenStyle.NoPen
                or item.brush().style() != Qt.BrushStyle.NoBrush
                or draws_coordination_sphere
            ):
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

        if selected_only:
            for atom_id in self._selected_atom_ids_for_transform():
                extend_atom_bounds(atom_id)
            for bond_id in self.state.selected_bonds:
                item = self.bond_items.get(bond_id)
                if item is None or item.scene() is not self.scene:
                    continue
                extend(item.sceneBoundingRect())
            for item in self.scene.selectedItems():
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
        else:
            for atom_id in self.atom_items.keys():
                extend_atom_bounds(atom_id)
            for item in self.bond_items.values():
                if item.scene() is not self.scene:
                    continue
                if item.isVisible():
                    extend(item.sceneBoundingRect())
            for item in self.aromatic_circles:
                if item.scene() is not self.scene:
                    continue
                if item.isVisible():
                    extend(item.sceneBoundingRect())
            for item in self.arrow_items:
                if item.scene() is not self.scene:
                    continue
                if item.isVisible():
                    extend(item.sceneBoundingRect())
            for item in self.bracket_items:
                if item.scene() is not self.scene:
                    continue
                if item.isVisible():
                    extend(item.sceneBoundingRect())
            for item in self.orbital_items:
                if item.scene() is not self.scene:
                    continue
                if item.isVisible():
                    extend(item.sceneBoundingRect())
            for item in self.semantic_diagram_items:
                if item.scene() is not self.scene:
                    continue
                if item.isVisible():
                    extend(item.sceneBoundingRect())
            for item in self.image_items:
                if item.scene() is not self.scene:
                    continue
                if item.isVisible():
                    extend(item.sceneBoundingRect())
            for item in self._wavy_anchors:
                if item.scene() is not self.scene:
                    continue
                if item.isVisible():
                    extend(item.sceneBoundingRect())
            for item in self.plate_items:
                if item.scene() is not self.scene:
                    continue
                if item.isVisible():
                    extend(item.sceneBoundingRect())
                    # Spots/bands are standalone scene items — include them.
                    for lane in item.lane_items:
                        if hasattr(lane, "rf_labels"):
                            for spot, _ in lane.rf_labels:
                                if spot.scene() is self.scene and spot.isVisible():
                                    extend(spot.sceneBoundingRect())
                        elif hasattr(lane, "bands"):
                            for band, _ in lane.bands:
                                if band.scene() is self.scene and band.isVisible():
                                    extend(band.sceneBoundingRect())
            for item in self.scene.items():
                if (
                    isinstance(item, TextAnnotationItem)
                    and item.isVisible()
                    and not bool(item.data(NUMBERING_TEXT_ROLE))
                ):
                    extend(item.sceneBoundingRect())
            if self.state.numbering_enabled and self.state.numbering_include_in_export:
                for item in self._numbering_overlay_items:
                    if item is None or item.scene() is not self.scene:
                        continue
                    if item.isVisible():
                        extend(item.sceneBoundingRect())

        if rect is None:
            return None
        pad = max(1.0, self.drawing_style.stroke_px)
        return rect.adjusted(-pad, -pad, pad, pad)

    def _hidden_render_items(self) -> list:
        """Método auxiliar para  hidden render items.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        items = list(self._paper_scene_items())
        overlay_attrs = [
            "_hover_atom_indicator",
            "_hover_bond_indicator",
            "_optimize_zone",
            "_preview_bond_item",
            "_preview_ring_item",
            "_preview_chain_item",
            "_preview_chain_label",
            "_preview_arrow_item",
            "_grid_minor_item",
            "_grid_major_item",
        ]
        for attr in overlay_attrs:
            items.append(getattr(self, attr, None))
        items.append(getattr(self, "_bracket_preview", None))
        items.append(getattr(self, "_select_preview_path", None))
        items.append(getattr(self, "_select_preview_rect", None))
        if not self.state.numbering_include_in_export:
            items.extend(self._numbering_overlay_items)
        return [
            item
            for item in items
            if item is not None and item.scene() is self.scene
        ]

    def _with_hidden_render_items(self, render_fn):
        """Método auxiliar para  with hidden render items.

        Args:
            render_fn: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        hidden = []
        disabled_effects: list[tuple[QGraphicsItem, object, bool]] = []
        for item in self._hidden_render_items():
            try:
                effect = item.graphicsEffect()
            except Exception:
                effect = None
            if effect is not None:
                try:
                    was_enabled = bool(effect.isEnabled())
                    effect.setEnabled(False)
                    disabled_effects.append((item, effect, was_enabled))
                except Exception:
                    pass
            if item.isVisible():
                hidden.append(item)
                item.setVisible(False)

        selected_items = list(self.scene.selectedItems())
        saved_anchor = self.bond_anchor_id
        hovered_atom_id = self.hovered_atom_id
        hovered_bond_id = self.hovered_bond_id

        if hovered_atom_id in self.atom_items:
            self.atom_items[hovered_atom_id].set_hover(False)
        self.hovered_atom_id = None
        self.hovered_bond_id = None

        if selected_items or saved_anchor is not None:
            self.scene.blockSignals(True)
            self.scene.clearSelection()
            self.bond_anchor_id = None
            self.scene.blockSignals(False)
            self._sync_selection_from_scene()

        try:
            return render_fn()
        finally:
            if selected_items or saved_anchor is not None:
                self.scene.blockSignals(True)
                for item in selected_items:
                    item.setSelected(True)
                self.bond_anchor_id = saved_anchor
                self.scene.blockSignals(False)
                self._sync_selection_from_scene()

            self.hovered_atom_id = hovered_atom_id
            self.hovered_bond_id = hovered_bond_id
            if hovered_atom_id in self.atom_items:
                self.atom_items[hovered_atom_id].set_hover(True)

            for item in hidden:
                item.setVisible(True)
            for _item, effect, was_enabled in disabled_effects:
                try:
                    effect.setEnabled(was_enabled)
                except Exception:
                    pass

    def _with_hidden_unselected(self, render_fn):
        """Método auxiliar para  with hidden unselected.

        Args:
            render_fn: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        hidden = []
        selected = set(self.scene.selectedItems())
        for atom_id in self.state.selected_atoms:
            item = self.atom_items.get(atom_id)
            if item is not None:
                selected.add(item)
        for bond_id in self.state.selected_bonds:
            item = self.bond_items.get(bond_id)
            if item is not None:
                selected.add(item)
        def should_keep_visible(item):
            # 1. Direct selection or parent of selection (for recursive children visibility)
            if item in selected:
                return True
            
            # 2. Standalone logical associates (spots/bands)
            if isinstance(item, (TLCSpotItem, GelBandItem)):
                if item.lane_ref is not None:
                    plat = item.lane_ref.parentItem()
                    if plat in selected:
                        return True
            
            # 3. Descendants of any item in selected (Recursive check)
            p = item.parentItem()
            while p is not None:
                if p in selected:
                    return True
                p = p.parentItem()
            
            return False

        for item in self.scene.items():
            if item.zValue() >= 45: # Do not hide selection handles/overlay itself
                continue
            if should_keep_visible(item):
                continue
            if item.isVisible():
                hidden.append(item)
                item.setVisible(False)
        try:
            return render_fn()
        finally:
            for item in hidden:
                item.setVisible(True)

    def _render_scene_image(
        self,
        scale: float = 1.0,
        selected_only: bool = False,
        background: Optional[QColor] = None,
    ) -> Optional[QImage]:
        """Método auxiliar para  render scene image.

        Args:
            scale: Descripción del parámetro.
            selected_only: Descripción del parámetro.
            background: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        rect = self._render_scene_bounds(selected_only=selected_only)
        if rect is None:
            return None
        scale = max(1.0, float(scale))
        width = max(1, math.ceil(rect.width() * scale))
        height = max(1, math.ceil(rect.height() * scale))

        def render():
            """Método auxiliar para render.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
            image = QImage(width, height, QImage.Format.Format_ARGB32)
            image.fill(Qt.GlobalColor.transparent)
            painter = QPainter(image)
            painter.setRenderHint(QPainter.RenderHint.Antialiasing)
            painter.setRenderHint(QPainter.RenderHint.TextAntialiasing)
            self.scene.render(painter, QRectF(0, 0, width, height), rect)
            painter.end()
            trimmed = self._trim_transparent_image(image)
            if trimmed is None:
                return None
            if background is not None:
                flattened = QImage(trimmed.width(), trimmed.height(), QImage.Format.Format_RGB32)
                flattened.fill(background)
                painter = QPainter(flattened)
                painter.setRenderHint(QPainter.RenderHint.Antialiasing)
                painter.setRenderHint(QPainter.RenderHint.TextAntialiasing)
                painter.drawImage(0, 0, trimmed)
                painter.end()
                self._apply_image_dpi(flattened, scale)
                return flattened
            self._apply_image_dpi(trimmed, scale)
            return trimmed

        if selected_only:
            return self._with_hidden_unselected(
                lambda: self._with_hidden_render_items(render)
            )
        return self._with_hidden_render_items(render)

    def _apply_image_dpi(self, image: QImage, scale: float) -> None:
        """Método auxiliar para  apply image dpi.

        Args:
            image: Descripción del parámetro.
            scale: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        base_dpi_x = self.logicalDpiX() or 96.0
        base_dpi_y = self.logicalDpiY() or 96.0
        dpi_x = base_dpi_x * scale
        dpi_y = base_dpi_y * scale
        dpm_x = max(1, round(dpi_x * 1000.0 / 25.4))
        dpm_y = max(1, round(dpi_y * 1000.0 / 25.4))
        image.setDotsPerMeterX(dpm_x)
        image.setDotsPerMeterY(dpm_y)

    def _trim_transparent_image(self, image: QImage) -> Optional[QImage]:
        """Método auxiliar para  trim transparent image.

        Args:
            image: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        width = image.width()
        height = image.height()
        if width <= 0 or height <= 0:
            return None
        left = width
        right = -1
        top = height
        bottom = -1
        for y in range(height):
            for x in range(width):
                if (image.pixel(x, y) >> 24) & 0xFF:
                    if x < left:
                        left = x
                    if x > right:
                        right = x
                    if y < top:
                        top = y
                    if y > bottom:
                        bottom = y
        if right < left or bottom < top:
            return None
        pad = 1
        left = max(0, left - pad)
        top = max(0, top - pad)
        right = min(width - 1, right + pad)
        bottom = min(height - 1, bottom + pad)
        return image.copy(QRect(left, top, right - left + 1, bottom - top + 1))

    def _render_scene_svg(self, selected_only: bool = False) -> Optional[bytes]:
        """Método auxiliar para  render scene svg.

        Args:
            selected_only: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if QSvgGenerator is None:
            return None
        rect = self._render_scene_bounds(selected_only=selected_only)
        if rect is None:
            return None
        width = max(1, math.ceil(rect.width()))
        height = max(1, math.ceil(rect.height()))

        def render():
            """Método auxiliar para render.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la escena.
            """
            buffer = QBuffer()
            buffer.open(QBuffer.OpenModeFlag.WriteOnly)
            generator = QSvgGenerator()
            generator.setOutputDevice(buffer)
            generator.setSize(QSize(width, height))
            generator.setViewBox(QRect(0, 0, width, height))
            painter = QPainter(generator)
            self.scene.render(painter, QRectF(0, 0, width, height), rect)
            painter.end()
            return bytes(buffer.data())

        if selected_only:
            return self._with_hidden_unselected(
                lambda: self._with_hidden_render_items(render)
            )
        return self._with_hidden_render_items(render)

    @staticmethod
    def _encode_qimage_png(image: QImage) -> bytes:
        """Serializa un QImage a PNG para persistencia interna."""
        buffer = QBuffer()
        buffer.open(QBuffer.OpenModeFlag.WriteOnly)
        image.save(buffer, "PNG")
        return bytes(buffer.data())

    def _image_entries_from_clipboard(self, mime: QMimeData) -> list[dict]:
        """Extrae imágenes pegables desde URLs locales o datos crudos del portapapeles."""
        entries: list[dict] = []

        if mime.hasUrls():
            for url in mime.urls():
                if not url.isLocalFile():
                    continue
                path = url.toLocalFile()
                ext = os.path.splitext(path)[1].lower()
                mime_type = SUPPORTED_IMAGE_FILE_MIME_TYPES.get(ext)
                if mime_type is None:
                    continue
                try:
                    with open(path, "rb") as fh:
                        data = fh.read()
                except Exception:
                    continue
                item = ImageAnnotationItem(data, mime_type, source_name=os.path.basename(path))
                if not item.is_valid_image():
                    continue
                entries.append(
                    {
                        "data": data,
                        "mime_type": mime_type,
                        "source_name": os.path.basename(path),
                        "natural_width": item.boundingRect().width(),
                        "natural_height": item.boundingRect().height(),
                    }
                )
            if entries:
                return entries

        if mime.hasFormat("image/svg+xml"):
            data = bytes(mime.data("image/svg+xml"))
            item = ImageAnnotationItem(data, "image/svg+xml")
            if item.is_valid_image():
                entries.append(
                    {
                        "data": data,
                        "mime_type": "image/svg+xml",
                        "source_name": "",
                        "natural_width": item.boundingRect().width(),
                        "natural_height": item.boundingRect().height(),
                    }
                )
                return entries

        raw_image = None
        if mime.hasImage():
            image = mime.imageData()
            if isinstance(image, QImage):
                raw_image = image
            elif isinstance(image, QPixmap):
                raw_image = image.toImage()
        elif mime.hasFormat("image/png"):
            png_data = bytes(mime.data("image/png"))
            if png_data:
                item = ImageAnnotationItem(png_data, "image/png")
                if item.is_valid_image():
                    entries.append(
                        {
                            "data": png_data,
                            "mime_type": "image/png",
                            "source_name": "",
                            "natural_width": item.boundingRect().width(),
                            "natural_height": item.boundingRect().height(),
                        }
                    )
                    return entries
        if raw_image is not None and not raw_image.isNull():
            png_data = self._encode_qimage_png(raw_image)
            item = ImageAnnotationItem(png_data, "image/png")
            if item.is_valid_image():
                entries.append(
                    {
                        "data": png_data,
                        "mime_type": "image/png",
                        "source_name": "",
                        "natural_width": item.boundingRect().width(),
                        "natural_height": item.boundingRect().height(),
                    }
                )
        return entries

    def _insert_images_from_clipboard(self, mime: QMimeData) -> bool:
        """Inserta imágenes del portapapeles y las deja seleccionadas."""
        entries = self._image_entries_from_clipboard(mime)
        if not entries:
            return False

        target = self._last_scene_pos
        if target is None:
            target = self.mapToScene(self.viewport().rect().center())

        self.undo_stack.beginMacro("Paste image")
        created: list[ImageAnnotationItem] = []
        try:
            for index, entry in enumerate(entries):
                natural_width = max(1.0, float(entry["natural_width"]))
                natural_height = max(1.0, float(entry["natural_height"]))
                max_width = max(32.0, (self.paper_width - 2.0 * PAPER_MARGIN) * PASTE_IMAGE_MAX_PAPER_FRACTION)
                max_height = max(32.0, (self.paper_height - 2.0 * PAPER_MARGIN) * PASTE_IMAGE_MAX_PAPER_FRACTION)
                fit_scale = min(max_width / natural_width, max_height / natural_height, 1.0)
                width = natural_width * fit_scale
                height = natural_height * fit_scale
                item = ImageAnnotationItem(
                    entry["data"],
                    entry["mime_type"],
                    width=width,
                    height=height,
                    source_name=entry.get("source_name"),
                )
                if not item.is_valid_image():
                    continue
                offset = float(index) * PASTE_IMAGE_OFFSET_PX
                item.set_display_rect(
                    QRectF(
                        target.x() - width / 2.0 + offset,
                        target.y() - height / 2.0 + offset,
                        width,
                        height,
                    )
                )
                self.undo_stack.push(AddImageItemCommand(self, item))
                created.append(item)
        finally:
            self.undo_stack.endMacro()

        if created:
            if self.current_tool != "tool_select":
                self.set_current_tool("tool_select")
            self._select_inserted_items(items=created)
            self._refresh_scene_after_image_insert()
            return True
        return False

    def add_arrow_item(
        self,
        start: QPointF,
        end: QPointF,
        kind: str,
        curve_factor: float | None = None,
        stroke_px: float | None = None,
        opacity: Optional[float] | object = ITEM_OPACITY_UNSET,
    ) -> ArrowItem:
        """Añade arrow item.

        Args:
            start: Descripción del parámetro.
            end: Descripción del parámetro.
            kind: Descripción del parámetro.
            curve_factor: Curvatura normalizada para flechas curvas.
            stroke_px: Grosor personalizado de flecha.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = ArrowItem(
            start,
            end,
            kind=kind,
            curve_factor=curve_factor,
            stroke_px=stroke_px,
            style=self.drawing_style,
        )
        self.ensure_graphics_item_opacity(item, opacity)
        self.scene.addItem(item)
        self.arrow_items.append(item)
        return item

    def readd_arrow_item(
        self,
        item: ArrowItem,
        start: QPointF,
        end: QPointF,
        kind: str,
        curve_factor: float | None = None,
        opacity: Optional[float] | object = ITEM_OPACITY_UNSET,
    ) -> None:
        """Método auxiliar para readd arrow item.

        Args:
            item: Descripción del parámetro.
            start: Descripción del parámetro.
            end: Descripción del parámetro.
            kind: Descripción del parámetro.
            curve_factor: Curvatura normalizada para flechas curvas.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item.set_kind(kind)
        item.update_positions(start, end, curve_factor=curve_factor)
        self.ensure_graphics_item_opacity(item, opacity)
        if item.scene() is not self.scene:
            self.scene.addItem(item)
        if item not in self.arrow_items:
            self.arrow_items.append(item)

    def remove_arrow_item(self, item: ArrowItem) -> None:
        """Elimina arrow item.

        Args:
            item: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if item in self.arrow_items:
            self.arrow_items.remove(item)
        if item.scene() is self.scene:
            self.scene.removeItem(item)

    def add_bracket_item(
        self,
        rect: QRectF,
        kind: str,
        *,
        padding: float = 8.0,
        stroke_px: float | None = None,
        opacity: Optional[float] | object = ITEM_OPACITY_UNSET,
    ) -> BracketItem:
        """Añade bracket item.

        Args:
            rect: Descripción del parámetro.
            kind: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item = BracketItem(
            rect,
            kind=kind,
            padding=padding,
            stroke_px=stroke_px,
            style=self.drawing_style,
        )
        self.ensure_graphics_item_opacity(item, opacity)
        self.scene.addItem(item)
        self.bracket_items.append(item)
        return item

    def readd_bracket_item(
        self,
        item: BracketItem,
        rect: QRectF,
        kind: str,
        padding: Optional[float] = None,
        stroke_px: float | None = None,
        opacity: Optional[float] | object = ITEM_OPACITY_UNSET,
    ) -> None:
        """Método auxiliar para readd bracket item.

        Args:
            item: Descripción del parámetro.
            rect: Descripción del parámetro.
            kind: Descripción del parámetro.
            padding: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        item.set_rect(rect, padding=padding)
        item._kind = kind
        item.set_stroke_px(stroke_px)
        self.ensure_graphics_item_opacity(item, opacity)
        if item.scene() is not self.scene:
            self.scene.addItem(item)
        if item not in self.bracket_items:
            self.bracket_items.append(item)

    def remove_bracket_item(self, item: BracketItem) -> None:
        """Elimina bracket item.

        Args:
            item: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if item in self.bracket_items:
            self.bracket_items.remove(item)
        if item.scene() is self.scene:
            self.scene.removeItem(item)

    def refresh_aromatic_circles(self) -> None:
        """Add/remove aromatic circle items based on current state."""
        use_circles = self.state.use_aromatic_circles
        rings = self._detect_aromatic_rings(with_atoms=True) if use_circles else []
        ring_pairs_with_circle: set[frozenset[int]] = set()
        for ring in rings:
            ring_pairs_with_circle.update(ring.get("bond_pairs", set()))

        # Update aromatic bond items:
        # - Render as single-edge only when a visible aromatic circle exists
        #   for that ring bond.
        # - Fall back to normal bond rendering otherwise.
        for bond_id, item in self.bond_items.items():
            if item.is_aromatic:
                bond = self.model.bonds.get(bond_id)
                render_as_circle = (
                    use_circles
                    and bond is not None
                    and frozenset({bond.a1_id, bond.a2_id}) in ring_pairs_with_circle
                )
                item.set_render_aromatic_as_circle(render_as_circle)
                self.update_bond_item(bond_id)

        # Remove existing circles
        for circle in self.aromatic_circles:
            self.scene.removeItem(circle)
        self.aromatic_circles.clear()
        
        if not use_circles:
            return
        
        # Draw circles for complete aromatic rings.
        for ring in rings:
            atom_ids = sorted(ring.get("atoms", set()))
            self._draw_aromatic_circle(atom_ids)

    def _trackball_affine_context(
        self,
    ) -> Optional[tuple[float, float, float, float, float, float]]:
        """Devuelve el contexto afín del trackball si la referencia sigue válida."""
        atom_ids = self._rotation_3d_ref_atom_ids
        if not atom_ids or not self._rotation_3d_ref_positions:
            return None
        if any(atom_id not in self.model.atoms for atom_id in atom_ids):
            return None
        if any(atom_id not in self._rotation_3d_ref_positions for atom_id in atom_ids):
            return None

        if not self._is_rotating_3d:
            projected = self._project_trackball_reference(
                atom_ids,
                self._rotation_3d_pitch_deg,
                self._rotation_3d_yaw_deg,
            )
            max_delta = 0.0
            for atom_id in atom_ids:
                px, py = projected[atom_id]
                atom = self.model.get_atom(atom_id)
                max_delta = max(max_delta, math.hypot(px - atom.x, py - atom.y))
            if max_delta > TRACKBALL_REFERENCE_MATCH_TOLERANCE_PX:
                return None

        count = float(len(atom_ids))
        center_x = sum(self._rotation_3d_ref_positions[atom_id][0] for atom_id in atom_ids) / count
        center_y = sum(self._rotation_3d_ref_positions[atom_id][1] for atom_id in atom_ids) / count
        pitch_rad = math.radians(self._rotation_3d_pitch_deg)
        yaw_rad = math.radians(self._rotation_3d_yaw_deg)
        return (
            center_x,
            center_y,
            math.cos(pitch_rad),
            math.sin(pitch_rad),
            math.cos(yaw_rad),
            math.sin(yaw_rad),
        )

    def _aromatic_circle_geometry_from_trackball_reference(
        self,
        atom_ids: Iterable[int],
        affine_ctx: tuple[float, float, float, float, float, float],
    ) -> Optional[tuple[float, float, float, float, float]]:
        """Calcula la elipse aromática aplicando la misma afín pseudo-3D."""
        ref_points: list[tuple[float, float]] = []
        for atom_id in atom_ids:
            point = self._rotation_3d_ref_positions.get(int(atom_id))
            if point is None:
                return None
            ref_points.append(point)
        if len(ref_points) < 3:
            return None

        ring_cx = sum(x for x, _ in ref_points) / len(ref_points)
        ring_cy = sum(y for _, y in ref_points) / len(ref_points)
        avg_dist = (
            sum(math.hypot(x - ring_cx, y - ring_cy) for x, y in ref_points)
            / len(ref_points)
        )
        base_radius = max(2.0, avg_dist * 0.65)

        center_x, center_y, cos_x, sin_x, cos_y, sin_y = affine_ctx
        x0 = ring_cx - center_x
        y0 = ring_cy - center_y
        y_new = y0 * cos_x
        z_new = y0 * sin_x
        proj_cx = x0 * cos_y + z_new * sin_y + center_x
        proj_cy = y_new + center_y

        # Imagen afín del círculo base de radio r en el plano XY.
        v1x = base_radius * cos_y
        v1y = 0.0
        v2x = base_radius * sin_x * sin_y
        v2y = base_radius * cos_x

        c11 = v1x * v1x + v2x * v2x
        c22 = v1y * v1y + v2y * v2y
        c12 = v1x * v1y + v2x * v2y
        trace = c11 + c22
        disc = max(0.0, (c11 - c22) * (c11 - c22) + 4.0 * c12 * c12)
        root = math.sqrt(disc)
        eig_1 = max(0.0, 0.5 * (trace + root))
        eig_2 = max(0.0, 0.5 * (trace - root))

        radius_x = max(2.0, math.sqrt(eig_1))
        radius_y = max(2.0, math.sqrt(eig_2))
        angle_deg = 0.5 * math.degrees(math.atan2(2.0 * c12, c11 - c22))
        return proj_cx, proj_cy, radius_x, radius_y, angle_deg

    def _aromatic_circle_geometry_for_atom_ids(
        self,
        atom_ids: Iterable[int],
        trackball_affine: Optional[tuple[float, float, float, float, float, float]] = None,
    ) -> Optional[tuple[float, float, float, float, float]]:
        """Calcula la geometría del marcador aromático."""
        atom_id_list = [int(atom_id) for atom_id in atom_ids]
        if len(atom_id_list) < 3:
            return None

        ring_atoms = []
        for atom_id in atom_id_list:
            atom = self.model.atoms.get(atom_id)
            if atom is not None:
                ring_atoms.append(atom)
        if len(ring_atoms) < 3:
            return None

        cx = sum(atom.x for atom in ring_atoms) / len(ring_atoms)
        cy = sum(atom.y for atom in ring_atoms) / len(ring_atoms)
        avg_dist = (
            sum(math.hypot(atom.x - cx, atom.y - cy) for atom in ring_atoms)
            / len(ring_atoms)
        )
        circular_radius = max(2.0, avg_dist * 0.65)

        centered = [(atom.x - cx, atom.y - cy) for atom in ring_atoms]
        c11 = sum(dx * dx for dx, _ in centered) / len(centered)
        c22 = sum(dy * dy for _, dy in centered) / len(centered)
        c12 = sum(dx * dy for dx, dy in centered) / len(centered)
        trace = c11 + c22
        disc = max(0.0, (c11 - c22) * (c11 - c22) + 4.0 * c12 * c12)
        root = math.sqrt(disc)
        eig_1 = max(0.0, 0.5 * (trace + root))
        eig_2 = max(0.0, 0.5 * (trace - root))

        max_eig = max(eig_1, eig_2)
        if max_eig <= 1e-9:
            return cx, cy, circular_radius, circular_radius, 0.0

        # Evita jitter angular en anillos prácticamente isotrópicos.
        if abs(eig_1 - eig_2) <= max(1e-6, max_eig * 1e-4):
            return cx, cy, circular_radius, circular_radius, 0.0

        # Para puntos distribuidos sobre una elipse, var ~= r^2 / 2.
        scale = 0.65 * math.sqrt(2.0)
        radius_x = max(2.0, math.sqrt(eig_1) * scale)
        radius_y = max(2.0, math.sqrt(eig_2) * scale)
        angle_deg = 0.5 * math.degrees(math.atan2(2.0 * c12, c11 - c22))
        return cx, cy, radius_x, radius_y, angle_deg

    def _draw_aromatic_circle(
        self,
        atom_ids: Iterable[int],
        trackball_affine: Optional[tuple[float, float, float, float, float, float]] = None,
    ) -> None:
        """Crea y añade a escena el círculo interno de un anillo aromático."""
        atom_ids_tuple = tuple(sorted(int(atom_id) for atom_id in atom_ids))
        if not atom_ids_tuple:
            return
        geometry = self._aromatic_circle_geometry_for_atom_ids(
            atom_ids_tuple,
            trackball_affine=trackball_affine,
        )
        if geometry is None:
            return
        cx, cy, radius_x, radius_y, angle_deg = geometry
        circle = AromaticCircleItem(QRectF(-radius_x, -radius_y, 2.0 * radius_x, 2.0 * radius_y))
        circle.set_geometry(cx, cy, radius_x, radius_y, angle_deg)
        circle.setData(AROMATIC_CIRCLE_ATOMS_ROLE, atom_ids_tuple)
        circle.setOpacity(self._aromatic_circle_effective_opacity(atom_ids_tuple))
        self.scene.addItem(circle)
        self.aromatic_circles.append(circle)

    def _update_aromatic_circles_for_atoms(self, atom_ids: set[int]) -> None:
        """Actualiza posición/tamaño de círculos aromáticos afectados por átomos movidos."""
        if not self.state.use_aromatic_circles or not self.aromatic_circles:
            return
        if not atom_ids:
            return

        changed_atoms = set(atom_ids)
        force_refresh = False
        for circle in list(self.aromatic_circles):
            try:
                ring_atom_ids = circle.data(AROMATIC_CIRCLE_ATOMS_ROLE)
            except RuntimeError:
                force_refresh = True
                break
            if not ring_atom_ids:
                force_refresh = True
                break
            ring_set = {int(atom_id) for atom_id in ring_atom_ids}
            if changed_atoms.isdisjoint(ring_set):
                continue
            geometry = self._aromatic_circle_geometry_for_atom_ids(
                ring_set,
            )
            if geometry is None:
                force_refresh = True
                break
            cx, cy, radius_x, radius_y, angle_deg = geometry
            circle.set_geometry(cx, cy, radius_x, radius_y, angle_deg)
            circle.setOpacity(self._aromatic_circle_effective_opacity(ring_set))

        if force_refresh:
            self.refresh_aromatic_circles()

    def _detect_aromatic_rings(self, with_atoms: bool = False) -> list:
        """
        Detect complete aromatic rings and return their centers.
        Returns list of dicts with center and radius (and optionally atoms).
        """
        view = MolView(self.model)
        rings = find_rings_simple(view)
        if not rings:
            return []

        bond_by_pair = {
            frozenset({bond.a1_id, bond.a2_id}): bond
            for bond in self.model.bonds.values()
        }
        aromatic_pairs = {
            pair
            for pair, bond in bond_by_pair.items()
            if bond.is_aromatic
        }
        if not aromatic_pairs:
            return []
        
        # Calculate centers and radii for each ring
        result = []
        for ring in rings:
            ring_pairs = ring_bonds(view, ring)
            if not ring_pairs:
                continue
            # Modo estricto: todos los enlaces del ciclo son aromáticos.
            missing_pairs = [pair for pair in ring_pairs if pair not in aromatic_pairs]
            if missing_pairs:
                # Tolerancia visual: permitir hasta dos enlaces no aromáticos
                # cuando fueron estilizados para proyección 2D
                # (cuña/hash/bold). Esto evita perder el círculo aromático por
                # un ajuste de estilo local en anillos aromáticos.
                if len(missing_pairs) > 2:
                    continue
                missing_bonds = [bond_by_pair.get(pair) for pair in missing_pairs]
                if any(bond is None for bond in missing_bonds):
                    continue
                if any(
                    bond.style not in {BondStyle.WEDGE, BondStyle.HASHED, BondStyle.BOLD}
                    for bond in missing_bonds
                ):
                    continue
                aromatic_count = len(ring_pairs) - len(missing_pairs)
                if aromatic_count < max(1, len(ring_pairs) - 2):
                    continue
            if not ring_pairs:
                continue

            atom_ids = list(ring)
            atoms = [self.model.get_atom(aid) for aid in atom_ids if aid in self.model.atoms]
            if len(atoms) < 3:
                continue
            cx = sum(a.x for a in atoms) / len(atoms)
            cy = sum(a.y for a in atoms) / len(atoms)
            avg_dist = sum(((a.x - cx)**2 + (a.y - cy)**2)**0.5 for a in atoms) / len(atoms)
            radius = max(2.0, avg_dist * 0.65)
            entry = {
                "center": QPointF(cx, cy),
                "radius": radius,
            }
            if with_atoms:
                entry["atoms"] = set(atom_ids)
                entry["bond_pairs"] = set(ring_pairs)
            result.append(entry)
        
        return result

    def _aromatic_ring_center_for_bond(self, bond: Bond) -> Optional[QPointF]:
        """Método auxiliar para  aromatic ring center for bond.

        Args:
            bond: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if not bond.is_aromatic:
            return None
        adjacency: dict[int, list[tuple[int, int]]] = {}
        for b in self.model.bonds.values():
            if not b.is_aromatic:
                continue
            adjacency.setdefault(b.a1_id, []).append((b.a2_id, b.id))
            adjacency.setdefault(b.a2_id, []).append((b.a1_id, b.id))
        cycle = self._find_cycle_for_bond(bond, adjacency)
        if not cycle:
            return None
        return self._center_for_atoms(cycle["atom_ids"])

    def _refresh_aromatic_ring_contexts(self) -> None:
        """Método auxiliar para  refresh aromatic ring contexts.

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
        if not adjacency:
            for bond in self.model.bonds.values():
                if bond.is_aromatic and bond.ring_id is None:
                    item = self.bond_items.get(bond.id)
                    if item is not None:
                        item.set_ring_context(None)
                        atom1 = self.model.get_atom(bond.a1_id)
                        atom2 = self.model.get_atom(bond.a2_id)
                        if atom1 and atom2:
                            item.update_positions(atom1, atom2)
            return
        for bond in self.model.bonds.values():
            if not bond.is_aromatic or bond.ring_id is not None:
                continue
            item = self.bond_items.get(bond.id)
            if item is None:
                continue
            cycle = self._find_cycle_for_bond(bond, adjacency)
            center = self._center_for_atoms(cycle["atom_ids"]) if cycle else None
            item.set_ring_context(center)
            atom1 = self.model.get_atom(bond.a1_id)
            atom2 = self.model.get_atom(bond.a2_id)
            if atom1 and atom2:
                item.update_positions(atom1, atom2)

    def _build_aromatic_adjacency(self) -> dict[int, list[int]]:
        """Método auxiliar para  build aromatic adjacency.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        adjacency: dict[int, list[int]] = {}
        for bond in self.model.bonds.values():
            if not bond.is_aromatic:
                continue
            adjacency.setdefault(bond.a1_id, []).append(bond.a2_id)
            adjacency.setdefault(bond.a2_id, []).append(bond.a1_id)
        return adjacency

    def _center_for_atoms(self, atom_ids: list[int]) -> Optional[QPointF]:
        """Método auxiliar para  center for atoms.

        Args:
            atom_ids: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la escena.
        """
        if not atom_ids:
            return None
        xs = []
        ys = []
        for aid in atom_ids:
            atom = self.model.get_atom(aid)
            if atom is None:
                continue
            xs.append(atom.x)
            ys.append(atom.y)
        if not xs:
            return None
        return QPointF(sum(xs) / len(xs), sum(ys) / len(ys))

    def apply_drawing_style(self, style: DrawingStyle) -> None:
        """Apply a new drawing style and refresh items."""
        self.drawing_style = style
        self.state.bond_length = style.bond_length_px
        for item in self.atom_items.values():
            item.set_style(style)
        self.update_bond_items_for_atoms(set(self.model.atoms.keys()))
        for arrow in self.arrow_items:
            arrow.set_style(style)
        for bracket in self.bracket_items:
            bracket.set_style(style)
        self.refresh_atom_visibility()
