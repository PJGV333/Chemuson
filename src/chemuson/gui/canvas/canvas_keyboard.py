from __future__ import annotations

from ._shared import (
    ArrowItem,
    AtomItem,
    Bond,
    BondItem,
    BracketItem,
    ChangeAtomCommand,
    ChangeDoubleBondOrientationCommand,
    CompositeDiagramItem,
    EnergyDiagramItem,
    GelBandItem,
    GelElectrophoresisItem,
    ImageAnnotationItem,
    MoveArrowItemsCommand,
    MoveAtomsCommand,
    MoveBracketItemsCommand,
    MoveTextItemsCommand,
    OrbitalAnnotationItem,
    QPointF,
    QWheelEvent,
    Qt,
    TLCPlateItem,
    TLCSpotItem,
    TextAnnotationItem,
    TransformEnergyDiagramItemsCommand,
    TransformImageItemsCommand,
    TransformOrbitalItemsCommand,
    WavyAnchorItem,
)

class CanvasKeyboardMixin:
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
