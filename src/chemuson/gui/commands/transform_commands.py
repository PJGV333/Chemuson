from __future__ import annotations

from ._shared import *

class ChangeCanvasOpacityCommand(QUndoCommand):
    """Comando para ajustar opacidad local o global de elementos del canvas."""

    def __init__(
        self,
        model: MolGraph,
        view,
        *,
        atom_values: Optional[Dict[int, Optional[float]]] = None,
        bond_values: Optional[Dict[int, Optional[float]]] = None,
        item_values: Optional[Dict[object, Optional[float]]] = None,
        canvas_opacity: Optional[float] = None,
        text: str = "Change opacity",
    ) -> None:
        super().__init__(text)
        self._model = model
        self._view = view
        self._atom_after = dict(atom_values or {})
        self._bond_after = dict(bond_values or {})
        self._item_after = dict(item_values or {})
        self._canvas_after = canvas_opacity
        self._applies_canvas_default = canvas_opacity is not None
        self._canvas_before = float(getattr(getattr(view, "state", None), "canvas_opacity", 1.0))
        self._atom_before = {
            atom_id: getattr(model.get_atom(atom_id), "opacity", None)
            for atom_id in self._atom_after
            if atom_id in model.atoms
        }
        self._bond_before = {
            bond_id: getattr(model.get_bond(bond_id), "opacity", None)
            for bond_id in self._bond_after
            if bond_id in model.bonds
        }
        raw_item_opacity = getattr(view, "item_raw_opacity", None)
        self._item_before: Dict[object, Optional[float]] = {}
        if callable(raw_item_opacity):
            for item in self._item_after:
                try:
                    self._item_before[item] = raw_item_opacity(item)
                except Exception:
                    continue

    def _apply(
        self,
        atom_values: Dict[int, Optional[float]],
        bond_values: Dict[int, Optional[float]],
        item_values: Dict[object, Optional[float]],
        canvas_opacity: Optional[float],
    ) -> None:
        if self._applies_canvas_default and canvas_opacity is not None and hasattr(self._view, "set_canvas_opacity_default"):
            self._view.set_canvas_opacity_default(float(canvas_opacity))

        refresh_aromatic = False
        for atom_id, opacity in atom_values.items():
            if atom_id not in self._model.atoms:
                continue
            self._model.update_atom_opacity(atom_id, opacity)
            if hasattr(self._view, "update_atom_item_opacity"):
                self._view.update_atom_item_opacity(atom_id)

        for bond_id, opacity in bond_values.items():
            if bond_id not in self._model.bonds:
                continue
            self._model.update_bond(bond_id, opacity=opacity)
            if hasattr(self._view, "update_bond_item"):
                self._view.update_bond_item(bond_id)
            refresh_aromatic = True

        if hasattr(self._view, "set_graphics_item_opacity"):
            for item, opacity in item_values.items():
                try:
                    self._view.set_graphics_item_opacity(item, opacity)
                except RuntimeError:
                    continue

        if self._applies_canvas_default and hasattr(self._view, "refresh_numbering_opacity"):
            self._view.refresh_numbering_opacity()
        if (refresh_aromatic or self._applies_canvas_default) and hasattr(
            self._view,
            "refresh_aromatic_circle_opacities",
        ):
            self._view.refresh_aromatic_circle_opacities()
        if hasattr(self._view, "_update_selection_overlay"):
            self._view._update_selection_overlay()

    def redo(self) -> None:
        self._apply(
            self._atom_after,
            self._bond_after,
            self._item_after,
            self._canvas_after,
        )

    def undo(self) -> None:
        canvas_before = self._canvas_before if self._applies_canvas_default else None
        self._apply(
            self._atom_before,
            self._bond_before,
            self._item_before,
            canvas_before,
        )

class MoveAtomsCommand(QUndoCommand):
    """Comando para mover átomos y actualizar enlaces."""

    def __init__(
        self,
        model: MolGraph,
        view,
        before: Dict[int, Tuple[float, float]],
        after: Dict[int, Tuple[float, float]],
        skip_first_redo: bool = False,
    ) -> None:
        """Inicializa el comando de movimiento de átomos."""
        super().__init__("Move atoms")
        self._model = model
        self._view = view
        self._before = before
        self._after = after
        self._skip_first_redo = skip_first_redo
        self._first_redo = True

    def redo(self) -> None:
        """Aplica el movimiento (redo) considerando la primera ejecución."""
        if self._skip_first_redo and self._first_redo:
            self._first_redo = False
            return
        self._apply_positions(self._after)

    def undo(self) -> None:
        """Revierte las posiciones anteriores."""
        self._apply_positions(self._before)

    def _apply_positions(self, positions: Dict[int, Tuple[float, float]]) -> None:
        """Aplica un conjunto de posiciones a átomos y refresca enlaces."""
        for atom_id, (x, y) in positions.items():
            self._model.update_atom_position(atom_id, x, y)
            self._view.update_atom_item(atom_id, x, y)
        self._view.update_bond_items_for_atoms(set(positions.keys()))
        if hasattr(self._view, "recompute_numbering"):
            self._view.recompute_numbering()
        if hasattr(self._view, "_update_selection_overlay"):
            self._view._update_selection_overlay()

class MoveArrowItemsCommand(QUndoCommand):
    """Comando para mover elementos de flecha."""

    def __init__(self, view, before: dict, after: dict, text: str = "Move arrows"):
        """Inicializa el comando de movimiento de flechas."""
        super().__init__(text)
        self._view = view
        self._before = before  # {item: (start, end)}
        self._after = after    # {item: (start, end)}

    def redo(self) -> None:
        """Aplica nuevas posiciones de flechas."""
        for item, (start, end) in self._after.items():
            item.update_positions(start, end)
        self._view._update_selection_overlay()

    def undo(self) -> None:
        """Revierte posiciones de flechas."""
        for item, (start, end) in self._before.items():
            item.update_positions(start, end)
        self._view._update_selection_overlay()

class MoveBracketItemsCommand(QUndoCommand):
    """Comando para mover elementos de corchetes."""

    def __init__(self, view, before: dict, after: dict):
        """Inicializa el comando de movimiento de corchetes."""
        super().__init__("Move brackets")
        self._view = view
        self._before = before  # {item: QRectF}
        self._after = after    # {item: QRectF}

    def redo(self) -> None:
        """Aplica nuevos rectángulos de corchetes."""
        for item, rect in self._after.items():
            item.set_rect(rect)
        self._view._update_selection_overlay()

    def undo(self) -> None:
        """Revierte los rectángulos de corchetes."""
        for item, rect in self._before.items():
            item.set_rect(rect)
        self._view._update_selection_overlay()

class TransformImageItemsCommand(QUndoCommand):
    """Comando para mover/escalar imágenes anotadas."""

    def __init__(self, view, before: dict, after: dict, text: str = "Transform images") -> None:
        super().__init__(text)
        self._view = view
        self._before = before
        self._after = after

    @staticmethod
    def _apply_snapshot(item, snapshot) -> None:
        pos, width, height, rotation = snapshot
        item.set_display_rect(QRectF(float(pos.x()), float(pos.y()), float(width), float(height)))
        item.setRotation(float(rotation))

    def redo(self) -> None:
        for item, snapshot in self._after.items():
            self._apply_snapshot(item, snapshot)
        self._view._update_selection_overlay()

    def undo(self) -> None:
        for item, snapshot in self._before.items():
            self._apply_snapshot(item, snapshot)
        self._view._update_selection_overlay()

class ScaleArrowItemsCommand(QUndoCommand):
    """Comando para escalar geometría y grosor de flechas."""

    def __init__(self, view, before: dict, after: dict) -> None:
        super().__init__("Scale arrows")
        self._view = view
        self._before = before
        self._after = after

    @staticmethod
    def _apply_snapshot(item, snapshot) -> None:
        start, end, stroke_px = snapshot
        item.update_positions(start, end)
        item.set_stroke_px(stroke_px)

    def redo(self) -> None:
        for item, snapshot in self._after.items():
            self._apply_snapshot(item, snapshot)
        self._view._update_selection_overlay()

    def undo(self) -> None:
        for item, snapshot in self._before.items():
            self._apply_snapshot(item, snapshot)
        self._view._update_selection_overlay()

class ScaleBracketItemsCommand(QUndoCommand):
    """Comando para escalar corchetes, padding y grosor."""

    def __init__(self, view, before: dict, after: dict) -> None:
        super().__init__("Scale brackets")
        self._view = view
        self._before = before
        self._after = after

    @staticmethod
    def _apply_snapshot(item, snapshot) -> None:
        rect, padding, stroke_px = snapshot
        item.set_rect(rect, padding=padding)
        item.set_stroke_px(stroke_px)

    def redo(self) -> None:
        for item, snapshot in self._after.items():
            self._apply_snapshot(item, snapshot)
        self._view._update_selection_overlay()

    def undo(self) -> None:
        for item, snapshot in self._before.items():
            self._apply_snapshot(item, snapshot)
        self._view._update_selection_overlay()

class DeleteSelectionCommand(QUndoCommand):
    """Comando para eliminar selección (átomos, enlaces y anotaciones)."""

    def __init__(
        self,
        model: MolGraph,
        view,
        atom_ids: Iterable[int],
        bond_ids: Iterable[int],
        arrow_items: Iterable = (),
        bracket_items: Iterable = (),
        text_items: Iterable = (),
        wavy_items: Iterable = (),
        image_items: Iterable = (),
        orbital_items: Iterable = (),
        energy_diagram_items: Iterable = (),
        semantic_diagram_items: Iterable = (),
        plate_items: Iterable = (),
    ) -> None:
        """Inicializa el comando de borrado de selección."""
        super().__init__("Delete selection")
        self._model = model
        self._view = view
        self._atom_ids = sorted(set(atom_ids))
        self._bond_ids = sorted(set(bond_ids))
        self._arrow_items = list(arrow_items)
        self._bracket_items = list(bracket_items)
        self._text_items = list(text_items)
        self._wavy_items = list(wavy_items)
        self._image_items = list(image_items)
        self._orbital_items = list(orbital_items)
        self._energy_diagram_items = list(energy_diagram_items)
        self._semantic_diagram_items = list(semantic_diagram_items)
        self._plate_items = list(plate_items)
        self._removed_atoms = []
        self._removed_bonds = []
        self._removed_arrows = []
        self._removed_brackets = []
        self._removed_texts = []
        self._removed_wavy = []
        self._removed_images = []
        self._removed_orbitals = []
        self._removed_energy_diagrams = []
        self._removed_semantic_diagrams = []
        self._removed_plates = []

    def _refresh_view_after_structure_change(self) -> None:
        """Sincroniza overlays derivados tras borrar/restaurar estructura."""
        if hasattr(self._view, "refresh_ring_centers"):
            self._view.refresh_ring_centers()
        if hasattr(self._view, "refresh_aromatic_circles"):
            self._view.refresh_aromatic_circles()
        if hasattr(self._view, "recompute_numbering"):
            self._view.recompute_numbering()
        _validate_view_structure(self._view)

    def redo(self) -> None:
        """Aplica el borrado y guarda copias para restaurar."""
        if not self._removed_atoms and not self._removed_bonds:
            captured_atom_ids: set[int] = set()
            for atom_id in self._atom_ids:
                if atom_id in self._model.atoms:
                    self._removed_atoms.append(replace(self._model.atoms[atom_id]))
                    captured_atom_ids.add(atom_id)
            for bond in self._model.bonds.values():
                if (
                    bond.id in self._bond_ids
                    or bond.a1_id in self._atom_ids
                    or bond.a2_id in self._atom_ids
                ):
                    self._removed_bonds.append(replace(bond))
            candidate_orphans: set[int] = set()
            for bond in self._removed_bonds:
                candidate_orphans.add(bond.a1_id)
                candidate_orphans.add(bond.a2_id)
            candidate_orphans.difference_update(captured_atom_ids)
            for atom_id in sorted(candidate_orphans):
                if atom_id not in self._model.atoms:
                    continue
                if not _is_hidden_carbon_placeholder(self._model, atom_id):
                    continue
                removed_incident = sum(
                    1
                    for bond in self._removed_bonds
                    if bond.a1_id == atom_id or bond.a2_id == atom_id
                )
                if _atom_degree(self._model, atom_id) - removed_incident > 0:
                    continue
                self._removed_atoms.append(replace(self._model.atoms[atom_id]))
                captured_atom_ids.add(atom_id)
            for item in self._arrow_items:
                self._removed_arrows.append(
                    (
                        item,
                        item.start_point(),
                        item.end_point(),
                        item.kind(),
                        item.curve_factor(),
                    )
                )
            for item in self._bracket_items:
                self._removed_brackets.append(
                    (item, item.base_rect(), item._padding, item._kind, item.stroke_px())
                )
            for item in self._text_items:
                self._removed_texts.append(item)
            for item in self._wavy_items:
                self._removed_wavy.append(item)
            for item in self._image_items:
                rect = item.display_rect()
                self._removed_images.append(
                    (
                        item,
                        QPointF(rect.topLeft()),
                        float(rect.width()),
                        float(rect.height()),
                        float(item.rotation()),
                    )
                )
            for item in self._orbital_items:
                self._removed_orbitals.append((item, item.anchor0(), item.anchor1()))
            for item in self._energy_diagram_items:
                rect = item.display_rect()
                self._removed_energy_diagrams.append(
                    (
                        item,
                        QPointF(rect.topLeft()),
                        float(rect.width()),
                        float(rect.height()),
                        float(item.rotation()),
                        item.config_payload(),
                    )
                )
            for item in self._semantic_diagram_items:
                self._removed_semantic_diagrams.append((item, item.to_json()))
            for item in self._plate_items:
                self._removed_plates.append((item, item.to_dict()))

        for bond in list(self._removed_bonds):
            if bond.id in self._model.bonds:
                self._model.remove_bond(bond.id)
                self._view.remove_bond_item(bond.id)
        for atom in list(self._removed_atoms):
            if atom.id in self._model.atoms:
                self._model.remove_atom(atom.id)
                self._view.remove_atom_item(atom.id)
        for item, _start, _end, _kind, _curve_factor in list(self._removed_arrows):
            self._view.remove_arrow_item(item)
        for item, _rect, _padding, _kind, _stroke_px in list(self._removed_brackets):
            self._view.remove_bracket_item(item)
        for item in self._removed_texts:
            self._view.remove_text_item(item)
        for item in self._removed_wavy:
            self._view.remove_wavy_anchor_item(item)
        for item, _pos, _width, _height, _rotation in self._removed_images:
            self._view.remove_image_item(item)
        for item, _anchor0, _anchor1 in self._removed_orbitals:
            self._view.remove_orbital_item(item)
        for item, _pos, _width, _height, _rotation, _payload in self._removed_energy_diagrams:
            self._view.remove_energy_diagram_item(item)
        for item, _payload in self._removed_semantic_diagrams:
            self._view.remove_semantic_diagram_item(item)
        for item, _data in self._removed_plates:
            self._view.remove_plate_item(item)
        self._refresh_view_after_structure_change()

    def undo(self) -> None:
        """Restaura los elementos eliminados."""
        for atom in self._removed_atoms:
            restored_atom = self._model.add_atom(
                atom.element,
                atom.x,
                atom.y,
                atom_id=atom.id,
                charge=atom.charge,
                isotope=atom.isotope,
                radical_electrons=int(getattr(atom, "radical_electrons", 0) or 0),
                oxidation_state=getattr(atom, "oxidation_state", None),
                explicit_h=atom.explicit_h,
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
                opacity=getattr(atom, "opacity", None),
            )
            self._view.add_atom_item(restored_atom)
        for bond in self._removed_bonds:
            self._model.add_bond(
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
                pi_offset_sign=getattr(bond, "pi_offset_sign", None),
                opacity=getattr(bond, "opacity", None),
            )
            self._view.add_bond_item(bond)
        for item, start, end, kind, curve_factor in self._removed_arrows:
            self._view.readd_arrow_item(item, start, end, kind, curve_factor=curve_factor)
        for item, rect, padding, kind, stroke_px in self._removed_brackets:
            self._view.readd_bracket_item(item, rect, kind, padding=padding, stroke_px=stroke_px)
        for item in self._removed_texts:
            self._view.readd_text_item(item)
        for item in self._removed_wavy:
            self._view.readd_wavy_anchor_item(item)
        for item, pos, width, height, rotation in self._removed_images:
            item.set_display_rect(QRectF(pos.x(), pos.y(), width, height))
            item.setRotation(rotation)
            self._view.readd_image_item(item)
        for item, anchor0, anchor1 in self._removed_orbitals:
            item.set_anchors(anchor0, anchor1)
            self._view.readd_orbital_item(item)
        for item, pos, width, height, rotation, payload in self._removed_energy_diagrams:
            item.apply_config_payload(payload)
            item.set_display_rect(QRectF(pos.x(), pos.y(), width, height))
            item.setRotation(rotation)
            self._view.readd_energy_diagram_item(item)
        for item, payload in self._removed_semantic_diagrams:
            item.set_display_rect(
                QRectF(
                    float(payload.get("x", 0.0)),
                    float(payload.get("y", 0.0)),
                    float(payload.get("width", 1.0)),
                    float(payload.get("height", 1.0)),
                )
            )
            item.setRotation(float(payload.get("rotation", 0.0)))
            self._view.readd_semantic_diagram_item(item)
        for item, data in self._removed_plates:
            item.load_dict(data, scene=self._view.scene)
            for lane in item.lane_items:
                if hasattr(lane, "rf_labels"):
                    for spot, _ in lane.rf_labels:
                        spot.setPos(spot.pos() + item.pos())
                        spot.lane_ref = lane
                elif hasattr(lane, "bands"):
                    for band, _ in lane.bands:
                        band.setPos(band.pos() + item.pos())
                        band.lane_ref = lane
            self._view.readd_plate_item(item)
        self._refresh_view_after_structure_change()
