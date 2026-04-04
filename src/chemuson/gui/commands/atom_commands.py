from __future__ import annotations

from ._shared import *

class AddAtomCommand(QUndoCommand):
    """Comando para añadir un átomo al modelo y a la vista."""

    def __init__(
        self,
        model: MolGraph,
        view,
        element: str,
        x: float,
        y: float,
        is_explicit: Optional[bool] = None,
        charge: int | None = None,
        isotope: Optional[int] = None,
        radical_electrons: int = 0,
        oxidation_state: Optional[int] = None,
        explicit_h: Optional[int] = None,
        group_h_cap: Optional[int] = None,
        mapping: Optional[int] = None,
        is_query: bool = False,
        anchor_override: Optional[str] = None,
        auto_hydrogens: bool = True,
        expected_bonds: int = 0,
        no_implicit: bool = False,
        label_scale: Optional[float] = None,
        is_coordination_center: bool = False,
        sphere_radius: Optional[float] = None,
        sphere_color: Optional[str] = None,
        sphere_filled: bool = True,
        sphere_transparent: bool = False,
        opacity: Optional[float] = None,
    ) -> None:
        """Inicializa el comando de adición de átomo.

        Args:
            model: Modelo molecular a modificar.
            view: Vista/lienzo asociado.
            element: Símbolo del elemento.
            x: Coordenada X.
            y: Coordenada Y.
            is_explicit: Si el símbolo debe mostrarse explícitamente.
            charge: Carga formal.
            isotope: Isótopo (número másico).
            radical_electrons: Número de electrones desapareados.
            oxidation_state: Estado de oxidación si aplica.
            explicit_h: Hidrógenos explícitos.
            mapping: Índice de mapeo.
            is_query: Marca de átomo de consulta.
            anchor_override: Ancla visual alternativa.
            auto_hydrogens: Si se auto-generan hidrógenos.
            expected_bonds: Número de enlaces esperados (para H implícitos).
            no_implicit: Si se desactivan hidrógenos implícitos en el átomo.
            label_scale: Escala local de etiqueta o `None` para heredar.
            is_coordination_center: Si el átomo se dibuja como esfera de coordinación.
            sphere_radius: Radio visual de la esfera de coordinación.
            sphere_color: Color base de la esfera (hex) o `None`.
            sphere_filled: Si la esfera se dibuja con relleno.
            sphere_transparent: Si la esfera se dibuja transparente (sin borde/fondo).
            opacity: Opacidad local del átomo o `None` para heredar del canvas.
        """
        super().__init__("Add atom")
        resolved_label = _resolve_atom_label_spec(view, element)
        self._model = model
        self._view = view
        self._element = str(resolved_label.get("element") or element)
        self._x = x
        self._y = y
        self._is_explicit = is_explicit
        self._charge = charge
        self._isotope = isotope
        self._radical_electrons = int(radical_electrons or 0)
        self._oxidation_state = (
            int(oxidation_state)
            if oxidation_state is not None
            else None
        )
        resolved_explicit_h = resolved_label.get("explicit_h")
        self._explicit_h = explicit_h if explicit_h is not None else resolved_explicit_h
        resolved_group_h_cap = (
            group_h_cap
            if group_h_cap is not None
            else resolved_label.get("group_h_cap")
        )
        self._group_h_cap = (
            None if resolved_group_h_cap is None else int(resolved_group_h_cap)
        )
        self._mapping = mapping
        self._is_query = is_query
        self._anchor_override = anchor_override
        self._auto_hydrogens = auto_hydrogens
        self._expected_bonds = expected_bonds
        self._no_implicit = bool(no_implicit or resolved_label.get("no_implicit", False))
        self._label_scale = label_scale
        self._is_coordination_center = bool(is_coordination_center)
        self._sphere_radius = sphere_radius
        self._sphere_color = sphere_color
        self._sphere_filled = bool(sphere_filled)
        self._sphere_transparent = bool(sphere_transparent)
        self._opacity = opacity
        self._atom_id: Optional[int] = None
        self._hydrogen_specs: list[tuple[int, float, float, int]] = []

    def redo(self) -> None:
        """Ejecuta la adición del átomo (redo)."""
        is_explicit = self._is_explicit
        if is_explicit is None:
            is_explicit = _default_is_explicit(self._element)
        charge = 0 if self._charge is None else self._charge
        if self._atom_id is None:
            atom = self._model.add_atom(
                self._element,
                self._x,
                self._y,
                is_explicit=is_explicit,
                charge=charge,
                isotope=self._isotope,
                radical_electrons=self._radical_electrons,
                oxidation_state=self._oxidation_state,
                explicit_h=self._explicit_h,
                group_h_cap=self._group_h_cap,
                mapping=self._mapping,
                is_query=self._is_query,
                no_implicit=self._no_implicit,
                label_scale=self._label_scale,
                is_coordination_center=self._is_coordination_center,
                sphere_radius=self._sphere_radius,
                sphere_color=self._sphere_color,
                sphere_filled=self._sphere_filled,
                sphere_transparent=self._sphere_transparent,
                opacity=self._opacity,
            )
            self._atom_id = atom.id
        else:
            atom = self._model.add_atom(
                self._element,
                self._x,
                self._y,
                atom_id=self._atom_id,
                is_explicit=is_explicit,
                charge=charge,
                isotope=self._isotope,
                radical_electrons=self._radical_electrons,
                oxidation_state=self._oxidation_state,
                explicit_h=self._explicit_h,
                group_h_cap=self._group_h_cap,
                mapping=self._mapping,
                is_query=self._is_query,
                no_implicit=self._no_implicit,
                label_scale=self._label_scale,
                is_coordination_center=self._is_coordination_center,
                sphere_radius=self._sphere_radius,
                sphere_color=self._sphere_color,
                sphere_filled=self._sphere_filled,
                sphere_transparent=self._sphere_transparent,
                opacity=self._opacity,
            )
        self._view.add_atom_item(atom)
        if self._anchor_override:
            self._view.set_anchor_override(atom.id, self._anchor_override)
            self._view._refresh_atom_label(atom.id)
        _validate_view_structure(self._view)

    def undo(self) -> None:
        """Deshace la adición del átomo y sus hidrógenos."""
        if self._hydrogen_specs:
            _remove_hydrogen_specs(self._model, self._view, self._hydrogen_specs)
        atom, removed_bonds = self._model.remove_atom(self._atom_id)
        for bond in removed_bonds:
            self._view.remove_bond_item(bond.id)
        self._view.remove_atom_item(atom.id)
        _validate_view_structure(self._view)

    @property
    def atom_id(self) -> Optional[int]:
        """Devuelve el ID del átomo creado (si existe)."""
        return self._atom_id

class ChangeAtomCommand(QUndoCommand):
    """Comando para cambiar el elemento de un átomo existente."""

    def __init__(
        self,
        model: MolGraph,
        view,
        atom_id: int,
        new_element: str,
        anchor_override=_ANCHOR_UNSET,
    ) -> None:
        """Inicializa el comando de cambio de átomo.

        Args:
            model: Modelo molecular.
            view: Vista/lienzo asociado.
            atom_id: ID del átomo a modificar.
            new_element: Nuevo símbolo de elemento.
            anchor_override: Ancla visual opcional.
        """
        super().__init__("Change atom")
        resolved_label = _resolve_atom_label_spec(view, new_element)
        self._model = model
        self._view = view
        self._atom_id = atom_id
        self._old_element = model.get_atom(atom_id).element
        self._old_is_explicit = model.get_atom(atom_id).is_explicit
        self._old_explicit_h = model.get_atom(atom_id).explicit_h
        self._old_group_h_cap = getattr(model.get_atom(atom_id), "group_h_cap", None)
        self._old_no_implicit = bool(getattr(model.get_atom(atom_id), "no_implicit", False))
        self._new_element = str(resolved_label.get("element") or new_element)
        self._new_explicit_h = resolved_label.get("explicit_h")
        resolved_group_h_cap = resolved_label.get("group_h_cap")
        self._new_group_h_cap = (
            None if resolved_group_h_cap is None else int(resolved_group_h_cap)
        )
        self._new_no_implicit = bool(resolved_label.get("no_implicit", False))
        self._new_is_explicit = _default_is_explicit(self._new_element)
        self._anchor_override = anchor_override
        if anchor_override is _ANCHOR_UNSET:
            self._old_anchor_override = _ANCHOR_UNSET
        else:
            self._old_anchor_override = view.get_anchor_override(atom_id)
        self._added_hydrogen_specs: list[tuple[int, float, float, int]] = []
        self._removed_hydrogens = []
        self._removed_hydrogen_bonds = []
        self._removed_hydrogen_specs: list[tuple[int, float, float, int]] = []

    def redo(self) -> None:
        """Aplica el cambio de elemento y actualiza la vista."""
        # Check if we need to remove hydrogens due to valence
        if not self._removed_hydrogen_specs:
            self._check_and_remove_hydrogens()
        else:
            # Re-apply removal if we already calculated it
            if self._removed_hydrogen_specs:
                _remove_hydrogen_specs(self._model, self._view, self._removed_hydrogen_specs)
        if self._anchor_override is not _ANCHOR_UNSET:
            self._view.set_anchor_override(self._atom_id, self._anchor_override)
        self._model.update_atom_element(
            self._atom_id,
            self._new_element,
            is_explicit=self._new_is_explicit,
        )
        atom = self._model.get_atom(self._atom_id)
        atom.explicit_h = self._new_explicit_h
        atom.group_h_cap = self._new_group_h_cap
        atom.no_implicit = self._new_no_implicit
        self._view.update_atom_item_element(
            self._atom_id,
            self._new_element,
            is_explicit=self._new_is_explicit,
        )
        _validate_view_structure(self._view)

    def undo(self) -> None:
        """Revierte el cambio de elemento y restaura hidrógenos."""
        if self._anchor_override is not _ANCHOR_UNSET:
            self._view.set_anchor_override(self._atom_id, self._old_anchor_override)
        self._model.update_atom_element(
            self._atom_id,
            self._old_element,
            is_explicit=self._old_is_explicit,
        )
        atom = self._model.get_atom(self._atom_id)
        atom.explicit_h = self._old_explicit_h
        atom.group_h_cap = self._old_group_h_cap
        atom.no_implicit = self._old_no_implicit
        self._view.update_atom_item_element(
            self._atom_id,
            self._old_element,
            is_explicit=self._old_is_explicit,
        )
        
        # Restore hydrogens
        if self._removed_hydrogen_specs:
            _readd_hydrogen_specs(self._model, self._view, self._atom_id, self._removed_hydrogen_specs)
        _validate_view_structure(self._view)

    def _check_and_remove_hydrogens(self) -> None:
        """Elimina hidrógenos explícitos si exceden la valencia permitida."""
        from chemuson.core.model import VALENCE_MAP
        
        # 1. Get current bonds (excluding explicit H)
        bonds = self._model.bonds.values()
        non_h_bonds = 0
        attached_hydrogens = []
        
        for bond in bonds:
            if bond.a1_id == self._atom_id:
                other_id = bond.a2_id
            elif bond.a2_id == self._atom_id:
                other_id = bond.a1_id
            else:
                continue
                
            other = self._model.get_atom(other_id)
            if other.element == "H" and other.is_explicit:
                # Check if it's a terminal H
                if _atom_degree(self._model, other_id) == 1:
                    attached_hydrogens.append((other, bond))
                    continue
            
            non_h_bonds += bond.order

        # 2. Get max valence for new element
        max_valence = VALENCE_MAP.get(self._new_element, 0)
        
        # 3. Calculate allowed H
        allowed_h = max(0, max_valence - non_h_bonds)
        
        # 4. If we have more H than allowed, mark for removal
        if len(attached_hydrogens) > allowed_h:
            # Remove excess hydrogens
            # We remove all if 0 allowed, or just the excess
            # Usually users expect all attached H to be recalculated or kept if they fit?
            # If I change C to O in a ring (2 bonds). Valence O=2. allowed_h = 0.
            # If there was an H attached (making it CH-), degree was 3.
            # Now degree is 2 (bonds) + 1 (H) = 3 > 2. So H must go.
            
            # Sort hydrogens by id or position to be deterministic?
            # Just take the excess
            excess_count = len(attached_hydrogens) - allowed_h
            to_remove = attached_hydrogens[:excess_count]
            
            specs = []
            for h_atom, h_bond in to_remove:
                specs.append((h_atom.id, h_atom.x, h_atom.y, h_bond.id))
            
            self._removed_hydrogen_specs = specs
            _remove_hydrogen_specs(self._model, self._view, self._removed_hydrogen_specs)

class ChangeChargeCommand(QUndoCommand):
    """Comando para modificar la carga formal de un átomo."""

    def __init__(self, model: MolGraph, view, atom_id: int, new_charge: int) -> None:
        """Inicializa el comando de cambio de carga.

        Args:
            model: Modelo molecular.
            view: Vista/lienzo asociado.
            atom_id: ID del átomo a modificar.
            new_charge: Nueva carga formal.
        """
        super().__init__("Change charge")
        self._model = model
        self._view = view
        self._atom_id = atom_id
        self._old_charge = model.get_atom(atom_id).charge
        self._new_charge = new_charge

    def redo(self) -> None:
        """Aplica el cambio de carga."""
        self._model.update_atom_charge(self._atom_id, self._new_charge)
        self._view.update_atom_item_charge(self._atom_id, self._new_charge)
        _validate_view_structure(self._view)

    def undo(self) -> None:
        """Revierte el cambio de carga."""
        self._model.update_atom_charge(self._atom_id, self._old_charge)
        self._view.update_atom_item_charge(self._atom_id, self._old_charge)
        _validate_view_structure(self._view)

class ChangeNoImplicitCommand(QUndoCommand):
    """Comando para activar/desactivar H implícitos en un átomo."""

    def __init__(self, model: MolGraph, view, atom_id: int, enabled: bool) -> None:
        super().__init__("Change no implicit H")
        self._model = model
        self._view = view
        self._atom_id = atom_id
        atom = model.get_atom(atom_id)
        self._old_enabled = bool(getattr(atom, "no_implicit", False))
        self._new_enabled = bool(enabled)

    def redo(self) -> None:
        self._model.update_atom_no_implicit(self._atom_id, self._new_enabled)
        self._view.refresh_atom_labels([self._atom_id])
        _validate_view_structure(self._view)

    def undo(self) -> None:
        self._model.update_atom_no_implicit(self._atom_id, self._old_enabled)
        self._view.refresh_atom_labels([self._atom_id])
        _validate_view_structure(self._view)

class ChangeAtomLabelScaleCommand(QUndoCommand):
    """Comando para ajustar la escala local de etiqueta de un átomo."""

    def __init__(
        self,
        model: MolGraph,
        view,
        atom_id: int,
        new_scale: Optional[float],
        old_scale: Optional[float] | object = _LABEL_SCALE_UNSET,
    ) -> None:
        super().__init__("Change atom label scale")
        self._model = model
        self._view = view
        self._atom_id = atom_id
        atom = model.get_atom(atom_id)
        self._old_scale = (
            getattr(atom, "label_scale", None)
            if old_scale is _LABEL_SCALE_UNSET
            else old_scale
        )
        self._new_scale = new_scale

    def _apply(self, value: Optional[float]) -> None:
        if self._atom_id not in self._model.atoms:
            return
        self._model.update_atom_label_scale(self._atom_id, value)
        self._view.refresh_label_fonts([self._atom_id])
        if hasattr(self._view, "_update_selection_overlay"):
            self._view._update_selection_overlay()

    def redo(self) -> None:
        self._apply(self._new_scale)

    def undo(self) -> None:
        self._apply(self._old_scale)

class SetCoordinationCenterCommand(QUndoCommand):
    """Comando para activar/desactivar visualización de centro de coordinación."""

    def __init__(self, model: MolGraph, view, atom_id: int, enabled: bool) -> None:
        """Inicializa el comando de cambio de centro de coordinación."""
        super().__init__("Set coordination center")
        self._model = model
        self._view = view
        self._atom_id = atom_id
        atom = model.get_atom(atom_id)
        self._old_enabled = bool(getattr(atom, "is_coordination_center", False))
        self._new_enabled = bool(enabled)
        self._old_sphere_radius = getattr(atom, "sphere_radius", None)
        self._old_sphere_color = getattr(atom, "sphere_color", None)
        self._old_sphere_filled = bool(getattr(atom, "sphere_filled", True))
        self._old_sphere_transparent = bool(getattr(atom, "sphere_transparent", False))

    def _apply(self, enabled: bool) -> None:
        """Aplica el estado solicitado y sincroniza el item visual."""
        if self._atom_id not in self._model.atoms:
            return
        atom = self._model.get_atom(self._atom_id)
        atom.is_coordination_center = bool(enabled)
        if enabled and getattr(atom, "sphere_radius", None) is None:
            atom.sphere_radius = 16.0
        if enabled and getattr(atom, "sphere_color", None) is None:
            atom.sphere_color = "#D9DDE3"
        if enabled and not hasattr(atom, "sphere_filled"):
            atom.sphere_filled = True
        if enabled and not hasattr(atom, "sphere_transparent"):
            atom.sphere_transparent = False
        item = self._view.atom_items.get(self._atom_id)
        if item is not None:
            item.set_coordination_center(bool(enabled))
        else:
            self._view.update_atom_item(self._atom_id, atom.x, atom.y)
        self._view.update_bond_items_for_atoms({self._atom_id})

    def redo(self) -> None:
        """Aplica el cambio."""
        self._apply(self._new_enabled)

    def undo(self) -> None:
        """Revierte el cambio."""
        self._apply(self._old_enabled)
        atom = self._model.atoms.get(self._atom_id)
        if atom is None:
            return
        atom.sphere_radius = self._old_sphere_radius
        atom.sphere_color = self._old_sphere_color
        atom.sphere_filled = self._old_sphere_filled
        atom.sphere_transparent = self._old_sphere_transparent
        item = self._view.atom_items.get(self._atom_id)
        if item is not None and hasattr(item, "refresh_coordination_visual"):
            item.refresh_coordination_visual()
        self._view.update_bond_items_for_atoms({self._atom_id})

class ChangeCoordinationSphereStyleCommand(QUndoCommand):
    """Comando para cambiar estilo visual de una esfera de coordinación."""

    def __init__(
        self,
        model: MolGraph,
        view,
        atom_id: int,
        new_radius: Optional[float] | object = _SPHERE_STYLE_UNSET,
        new_color: Optional[str] | object = _SPHERE_STYLE_UNSET,
        new_filled: bool | object = _SPHERE_STYLE_UNSET,
        new_transparent: bool | object = _SPHERE_STYLE_UNSET,
        old_radius: Optional[float] | object = _SPHERE_STYLE_UNSET,
    ) -> None:
        """Inicializa el cambio de estilo de esfera."""
        super().__init__("Change coordination sphere style")
        self._model = model
        self._view = view
        self._atom_id = atom_id
        atom = model.get_atom(atom_id)
        self._old_radius = (
            getattr(atom, "sphere_radius", None)
            if old_radius is _SPHERE_STYLE_UNSET
            else old_radius
        )
        self._old_color = getattr(atom, "sphere_color", None)
        self._old_filled = bool(getattr(atom, "sphere_filled", True))
        self._old_transparent = bool(getattr(atom, "sphere_transparent", False))
        self._new_radius = (
            self._old_radius if new_radius is _SPHERE_STYLE_UNSET else new_radius
        )
        self._new_color = (
            self._old_color if new_color is _SPHERE_STYLE_UNSET else new_color
        )
        self._new_filled = (
            self._old_filled if new_filled is _SPHERE_STYLE_UNSET else bool(new_filled)
        )
        self._new_transparent = (
            self._old_transparent
            if new_transparent is _SPHERE_STYLE_UNSET
            else bool(new_transparent)
        )

    def _apply(
        self,
        radius: Optional[float] | object,
        color: Optional[str] | object,
        filled: bool | object,
        transparent: bool | object,
    ) -> None:
        """Aplica estado visual sobre átomo e item."""
        atom = self._model.atoms.get(self._atom_id)
        if atom is None:
            return
        if radius is not _SPHERE_STYLE_UNSET:
            atom.sphere_radius = None if radius is None else max(4.0, float(radius))
        if color is not _SPHERE_STYLE_UNSET:
            atom.sphere_color = None if color is None else str(color)
        if filled is not _SPHERE_STYLE_UNSET:
            atom.sphere_filled = bool(filled)
        if transparent is not _SPHERE_STYLE_UNSET:
            atom.sphere_transparent = bool(transparent)
        item = self._view.atom_items.get(self._atom_id)
        if item is not None and hasattr(item, "refresh_coordination_visual"):
            item.refresh_coordination_visual()
        else:
            self._view.update_atom_item(self._atom_id, atom.x, atom.y)
        self._view.update_bond_items_for_atoms({self._atom_id})

    def redo(self) -> None:
        """Aplica nuevo estilo."""
        self._apply(
            self._new_radius,
            self._new_color,
            self._new_filled,
            self._new_transparent,
        )

    def undo(self) -> None:
        """Restaura estilo anterior."""
        self._apply(
            self._old_radius,
            self._old_color,
            self._old_filled,
            self._old_transparent,
        )
