from __future__ import annotations

from typing import Optional, List, Tuple

from PyQt6.QtGui import QUndoCommand

from chemuson.core.model import BondStereo, BondStyle, MolGraph

from ._shared import (
    _BOND_LENGTH_UNSET,
    _BOND_STROKE_UNSET,
    _DONOR_UNSET,
    _FLEX_CURVE_UNSET,
    _atom_degree,
    _default_is_explicit,
    _remove_hydrogen_specs,
    _resolve_atom_label_spec,
    _validate_view_structure,
    replace,
)

class AddBondCommand(QUndoCommand):
    """Comando para añadir un enlace (y opcionalmente un átomo nuevo)."""

    def __init__(
        self,
        model: MolGraph,
        view,
        a1_id: int,
        a2_id: Optional[int],
        order: int = 1,
        style: BondStyle = BondStyle.PLAIN,
        stereo: BondStereo = BondStereo.NONE,
        is_aromatic: bool = False,
        display_order: Optional[int] = None,
        length_px: Optional[float] = None,
        ring_id: Optional[int] = None,
        stroke_px: Optional[float] = None,
        color: Optional[str] = None,
        donor_atom_id: Optional[int] = None,
        flex_curve_1: Optional[float] = None,
        flex_curve_2: Optional[float] = None,
        opacity: Optional[float] = None,
        new_atom_element: Optional[str] = None,
        new_atom_pos: Optional[Tuple[float, float]] = None,
    ) -> None:
        """Inicializa el comando de adición de enlace.

        Args:
            model: Modelo molecular.
            view: Vista/lienzo asociado.
            a1_id: ID del átomo inicial.
            a2_id: ID del átomo final (o `None` si se crea uno nuevo).
            order: Orden del enlace.
            style: Estilo visual del enlace.
            stereo: Estereoquímica del enlace.
            is_aromatic: Si el enlace se marca como aromático.
            display_order: Orden visual alternativo.
            length_px: Longitud de dibujo fija.
            ring_id: Identificador de anillo asociado.
            stroke_px: Grosor de línea.
            color: Color del enlace.
            donor_atom_id: ID del átomo donador (si es coordinativo).
            flex_curve_1: Curvatura normalizada de control 1 (estilo FLEX).
            flex_curve_2: Curvatura normalizada de control 2 (estilo FLEX).
            opacity: Opacidad local del enlace o `None` para heredar del canvas.
            new_atom_element: Elemento del átomo a crear si `a2_id` es `None`.
            new_atom_pos: Posición del átomo a crear.
        """
        super().__init__("Add bond")
        self._model = model
        self._view = view
        self._a1_id = a1_id
        self._a2_id = a2_id
        self._order = order
        self._style = style
        self._stereo = stereo
        self._is_aromatic = is_aromatic
        self._display_order = display_order
        self._length_px = length_px
        self._ring_id = ring_id
        self._stroke_px = stroke_px
        self._color = color
        self._donor_atom_id = donor_atom_id
        self._flex_curve_1 = flex_curve_1
        self._flex_curve_2 = flex_curve_2
        self._opacity = opacity
        self._bond_id: Optional[int] = None
        resolved_new_atom = (
            _resolve_atom_label_spec(view, new_atom_element)
            if new_atom_element is not None
            else None
        )
        self._new_atom_element = (
            str(resolved_new_atom.get("element") or new_atom_element)
            if resolved_new_atom is not None
            else new_atom_element
        )
        self._new_atom_explicit_h = (
            None if resolved_new_atom is None else resolved_new_atom.get("explicit_h")
        )
        resolved_group_h_cap = (
            None if resolved_new_atom is None else resolved_new_atom.get("group_h_cap")
        )
        self._new_atom_group_h_cap = (
            None if resolved_group_h_cap is None else int(resolved_group_h_cap)
        )
        self._new_atom_no_implicit = bool(
            resolved_new_atom is not None and resolved_new_atom.get("no_implicit", False)
        )
        self._new_atom_pos = new_atom_pos
        self._created_atom_id: Optional[int] = None
        self._hydrogen_specs: list[tuple[int, float, float, int]] = []
        self._demoted_explicit_atoms: Optional[list[int]] = None

    def _should_demote_explicit_carbon(self, atom_id: int) -> bool:
        """Decide si un carbono explícito debe ocultarse tras crear el enlace."""
        atom = self._model.get_atom(atom_id)
        if atom.element != "C" or not atom.is_explicit:
            return False
        if getattr(getattr(self._view, "state", None), "show_implicit_carbons", True):
            return False
        if (
            atom.charge != 0
            or atom.isotope is not None
            or int(getattr(atom, "radical_electrons", 0) or 0) > 0
            or atom.explicit_h is not None
        ):
            return False
        if atom.mapping is not None or atom.is_query:
            return False
        if _atom_degree(self._model, atom_id) <= 0:
            return False
        return True

    def redo(self) -> None:
        """Aplica la creación del enlace y del átomo opcional."""
        if self._a2_id is None:
            if self._created_atom_id is None:
                if self._new_atom_element is None or self._new_atom_pos is None:
                    return
                is_explicit = _default_is_explicit(self._new_atom_element)
                atom = self._model.add_atom(
                    self._new_atom_element,
                    self._new_atom_pos[0],
                    self._new_atom_pos[1],
                    is_explicit=is_explicit,
                    explicit_h=self._new_atom_explicit_h,
                    group_h_cap=self._new_atom_group_h_cap,
                    no_implicit=self._new_atom_no_implicit,
                )
                self._created_atom_id = atom.id
            else:
                is_explicit = _default_is_explicit(self._new_atom_element)
                atom = self._model.add_atom(
                    self._new_atom_element,
                    self._new_atom_pos[0],
                    self._new_atom_pos[1],
                    atom_id=self._created_atom_id,
                    is_explicit=is_explicit,
                    explicit_h=self._new_atom_explicit_h,
                    group_h_cap=self._new_atom_group_h_cap,
                    no_implicit=self._new_atom_no_implicit,
                )
            self._view.add_atom_item(atom)
            self._a2_id = self._created_atom_id

        if self._bond_id is None:
            bond = self._model.add_bond(
                self._a1_id,
                self._a2_id,
                self._order,
                style=self._style,
                stereo=self._stereo,
                is_aromatic=self._is_aromatic,
                display_order=self._display_order,
                ring_id=self._ring_id,
                length_px=self._length_px,
                stroke_px=self._stroke_px,
                color=self._color,
                donor_atom_id=self._donor_atom_id,
                flex_curve_1=self._flex_curve_1,
                flex_curve_2=self._flex_curve_2,
                opacity=self._opacity,
            )
            self._bond_id = bond.id
        else:
            bond = self._model.add_bond(
                self._a1_id,
                self._a2_id,
                self._order,
                bond_id=self._bond_id,
                style=self._style,
                stereo=self._stereo,
                is_aromatic=self._is_aromatic,
                display_order=self._display_order,
                ring_id=self._ring_id,
                length_px=self._length_px,
                stroke_px=self._stroke_px,
                color=self._color,
                donor_atom_id=self._donor_atom_id,
                flex_curve_1=self._flex_curve_1,
                flex_curve_2=self._flex_curve_2,
                opacity=self._opacity,
            )
        self._view.add_bond_item(bond)
        if self._demoted_explicit_atoms is None:
            demoted: list[int] = []
            for atom_id in (self._a1_id, self._a2_id):
                if atom_id is None:
                    continue
                if self._should_demote_explicit_carbon(atom_id):
                    demoted.append(atom_id)
            self._demoted_explicit_atoms = demoted
        for atom_id in self._demoted_explicit_atoms:
            if atom_id in self._model.atoms:
                atom = self._model.get_atom(atom_id)
                self._model.update_atom_element(atom_id, atom.element, is_explicit=False)
                self._view.update_atom_item_element(atom_id, atom.element, is_explicit=False)
        _validate_view_structure(self._view)

    def undo(self) -> None:
        """Deshace la creación del enlace y restaura estado previo."""
        if self._hydrogen_specs:
            _remove_hydrogen_specs(self._model, self._view, self._hydrogen_specs)
        bond = self._model.remove_bond(self._bond_id)
        self._view.remove_bond_item(bond.id)
        if self._created_atom_id is not None:
            atom, removed_bonds = self._model.remove_atom(self._created_atom_id)
            for removed in removed_bonds:
                self._view.remove_bond_item(removed.id)
            self._view.remove_atom_item(atom.id)
            self._a2_id = None
        if self._demoted_explicit_atoms:
            for atom_id in self._demoted_explicit_atoms:
                if atom_id in self._model.atoms:
                    atom = self._model.get_atom(atom_id)
                    self._model.update_atom_element(atom_id, atom.element, is_explicit=True)
                    self._view.update_atom_item_element(atom_id, atom.element, is_explicit=True)
        _validate_view_structure(self._view)

class ChangeBondCommand(QUndoCommand):
    """Comando para modificar propiedades de un enlace."""

    def __init__(
        self,
        model: MolGraph,
        view,
        bond_id: int,
        new_order: Optional[int] = None,
        new_style: Optional[BondStyle] = None,
        new_stereo: Optional[BondStereo] = None,
        new_is_aromatic: Optional[bool] = None,
        new_donor_atom_id: Optional[int] | object = _DONOR_UNSET,
        new_flex_curve_1: Optional[float] | object = _FLEX_CURVE_UNSET,
        new_flex_curve_2: Optional[float] | object = _FLEX_CURVE_UNSET,
    ) -> None:
        """Inicializa el comando de cambio de enlace.

        Args:
            model: Modelo molecular.
            view: Vista/lienzo asociado.
            bond_id: ID del enlace a modificar.
            new_order: Nuevo orden de enlace.
            new_style: Nuevo estilo de enlace.
            new_stereo: Nueva estereoquímica.
            new_is_aromatic: Nueva bandera de aromaticidad.
            new_donor_atom_id: Nuevo átomo donador (si aplica).
            new_flex_curve_1: Curvatura normalizada de control 1 (estilo FLEX).
            new_flex_curve_2: Curvatura normalizada de control 2 (estilo FLEX).
        """
        super().__init__("Change bond")
        self._model = model
        self._view = view
        self._bond_id = bond_id
        bond = model.get_bond(bond_id)
        self._old_order = bond.order
        self._old_style = bond.style
        self._old_stereo = bond.stereo
        self._old_is_aromatic = bond.is_aromatic
        self._old_donor_atom_id = getattr(bond, "donor_atom_id", None)
        self._old_flex_curve_1 = getattr(bond, "flex_curve_1", None)
        self._old_flex_curve_2 = getattr(bond, "flex_curve_2", None)
        self._new_order = new_order if new_order is not None else bond.order
        self._new_style = new_style if new_style is not None else bond.style
        self._new_stereo = new_stereo if new_stereo is not None else bond.stereo
        self._new_is_aromatic = (
            new_is_aromatic if new_is_aromatic is not None else bond.is_aromatic
        )
        self._new_donor_atom_id = (
            self._old_donor_atom_id
            if new_donor_atom_id is _DONOR_UNSET
            else new_donor_atom_id
        )
        self._new_flex_curve_1 = (
            self._old_flex_curve_1
            if new_flex_curve_1 is _FLEX_CURVE_UNSET
            else new_flex_curve_1
        )
        self._new_flex_curve_2 = (
            self._old_flex_curve_2
            if new_flex_curve_2 is _FLEX_CURVE_UNSET
            else new_flex_curve_2
        )

    def redo(self) -> None:
        """Aplica el cambio de propiedades del enlace."""
        self._model.update_bond(
            self._bond_id,
            order=self._new_order,
            style=self._new_style,
            stereo=self._new_stereo,
            is_aromatic=self._new_is_aromatic,
            donor_atom_id=self._new_donor_atom_id,
            flex_curve_1=self._new_flex_curve_1,
            flex_curve_2=self._new_flex_curve_2,
        )
        self._view.update_bond_item(self._bond_id)
        refresher = getattr(self._view, "refresh_aromatic_circles", None)
        if callable(refresher):
            refresher()
        _validate_view_structure(self._view)

    def undo(self) -> None:
        """Revierte el cambio de propiedades del enlace."""
        self._model.update_bond(
            self._bond_id,
            order=self._old_order,
            style=self._old_style,
            stereo=self._old_stereo,
            is_aromatic=self._old_is_aromatic,
            donor_atom_id=self._old_donor_atom_id,
            flex_curve_1=self._old_flex_curve_1,
            flex_curve_2=self._old_flex_curve_2,
        )
        self._view.update_bond_item(self._bond_id)
        refresher = getattr(self._view, "refresh_aromatic_circles", None)
        if callable(refresher):
            refresher()
        _validate_view_structure(self._view)

class ChangeBondLengthCommand(QUndoCommand):
    """Comando para cambiar la longitud visual de un enlace."""

    def __init__(
        self,
        model: MolGraph,
        view,
        bond_id: int,
        new_length: Optional[float],
        old_length: Optional[float] | object = _BOND_LENGTH_UNSET,
    ) -> None:
        """Inicializa el comando de cambio de longitud de enlace."""
        super().__init__("Change bond length")
        self._model = model
        self._view = view
        self._bond_id = bond_id
        bond = model.get_bond(bond_id)
        self._old_length = bond.length_px if old_length is _BOND_LENGTH_UNSET else old_length
        self._new_length = new_length

    def redo(self) -> None:
        """Aplica la nueva longitud."""
        self._model.update_bond_length(self._bond_id, self._new_length)
        self._view.update_bond_item(self._bond_id)

    def undo(self) -> None:
        """Revierte a la longitud anterior."""
        self._model.update_bond_length(self._bond_id, self._old_length)
        self._view.update_bond_item(self._bond_id)

class ChangeBondStrokeCommand(QUndoCommand):
    """Comando para cambiar el grosor de un enlace."""

    def __init__(
        self,
        model: MolGraph,
        view,
        bond_id: int,
        new_stroke_px: Optional[float],
        old_stroke_px: Optional[float] | object = _BOND_STROKE_UNSET,
    ) -> None:
        """Inicializa el comando de cambio de grosor."""
        super().__init__("Change bond thickness")
        self._model = model
        self._view = view
        self._bond_id = bond_id
        bond = model.get_bond(bond_id)
        self._old_stroke = (
            bond.stroke_px if old_stroke_px is _BOND_STROKE_UNSET else old_stroke_px
        )
        self._new_stroke = new_stroke_px

    def redo(self) -> None:
        """Aplica el nuevo grosor."""
        self._model.update_bond(self._bond_id, stroke_px=self._new_stroke)
        self._view.update_bond_item(self._bond_id)

    def undo(self) -> None:
        """Revierte al grosor anterior."""
        self._model.update_bond(self._bond_id, stroke_px=self._old_stroke)
        self._view.update_bond_item(self._bond_id)

class ChangeBondColorCommand(QUndoCommand):
    """Comando para cambiar el color de un enlace."""

    def __init__(
        self,
        model: MolGraph,
        view,
        bond_id: int,
        new_color: Optional[str],
    ) -> None:
        """Inicializa el comando de cambio de color."""
        super().__init__("Change bond color")
        self._model = model
        self._view = view
        self._bond_id = bond_id
        bond = model.get_bond(bond_id)
        self._old_color = bond.color
        self._new_color = new_color

    def redo(self) -> None:
        """Aplica el nuevo color."""
        self._model.update_bond(self._bond_id, color=self._new_color)
        self._view.update_bond_item(self._bond_id)

    def undo(self) -> None:
        """Revierte al color anterior."""
        self._model.update_bond(self._bond_id, color=self._old_color)
        self._view.update_bond_item(self._bond_id)

class ChangeDoubleBondOrientationCommand(QUndoCommand):
    """Comando para fijar orientación manual de línea pi en un doble enlace."""

    def __init__(
        self,
        model: MolGraph,
        view,
        bond_id: int,
        old_sign: Optional[int],
        new_sign: Optional[int],
    ) -> None:
        super().__init__("Toggle double bond orientation")
        self._model = model
        self._view = view
        self._bond_id = int(bond_id)
        self._old_sign = old_sign if old_sign in {-1, 1} else None
        self._new_sign = new_sign if new_sign in {-1, 1} else None

    def _apply(self, sign: Optional[int]) -> None:
        """Aplica orientación en modelo e item para mantener consistencia visual."""
        if self._bond_id not in self._model.bonds:
            return
        self._model.update_bond(self._bond_id, pi_offset_sign=sign)
        item = self._view.bond_items.get(self._bond_id)
        if item is not None and hasattr(item, "set_manual_pi_offset"):
            item.set_manual_pi_offset(sign)
        else:
            self._view.update_bond_item(self._bond_id)

    def redo(self) -> None:
        self._apply(self._new_sign)

    def undo(self) -> None:
        self._apply(self._old_sign)

class AddRingCommand(QUndoCommand):
    """Comando para añadir un anillo (posiblemente con átomos nuevos)."""

    def __init__(
        self,
        model: MolGraph,
        view,
        vertices: List[Tuple[Optional[int], float, float]],
        edges: List[Tuple],
        element: str = "C",
    ) -> None:
        """Inicializa el comando de creación de anillo."""
        super().__init__("Add ring")
        self._model = model
        self._view = view
        self._vertices = vertices
        self._edges = edges
        self._element = element
        self._created_atom_ids: List[Optional[int]] = []
        self._created_bonds = []
        self._ring_id: Optional[int] = None
        self._updated_existing: List[Tuple[int, int, bool, int]] = []

    def redo(self) -> None:
        """Crea o reintroduce el anillo con sus enlaces."""
        if not self._created_atom_ids:
            self._created_atom_ids = [v[0] for v in self._vertices]

        if self._ring_id is None:
            self._ring_id = self._view.allocate_ring_id()

        atom_ids: List[int] = []
        for idx, (existing_id, x, y) in enumerate(self._vertices):
            if existing_id is not None:
                atom_ids.append(existing_id)
                continue
            atom_id = self._created_atom_ids[idx]
            if atom_id is None:
                atom = self._model.add_atom(
                    self._element,
                    x,
                    y,
                    is_explicit=_default_is_explicit(self._element),
                )
                self._created_atom_ids[idx] = atom.id
            else:
                atom = self._model.add_atom(
                    self._element,
                    x,
                    y,
                    atom_id=atom_id,
                    is_explicit=_default_is_explicit(self._element),
                )
            self._view.add_atom_item(atom)
            atom_ids.append(self._created_atom_ids[idx])

        if self._ring_id is not None:
            xs = [self._model.get_atom(aid).x for aid in atom_ids]
            ys = [self._model.get_atom(aid).y for aid in atom_ids]
            center = (sum(xs) / len(xs), sum(ys) / len(ys))
            self._view.register_ring_center(self._ring_id, center)

        if not self._created_bonds:
            for edge in self._edges:
                i, j, order, style, stereo = edge[:5]
                is_aromatic = edge[5] if len(edge) > 5 else False
                a1_id = atom_ids[i]
                a2_id = atom_ids[j]
                existing = self._model.find_bond_between(a1_id, a2_id)
                if existing is not None:
                    if is_aromatic:
                        if not any(bid == existing.id for bid, _, _, _ in self._updated_existing):
                            self._updated_existing.append(
                                (existing.id, existing.order, existing.is_aromatic, order)
                            )
                        self._model.update_bond(
                            existing.id,
                            order=order,
                            style=existing.style,
                            stereo=existing.stereo,
                            is_aromatic=True,
                        )
                        self._view.update_bond_item(existing.id)
                    continue
                bond = self._model.add_bond(
                    a1_id,
                    a2_id,
                    order,
                    style=style,
                    stereo=stereo,
                    is_aromatic=is_aromatic,
                    display_order=None,
                    ring_id=self._ring_id,
                )
                self._created_bonds.append(replace(bond))
                self._view.add_bond_item(bond)
                self._view.update_bond_item(bond.id)
        else:
            if self._ring_id is not None:
                xs = [self._model.get_atom(aid).x for aid in atom_ids]
                ys = [self._model.get_atom(aid).y for aid in atom_ids]
                center = (sum(xs) / len(xs), sum(ys) / len(ys))
                self._view.register_ring_center(self._ring_id, center)
            for bond_id, old_order, old_aromatic, new_order in self._updated_existing:
                bond = self._model.get_bond(bond_id)
                self._model.update_bond(
                    bond_id,
                    order=new_order,
                    style=bond.style,
                    stereo=bond.stereo,
                    is_aromatic=True,
                )
                self._view.update_bond_item(bond_id)
            for bond in self._created_bonds:
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
                )
                self._view.add_bond_item(bond)
                self._view.update_bond_item(bond.id)

    def undo(self) -> None:
        """Elimina el anillo y restaura enlaces previos."""
        for bond in list(self._created_bonds):
            if bond.id in self._model.bonds:
                self._model.remove_bond(bond.id)
                self._view.remove_bond_item(bond.id)
        for idx, (existing_id, _x, _y) in enumerate(self._vertices):
            if existing_id is not None:
                continue
            atom_id = self._created_atom_ids[idx]
            if atom_id is not None and atom_id in self._model.atoms:
                self._model.remove_atom(atom_id)
                self._view.remove_atom_item(atom_id)
        if self._ring_id is not None:
            self._view.unregister_ring_center(self._ring_id)
        for bond_id, old_order, old_aromatic, new_order in self._updated_existing:
            if bond_id in self._model.bonds:
                bond = self._model.get_bond(bond_id)
                self._model.update_bond(
                    bond_id,
                    order=old_order,
                    style=bond.style,
                    stereo=bond.stereo,
                    is_aromatic=old_aromatic,
                )
                self._view.update_bond_item(bond_id)

class AddChainCommand(QUndoCommand):
    """Comando para añadir una cadena lineal de átomos."""

    def __init__(
        self,
        model: MolGraph,
        view,
        anchor_id: Optional[int],
        positions: List[Tuple[float, float]],
        element: str = "C",
        anchor_position: Optional[Tuple[float, float]] = None,
    ) -> None:
        """Inicializa el comando de creación de cadena."""
        super().__init__("Add chain")
        self._model = model
        self._view = view
        self._anchor_id = anchor_id
        self._anchor_position = anchor_position
        self._positions = positions
        self._element = element
        self._created_atom_ids: List[Optional[int]] = []
        self._created_bond_ids: List[Optional[int]] = []
        self._created_anchor_id: Optional[int] = None

    def redo(self) -> None:
        """Crea o reintroduce la cadena."""
        if not self._created_atom_ids:
            self._created_atom_ids = [None for _ in self._positions]
        if not self._created_bond_ids:
            self._created_bond_ids = [None for _ in self._positions]

        prev_id = self._anchor_id
        if prev_id is None:
            if self._anchor_position is None:
                raise RuntimeError("Anchor position required to place a free chain")
            if self._created_anchor_id is None:
                anchor_atom = self._model.add_atom(
                    self._element,
                    self._anchor_position[0],
                    self._anchor_position[1],
                    is_explicit=_default_is_explicit(self._element),
                )
                self._created_anchor_id = anchor_atom.id
            else:
                anchor_atom = self._model.add_atom(
                    self._element,
                    self._anchor_position[0],
                    self._anchor_position[1],
                    atom_id=self._created_anchor_id,
                    is_explicit=_default_is_explicit(self._element),
                )
            self._view.add_atom_item(anchor_atom)
            prev_id = self._created_anchor_id
        for idx, (x, y) in enumerate(self._positions):
            atom_id = self._created_atom_ids[idx]
            if atom_id is None:
                atom = self._model.add_atom(
                    self._element,
                    x,
                    y,
                    is_explicit=_default_is_explicit(self._element),
                )
                self._created_atom_ids[idx] = atom.id
            else:
                atom = self._model.add_atom(
                    self._element,
                    x,
                    y,
                    atom_id=atom_id,
                    is_explicit=_default_is_explicit(self._element),
                )
            self._view.add_atom_item(atom)
            bond_id = self._created_bond_ids[idx]
            bond = self._model.add_bond(
                prev_id,
                self._created_atom_ids[idx],
                bond_id=bond_id,
            )
            self._created_bond_ids[idx] = bond.id
            self._view.add_bond_item(bond)
            prev_id = self._created_atom_ids[idx]

    def undo(self) -> None:
        """Elimina la cadena y sus enlaces."""
        for bond_id in list(self._created_bond_ids):
            if bond_id in self._model.bonds:
                self._model.remove_bond(bond_id)
                self._view.remove_bond_item(bond_id)

        for atom_id in list(self._created_atom_ids):
            if atom_id is not None and atom_id in self._model.atoms:
                self._model.remove_atom(atom_id)
                self._view.remove_atom_item(atom_id)
        if self._anchor_id is None and self._created_anchor_id is not None:
            if self._created_anchor_id in self._model.atoms:
                self._model.remove_atom(self._created_anchor_id)
                self._view.remove_atom_item(self._created_anchor_id)
