from __future__ import annotations

from dataclasses import replace
from typing import Dict, Iterable, List, Optional, Tuple

from PyQt6.QtCore import QPointF, QRectF
from PyQt6.QtGui import QUndoCommand

from core.model import BondStyle, BondStereo, MolGraph
from gui.geom import angle_deg, angle_distance_deg, endpoint_from_angle_len

_IMPLICIT_ELEMENTS = {"C"}
_ANCHOR_UNSET = object()


def _default_is_explicit(element: str) -> bool:
    return element not in _IMPLICIT_ELEMENTS


def _atom_degree(model: MolGraph, atom_id: int) -> int:
    return sum(
        1
        for bond in model.bonds.values()
        if bond.a1_id == atom_id or bond.a2_id == atom_id
    )


def _neighbor_angles_deg(model: MolGraph, atom_id: int) -> list[float]:
    anchor = model.get_atom(atom_id)
    origin = QPointF(anchor.x, anchor.y)
    angles: list[float] = []
    for bond in model.bonds.values():
        if bond.a1_id == atom_id:
            other = model.get_atom(bond.a2_id)
        elif bond.a2_id == atom_id:
            other = model.get_atom(bond.a1_id)
        else:
            continue
        angles.append(angle_deg(origin, QPointF(other.x, other.y)))
    return angles


def _select_hydrogen_angles(existing_angles_deg: list[float], count: int) -> list[float]:
    base_angles = [0.0, 120.0, 240.0]
    if count <= 0:
        return []
    if not existing_angles_deg:
        return base_angles[:count]

    def min_distance(angle: float) -> float:
        return min(angle_distance_deg(angle, existing) for existing in existing_angles_deg)

    ordered = sorted(base_angles, key=min_distance, reverse=True)
    return ordered[:count]


def _bond_length_from_view(view) -> float:
    return getattr(getattr(view, "state", None), "bond_length", 40.0)


def _remove_hydrogen_specs(model: MolGraph, view, specs: list[tuple[int, float, float, int]]) -> None:
    for _atom_id, _x, _y, bond_id in specs:
        if bond_id in model.bonds:
            model.remove_bond(bond_id)
            view.remove_bond_item(bond_id)
    for atom_id, _x, _y, _bond_id in specs:
        if atom_id in model.atoms:
            model.remove_atom(atom_id)
            view.remove_atom_item(atom_id)


def _readd_hydrogen_specs(
    model: MolGraph,
    view,
    anchor_id: int,
    specs: list[tuple[int, float, float, int]],
) -> None:
    for atom_id, x, y, _bond_id in specs:
        atom = model.add_atom("H", x, y, atom_id=atom_id, is_explicit=True)
        view.add_atom_item(atom)
    for atom_id, _x, _y, bond_id in specs:
        bond = model.add_bond(anchor_id, atom_id, bond_id=bond_id, style=BondStyle.PLAIN)
        view.add_bond_item(bond)


def _create_hydrogen_specs(
    model: MolGraph,
    view,
    anchor_id: int,
    angles_deg: list[float],
    bond_length: float,
) -> list[tuple[int, float, float, int]]:
    anchor = model.get_atom(anchor_id)
    origin = QPointF(anchor.x, anchor.y)
    specs: list[tuple[int, float, float, int]] = []
    for angle_deg_value in angles_deg:
        pos = endpoint_from_angle_len(origin, angle_deg_value, bond_length)
        atom = model.add_atom("H", pos.x(), pos.y(), is_explicit=True)
        view.add_atom_item(atom)
        bond = model.add_bond(anchor_id, atom.id, style=BondStyle.PLAIN)
        view.add_bond_item(bond)
        specs.append((atom.id, pos.x(), pos.y(), bond.id))
    return specs


def _collect_attached_hydrogens(
    model: MolGraph,
    atom_id: int,
) -> tuple[list, list]:
    removed_atoms = []
    removed_bonds = []
    for bond in model.bonds.values():
        if bond.a1_id == atom_id:
            other_id = bond.a2_id
        elif bond.a2_id == atom_id:
            other_id = bond.a1_id
        else:
            continue
        other = model.atoms.get(other_id)
        if other is None or other.element != "H":
            continue
        if _atom_degree(model, other_id) != 1:
            continue
        removed_atoms.append(replace(other))
        removed_bonds.append(replace(bond))
    return removed_atoms, removed_bonds


class AddAtomCommand(QUndoCommand):
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
        explicit_h: Optional[int] = None,
        mapping: Optional[int] = None,
        is_query: bool = False,
        anchor_override: Optional[str] = None,
        auto_hydrogens: bool = True,
        expected_bonds: int = 0,
    ) -> None:
        super().__init__("Add atom")
        self._model = model
        self._view = view
        self._element = element
        self._x = x
        self._y = y
        self._is_explicit = is_explicit
        self._charge = charge
        self._isotope = isotope
        self._explicit_h = explicit_h
        self._mapping = mapping
        self._is_query = is_query
        self._anchor_override = anchor_override
        self._auto_hydrogens = auto_hydrogens
        self._expected_bonds = expected_bonds
        self._atom_id: Optional[int] = None
        self._hydrogen_specs: list[tuple[int, float, float, int]] = []

    def redo(self) -> None:
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
                explicit_h=self._explicit_h,
                mapping=self._mapping,
                is_query=self._is_query,
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
                explicit_h=self._explicit_h,
                mapping=self._mapping,
                is_query=self._is_query,
            )
        self._view.add_atom_item(atom)
        if self._anchor_override:
            self._view.set_anchor_override(atom.id, self._anchor_override)
            self._view._refresh_atom_label(atom.id)

    def undo(self) -> None:
        if self._hydrogen_specs:
            _remove_hydrogen_specs(self._model, self._view, self._hydrogen_specs)
        atom, removed_bonds = self._model.remove_atom(self._atom_id)
        for bond in removed_bonds:
            self._view.remove_bond_item(bond.id)
        self._view.remove_atom_item(atom.id)

    @property
    def atom_id(self) -> Optional[int]:
        return self._atom_id


class ChangeAtomCommand(QUndoCommand):
    def __init__(
        self,
        model: MolGraph,
        view,
        atom_id: int,
        new_element: str,
        anchor_override=_ANCHOR_UNSET,
    ) -> None:
        super().__init__("Change atom")
        self._model = model
        self._view = view
        self._atom_id = atom_id
        self._old_element = model.get_atom(atom_id).element
        self._old_is_explicit = model.get_atom(atom_id).is_explicit
        self._new_element = new_element
        self._new_is_explicit = _default_is_explicit(new_element)
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
        self._view.update_atom_item_element(
            self._atom_id,
            self._new_element,
            is_explicit=self._new_is_explicit,
        )

    def undo(self) -> None:
        if self._anchor_override is not _ANCHOR_UNSET:
            self._view.set_anchor_override(self._atom_id, self._old_anchor_override)
        self._model.update_atom_element(
            self._atom_id,
            self._old_element,
            is_explicit=self._old_is_explicit,
        )
        self._view.update_atom_item_element(
            self._atom_id,
            self._old_element,
            is_explicit=self._old_is_explicit,
        )
        
        # Restore hydrogens
        if self._removed_hydrogen_specs:
            _readd_hydrogen_specs(self._model, self._view, self._atom_id, self._removed_hydrogen_specs)

    def _check_and_remove_hydrogens(self) -> None:
        from core.model import VALENCE_MAP
        
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
    def __init__(self, model: MolGraph, view, atom_id: int, new_charge: int) -> None:
        super().__init__("Change charge")
        self._model = model
        self._view = view
        self._atom_id = atom_id
        self._old_charge = model.get_atom(atom_id).charge
        self._new_charge = new_charge

    def redo(self) -> None:
        self._model.update_atom_charge(self._atom_id, self._new_charge)
        self._view.update_atom_item_charge(self._atom_id, self._new_charge)

    def undo(self) -> None:
        self._model.update_atom_charge(self._atom_id, self._old_charge)
        self._view.update_atom_item_charge(self._atom_id, self._old_charge)


class AddBondCommand(QUndoCommand):
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
        new_atom_element: Optional[str] = None,
        new_atom_pos: Optional[Tuple[float, float]] = None,
    ) -> None:
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
        self._bond_id: Optional[int] = None
        self._new_atom_element = new_atom_element
        self._new_atom_pos = new_atom_pos
        self._created_atom_id: Optional[int] = None
        self._hydrogen_specs: list[tuple[int, float, float, int]] = []
        self._demoted_explicit_atoms: Optional[list[int]] = None

    def _should_demote_explicit_carbon(self, atom_id: int) -> bool:
        atom = self._model.get_atom(atom_id)
        if atom.element != "C" or not atom.is_explicit:
            return False
        if getattr(getattr(self._view, "state", None), "show_implicit_carbons", True):
            return False
        if atom.charge != 0 or atom.isotope is not None or atom.explicit_h is not None:
            return False
        if atom.mapping is not None or atom.is_query:
            return False
        if _atom_degree(self._model, atom_id) <= 0:
            return False
        return True

    def redo(self) -> None:
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

    def undo(self) -> None:
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


class ChangeBondCommand(QUndoCommand):
    def __init__(
        self,
        model: MolGraph,
        view,
        bond_id: int,
        new_order: Optional[int] = None,
        new_style: Optional[BondStyle] = None,
        new_stereo: Optional[BondStereo] = None,
        new_is_aromatic: Optional[bool] = None,
    ) -> None:
        super().__init__("Change bond")
        self._model = model
        self._view = view
        self._bond_id = bond_id
        bond = model.get_bond(bond_id)
        self._old_order = bond.order
        self._old_style = bond.style
        self._old_stereo = bond.stereo
        self._old_is_aromatic = bond.is_aromatic
        self._new_order = new_order if new_order is not None else bond.order
        self._new_style = new_style if new_style is not None else bond.style
        self._new_stereo = new_stereo if new_stereo is not None else bond.stereo
        self._new_is_aromatic = (
            new_is_aromatic if new_is_aromatic is not None else bond.is_aromatic
        )

    def redo(self) -> None:
        self._model.update_bond(
            self._bond_id,
            order=self._new_order,
            style=self._new_style,
            stereo=self._new_stereo,
            is_aromatic=self._new_is_aromatic,
        )
        self._view.update_bond_item(self._bond_id)

    def undo(self) -> None:
        self._model.update_bond(
            self._bond_id,
            order=self._old_order,
            style=self._old_style,
            stereo=self._old_stereo,
            is_aromatic=self._old_is_aromatic,
        )
        self._view.update_bond_item(self._bond_id)


class ChangeBondLengthCommand(QUndoCommand):
    def __init__(self, model: MolGraph, view, bond_id: int, new_length: Optional[float]) -> None:
        super().__init__("Change bond length")
        self._model = model
        self._view = view
        self._bond_id = bond_id
        bond = model.get_bond(bond_id)
        self._old_length = bond.length_px
        self._new_length = new_length

    def redo(self) -> None:
        self._model.update_bond_length(self._bond_id, self._new_length)
        self._view.update_bond_item(self._bond_id)

    def undo(self) -> None:
        self._model.update_bond_length(self._bond_id, self._old_length)
        self._view.update_bond_item(self._bond_id)

class MoveAtomsCommand(QUndoCommand):
    def __init__(
        self,
        model: MolGraph,
        view,
        before: Dict[int, Tuple[float, float]],
        after: Dict[int, Tuple[float, float]],
        skip_first_redo: bool = False,
    ) -> None:
        super().__init__("Move atoms")
        self._model = model
        self._view = view
        self._before = before
        self._after = after
        self._skip_first_redo = skip_first_redo
        self._first_redo = True

    def redo(self) -> None:
        if self._skip_first_redo and self._first_redo:
            self._first_redo = False
            return
        self._apply_positions(self._after)

    def undo(self) -> None:
        self._apply_positions(self._before)

    def _apply_positions(self, positions: Dict[int, Tuple[float, float]]) -> None:
        for atom_id, (x, y) in positions.items():
            self._model.update_atom_position(atom_id, x, y)
            self._view.update_atom_item(atom_id, x, y)
        self._view.update_bond_items_for_atoms(set(positions.keys()))


class MoveTextItemsCommand(QUndoCommand):
    def __init__(self, view, before: dict, after: dict):
        super().__init__("Move text items")
        self._view = view
        self._before = before # {item: (pos, rotation)}
        self._after = after   # {item: (pos, rotation)}

    def redo(self):
        for item, (pos, rot) in self._after.items():
            item.setPos(pos)
            item.setRotation(rot)
        self._view._update_selection_overlay()

    def undo(self):
        for item, (pos, rot) in self._before.items():
            item.setPos(pos)
            item.setRotation(rot)
        self._view._update_selection_overlay()


class MoveArrowItemsCommand(QUndoCommand):
    def __init__(self, view, before: dict, after: dict):
        super().__init__("Move arrows")
        self._view = view
        self._before = before  # {item: (start, end)}
        self._after = after    # {item: (start, end)}

    def redo(self) -> None:
        for item, (start, end) in self._after.items():
            item.update_positions(start, end)
        self._view._update_selection_overlay()

    def undo(self) -> None:
        for item, (start, end) in self._before.items():
            item.update_positions(start, end)
        self._view._update_selection_overlay()


class MoveBracketItemsCommand(QUndoCommand):
    def __init__(self, view, before: dict, after: dict):
        super().__init__("Move brackets")
        self._view = view
        self._before = before  # {item: QRectF}
        self._after = after    # {item: QRectF}

    def redo(self) -> None:
        for item, rect in self._after.items():
            item.set_rect(rect)
        self._view._update_selection_overlay()

    def undo(self) -> None:
        for item, rect in self._before.items():
            item.set_rect(rect)
        self._view._update_selection_overlay()


class DeleteSelectionCommand(QUndoCommand):
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
    ) -> None:
        super().__init__("Delete selection")
        self._model = model
        self._view = view
        self._atom_ids = sorted(set(atom_ids))
        self._bond_ids = sorted(set(bond_ids))
        self._arrow_items = list(arrow_items)
        self._bracket_items = list(bracket_items)
        self._text_items = list(text_items)
        self._wavy_items = list(wavy_items)
        self._removed_atoms = []
        self._removed_bonds = []
        self._removed_arrows = []
        self._removed_brackets = []
        self._removed_texts = []
        self._removed_wavy = []

    def redo(self) -> None:
        if not self._removed_atoms and not self._removed_bonds:
            for atom_id in self._atom_ids:
                if atom_id in self._model.atoms:
                    self._removed_atoms.append(replace(self._model.atoms[atom_id]))
            for bond in self._model.bonds.values():
                if (
                    bond.id in self._bond_ids
                    or bond.a1_id in self._atom_ids
                    or bond.a2_id in self._atom_ids
                ):
                    self._removed_bonds.append(replace(bond))
            for item in self._arrow_items:
                self._removed_arrows.append(
                    (item, item.start_point(), item.end_point(), item.kind())
                )
            for item in self._bracket_items:
                self._removed_brackets.append(
                    (item, item.base_rect(), item._padding, item._kind)
                )
            for item in self._text_items:
                self._removed_texts.append(item)
            for item in self._wavy_items:
                self._removed_wavy.append(item)

        for bond in list(self._removed_bonds):
            if bond.id in self._model.bonds:
                self._model.remove_bond(bond.id)
                self._view.remove_bond_item(bond.id)
        for atom in list(self._removed_atoms):
            if atom.id in self._model.atoms:
                self._model.remove_atom(atom.id)
                self._view.remove_atom_item(atom.id)
        for item, _start, _end, _kind in list(self._removed_arrows):
            self._view.remove_arrow_item(item)
        for item, _rect, _padding, _kind in list(self._removed_brackets):
            self._view.remove_bracket_item(item)
        for item in self._removed_texts:
            self._view.remove_text_item(item)
        for item in self._removed_wavy:
            self._view.remove_wavy_anchor_item(item)

    def undo(self) -> None:
        for atom in self._removed_atoms:
            self._model.add_atom(
                atom.element,
                atom.x,
                atom.y,
                atom_id=atom.id,
                charge=atom.charge,
                isotope=atom.isotope,
                explicit_h=atom.explicit_h,
                mapping=atom.mapping,
                is_query=atom.is_query,
                is_explicit=atom.is_explicit,
            )
            self._view.add_atom_item(atom)
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
            )
            self._view.add_bond_item(bond)
        for item, start, end, kind in self._removed_arrows:
            self._view.readd_arrow_item(item, start, end, kind)
        for item, rect, padding, kind in self._removed_brackets:
            self._view.readd_bracket_item(item, rect, kind, padding=padding)
        for item in self._removed_texts:
            self._view.readd_text_item(item)
        for item in self._removed_wavy:
            self._view.readd_wavy_anchor_item(item)


class AddArrowCommand(QUndoCommand):
    def __init__(self, view, start: QPointF, end: QPointF, kind: str) -> None:
        super().__init__("Add arrow")
        self._view = view
        self._start = QPointF(start)
        self._end = QPointF(end)
        self._kind = kind
        self._item = None

    def redo(self) -> None:
        if self._item is None:
            self._item = self._view.add_arrow_item(self._start, self._end, self._kind)
        else:
            self._view.readd_arrow_item(self._item, self._start, self._end, self._kind)

    def undo(self) -> None:
        if self._item is not None:
            self._view.remove_arrow_item(self._item)


class AddBracketCommand(QUndoCommand):
    def __init__(self, view, rect: QRectF, kind: str) -> None:
        super().__init__("Add brackets")
        self._view = view
        self._rect = QRectF(rect)
        self._kind = kind
        self._item = None

    def redo(self) -> None:
        if self._item is None:
            self._item = self._view.add_bracket_item(self._rect, self._kind)
        else:
            self._view.readd_bracket_item(self._item, self._rect, self._kind)

    def undo(self) -> None:
        if self._item is not None:
            self._view.remove_bracket_item(self._item)


class AddTextItemCommand(QUndoCommand):
    def __init__(self, view, item) -> None:
        super().__init__("Add text")
        self._view = view
        self._item = item

    def redo(self) -> None:
        if self._item is not None:
            self._view.readd_text_item(self._item)

    def undo(self) -> None:
        if self._item is not None:
            self._view.remove_text_item(self._item)


class AddWavyAnchorCommand(QUndoCommand):
    def __init__(self, view, item) -> None:
        super().__init__("Add wavy anchor")
        self._view = view
        self._item = item

    def redo(self) -> None:
        if self._item is not None:
            self._view.readd_wavy_anchor_item(self._item)

    def undo(self) -> None:
        if self._item is not None:
            self._view.remove_wavy_anchor_item(self._item)


class AddRingCommand(QUndoCommand):
    def __init__(
        self,
        model: MolGraph,
        view,
        vertices: List[Tuple[Optional[int], float, float]],
        edges: List[Tuple],
        element: str = "C",
    ) -> None:
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
                )
                self._view.add_bond_item(bond)
                self._view.update_bond_item(bond.id)

    def undo(self) -> None:
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
    def __init__(
        self,
        model: MolGraph,
        view,
        anchor_id: Optional[int],
        positions: List[Tuple[float, float]],
        element: str = "C",
        anchor_position: Optional[Tuple[float, float]] = None,
    ) -> None:
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
