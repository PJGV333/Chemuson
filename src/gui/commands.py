from __future__ import annotations

from dataclasses import replace
from typing import Dict, Iterable, List, Optional, Tuple

from PyQt6.QtGui import QUndoCommand

from core.model import BondStyle, BondStereo, MolGraph


class AddAtomCommand(QUndoCommand):
    def __init__(self, model: MolGraph, view, element: str, x: float, y: float) -> None:
        super().__init__("Add atom")
        self._model = model
        self._view = view
        self._element = element
        self._x = x
        self._y = y
        self._atom_id: Optional[int] = None

    def redo(self) -> None:
        if self._atom_id is None:
            atom = self._model.add_atom(self._element, self._x, self._y)
            self._atom_id = atom.id
        else:
            atom = self._model.add_atom(self._element, self._x, self._y, atom_id=self._atom_id)
        self._view.add_atom_item(atom)

    def undo(self) -> None:
        atom, removed_bonds = self._model.remove_atom(self._atom_id)
        for bond in removed_bonds:
            self._view.remove_bond_item(bond.id)
        self._view.remove_atom_item(atom.id)

    @property
    def atom_id(self) -> Optional[int]:
        return self._atom_id


class ChangeAtomCommand(QUndoCommand):
    def __init__(self, model: MolGraph, view, atom_id: int, new_element: str) -> None:
        super().__init__("Change atom")
        self._model = model
        self._view = view
        self._atom_id = atom_id
        self._old_element = model.get_atom(atom_id).element
        self._new_element = new_element

    def redo(self) -> None:
        self._model.update_atom_element(self._atom_id, self._new_element)
        self._view.update_atom_item_element(self._atom_id, self._new_element)

    def undo(self) -> None:
        self._model.update_atom_element(self._atom_id, self._old_element)
        self._view.update_atom_item_element(self._atom_id, self._old_element)


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
        self._ring_id = ring_id
        self._bond_id: Optional[int] = None
        self._new_atom_element = new_atom_element
        self._new_atom_pos = new_atom_pos
        self._created_atom_id: Optional[int] = None

    def redo(self) -> None:
        if self._a2_id is None:
            if self._created_atom_id is None:
                if self._new_atom_element is None or self._new_atom_pos is None:
                    return
                atom = self._model.add_atom(
                    self._new_atom_element,
                    self._new_atom_pos[0],
                    self._new_atom_pos[1],
                )
                self._created_atom_id = atom.id
            else:
                atom = self._model.add_atom(
                    self._new_atom_element,
                    self._new_atom_pos[0],
                    self._new_atom_pos[1],
                    atom_id=self._created_atom_id,
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
                ring_id=self._ring_id,
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
                ring_id=self._ring_id,
            )
        self._view.add_bond_item(bond)

    def undo(self) -> None:
        bond = self._model.remove_bond(self._bond_id)
        self._view.remove_bond_item(bond.id)
        if self._created_atom_id is not None:
            atom, removed_bonds = self._model.remove_atom(self._created_atom_id)
            for removed in removed_bonds:
                self._view.remove_bond_item(removed.id)
            self._view.remove_atom_item(atom.id)
            self._a2_id = None


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


class DeleteSelectionCommand(QUndoCommand):
    def __init__(
        self,
        model: MolGraph,
        view,
        atom_ids: Iterable[int],
        bond_ids: Iterable[int],
    ) -> None:
        super().__init__("Delete selection")
        self._model = model
        self._view = view
        self._atom_ids = sorted(set(atom_ids))
        self._bond_ids = sorted(set(bond_ids))
        self._removed_atoms = []
        self._removed_bonds = []

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

        for bond in list(self._removed_bonds):
            if bond.id in self._model.bonds:
                self._model.remove_bond(bond.id)
                self._view.remove_bond_item(bond.id)
        for atom in list(self._removed_atoms):
            if atom.id in self._model.atoms:
                self._model.remove_atom(atom.id)
                self._view.remove_atom_item(atom.id)

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
            )
            self._view.add_bond_item(bond)


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
                atom = self._model.add_atom(self._element, x, y)
                self._created_atom_ids[idx] = atom.id
            else:
                atom = self._model.add_atom(self._element, x, y, atom_id=atom_id)
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
        anchor_id: int,
        positions: List[Tuple[float, float]],
        element: str = "C",
    ) -> None:
        super().__init__("Add chain")
        self._model = model
        self._view = view
        self._anchor_id = anchor_id
        self._positions = positions
        self._element = element
        self._created_atom_ids: List[Optional[int]] = []
        self._created_bond_ids: List[Optional[int]] = []

    def redo(self) -> None:
        if not self._created_atom_ids:
            self._created_atom_ids = [None for _ in self._positions]
        if not self._created_bond_ids:
            self._created_bond_ids = [None for _ in self._positions]

        prev_id = self._anchor_id
        for idx, (x, y) in enumerate(self._positions):
            atom_id = self._created_atom_ids[idx]
            if atom_id is None:
                atom = self._model.add_atom(self._element, x, y)
                self._created_atom_ids[idx] = atom.id
            else:
                atom = self._model.add_atom(self._element, x, y, atom_id=atom_id)
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
