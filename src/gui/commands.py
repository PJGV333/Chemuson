from __future__ import annotations

from typing import Optional

from PyQt6.QtGui import QUndoCommand

from core.model import MolGraph


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


class AddBondCommand(QUndoCommand):
    def __init__(self, model: MolGraph, view, a1_id: int, a2_id: int, order: int = 1) -> None:
        super().__init__("Add bond")
        self._model = model
        self._view = view
        self._a1_id = a1_id
        self._a2_id = a2_id
        self._order = order
        self._bond_id: Optional[int] = None

    def redo(self) -> None:
        if self._bond_id is None:
            bond = self._model.add_bond(self._a1_id, self._a2_id, self._order)
            self._bond_id = bond.id
        else:
            bond = self._model.add_bond(
                self._a1_id,
                self._a2_id,
                self._order,
                bond_id=self._bond_id,
            )
        self._view.add_bond_item(bond)

    def undo(self) -> None:
        bond = self._model.remove_bond(self._bond_id)
        self._view.remove_bond_item(bond.id)
