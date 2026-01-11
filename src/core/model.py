from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set


@dataclass
class Atom:
    id: int
    element: str
    x: float
    y: float
    charge: int = 0
    isotope: Optional[int] = None
    explicit_h: Optional[int] = None
    mapping: Optional[int] = None
    is_query: bool = False


@dataclass
class Bond:
    id: int
    a1_id: int
    a2_id: int
    order: int = 1
    stereo: Optional[str] = None
    is_aromatic: bool = False
    is_query: bool = False


@dataclass
class ChemState:
    active_tool: str = "atom_C"
    bond_order: int = 1
    default_element: str = "C"
    selected_atoms: Set[int] = field(default_factory=set)
    selected_bonds: Set[int] = field(default_factory=set)


class MolGraph:
    def __init__(self) -> None:
        self.atoms: Dict[int, Atom] = {}
        self.bonds: Dict[int, Bond] = {}
        self._next_atom_id = 1
        self._next_bond_id = 1

    def add_atom(
        self,
        element: str,
        x: float,
        y: float,
        atom_id: Optional[int] = None,
        charge: int = 0,
        isotope: Optional[int] = None,
        explicit_h: Optional[int] = None,
        mapping: Optional[int] = None,
        is_query: bool = False,
    ) -> Atom:
        if atom_id is None:
            atom_id = self._next_atom_id
            self._next_atom_id += 1
        else:
            self._next_atom_id = max(self._next_atom_id, atom_id + 1)
        atom = Atom(
            id=atom_id,
            element=element,
            x=x,
            y=y,
            charge=charge,
            isotope=isotope,
            explicit_h=explicit_h,
            mapping=mapping,
            is_query=is_query,
        )
        self.atoms[atom_id] = atom
        return atom

    def remove_atom(self, atom_id: int) -> tuple[Atom, List[Bond]]:
        atom = self.atoms.pop(atom_id)
        removed_bonds: List[Bond] = []
        for bond_id, bond in list(self.bonds.items()):
            if bond.a1_id == atom_id or bond.a2_id == atom_id:
                removed_bonds.append(self.bonds.pop(bond_id))
        return atom, removed_bonds

    def add_bond(
        self,
        a1_id: int,
        a2_id: int,
        order: int = 1,
        bond_id: Optional[int] = None,
        stereo: Optional[str] = None,
        is_aromatic: bool = False,
        is_query: bool = False,
    ) -> Bond:
        if bond_id is None:
            bond_id = self._next_bond_id
            self._next_bond_id += 1
        else:
            self._next_bond_id = max(self._next_bond_id, bond_id + 1)
        bond = Bond(
            id=bond_id,
            a1_id=a1_id,
            a2_id=a2_id,
            order=order,
            stereo=stereo,
            is_aromatic=is_aromatic,
            is_query=is_query,
        )
        self.bonds[bond_id] = bond
        return bond

    def remove_bond(self, bond_id: int) -> Bond:
        return self.bonds.pop(bond_id)

    def get_atom(self, atom_id: int) -> Atom:
        return self.atoms[atom_id]

    def get_bond(self, bond_id: int) -> Bond:
        return self.bonds[bond_id]

    def find_bond_between(self, a1_id: int, a2_id: int) -> Optional[Bond]:
        for bond in self.bonds.values():
            if {bond.a1_id, bond.a2_id} == {a1_id, a2_id}:
                return bond
        return None

    def update_atom_position(self, atom_id: int, x: float, y: float) -> None:
        atom = self.atoms[atom_id]
        atom.x = x
        atom.y = y

    def clear(self) -> None:
        self.atoms.clear()
        self.bonds.clear()
        self._next_atom_id = 1
        self._next_bond_id = 1
