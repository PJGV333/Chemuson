from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Set


class BondStyle(str, Enum):
    PLAIN = "plain"
    WEDGE = "wedge"
    HASHED = "hashed"
    WAVY = "wavy"
    INTERACTION = "interaction"


class BondStereo(str, Enum):
    NONE = "none"
    UP = "up"
    DOWN = "down"
    EITHER = "either"


VALENCE_MAP = {
    "H": 1,
    "C": 4,
    "N": 3,
    "O": 2,
    "S": 2,
    "P": 3,
    "F": 1,
    "Cl": 1,
    "Br": 1,
    "I": 1,
}


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
    is_explicit: bool = False


@dataclass
class Bond:
    id: int
    a1_id: int
    a2_id: int
    order: int = 1
    style: BondStyle = BondStyle.PLAIN
    stereo: BondStereo = BondStereo.NONE
    is_aromatic: bool = False
    display_order: Optional[int] = None
    is_query: bool = False
    ring_id: Optional[int] = None
    length_px: Optional[float] = None


@dataclass
class ChemState:
    active_tool: str = "tool_atom"
    active_bond_order: int = 1
    active_bond_style: BondStyle = BondStyle.PLAIN
    active_bond_stereo: BondStereo = BondStereo.NONE
    active_bond_mode: str = "increment"
    active_bond_aromatic: bool = False
    fixed_angles: bool = True
    fixed_lengths: bool = True
    angle_step_deg: int = 30
    bond_length: float = 40.0
    active_ring_size: int = 6
    active_ring_aromatic: bool = False
    active_bracket_type: str = "[]"
    default_element: str = "C"
    selected_atoms: Set[int] = field(default_factory=set)
    selected_bonds: Set[int] = field(default_factory=set)
    # Display preferences
    show_implicit_carbons: bool = False  # Show C labels (False = hide like ChemDoodle)
    show_implicit_hydrogens: bool = False  # Show H labels
    use_aromatic_circles: bool = False  # Draw circles in aromatic rings
    label_font_family: str = "Arial"
    label_font_size: float = 11.0
    label_font_bold: bool = False
    label_font_italic: bool = False
    label_font_underline: bool = False
    use_element_colors: bool = True


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
        is_explicit: bool = False,
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
            is_explicit=is_explicit,
        )
        self.atoms[atom_id] = atom
        return atom

    def remove_atom(self, atom_id: int) -> tuple[Atom, List[Bond]]:
        atom = self.atoms.pop(atom_id)
        removed_bonds: List[Bond] = []
        for bond_id, bond in list(self.bonds.items()):
            if bond.a1_id == atom_id or bond.a2_id == atom_id:
                removed_bonds.append(self.remove_bond(bond_id))
        return atom, removed_bonds

    def add_bond(
        self,
        a1_id: int,
        a2_id: int,
        order: int = 1,
        bond_id: Optional[int] = None,
        style: BondStyle = BondStyle.PLAIN,
        stereo: BondStereo = BondStereo.NONE,
        is_aromatic: bool = False,
        display_order: Optional[int] = None,
        is_query: bool = False,
        ring_id: Optional[int] = None,
        length_px: Optional[float] = None,
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
            style=style,
            stereo=stereo,
            is_aromatic=is_aromatic,
            display_order=display_order,
            is_query=is_query,
            ring_id=ring_id,
            length_px=length_px,
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

    def update_atom_element(
        self,
        atom_id: int,
        element: str,
        is_explicit: Optional[bool] = None,
    ) -> None:
        atom = self.atoms[atom_id]
        atom.element = element
        if is_explicit is not None:
            atom.is_explicit = is_explicit

    def update_atom_charge(self, atom_id: int, charge: int) -> None:
        atom = self.atoms[atom_id]
        atom.charge = charge

    def update_bond(
        self,
        bond_id: int,
        order: Optional[int] = None,
        style: Optional[BondStyle] = None,
        stereo: Optional[BondStereo] = None,
        is_aromatic: Optional[bool] = None,
        display_order: Optional[int] = None,
    ) -> Bond:
        bond = self.bonds[bond_id]
        if order is not None:
            bond.order = order
        if style is not None:
            bond.style = style
        if stereo is not None:
            bond.stereo = stereo
        if is_aromatic is not None:
            bond.is_aromatic = is_aromatic
        if display_order is not None:
            bond.display_order = display_order
        return bond

    def update_bond_length(self, bond_id: int, length_px: Optional[float]) -> None:
        bond = self.bonds[bond_id]
        bond.length_px = length_px

    def clear(self) -> None:
        self.atoms.clear()
        self.bonds.clear()
        self._next_atom_id = 1
        self._next_bond_id = 1

    def validate(self) -> List[int]:
        bond_order_sum: Dict[int, int] = {atom_id: 0 for atom_id in self.atoms}
        for bond in self.bonds.values():
            if bond.a1_id in bond_order_sum:
                bond_order_sum[bond.a1_id] += bond.order
            if bond.a2_id in bond_order_sum:
                bond_order_sum[bond.a2_id] += bond.order

        errors: List[int] = []
        for atom_id, atom in self.atoms.items():
            expected = VALENCE_MAP.get(atom.element)
            if expected is None:
                continue
            if bond_order_sum.get(atom_id, 0) > expected:
                errors.append(atom_id)
        return errors
