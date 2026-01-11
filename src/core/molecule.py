from typing import List, Tuple

class Atom:
    """
    Representa un átomo en la molécula.
    """
    def __init__(self, symbol: str, x: float, y: float) -> None:
        self.symbol = symbol
        self.x = x
        self.y = y

class Bond:
    """
    Representa un enlace entre dos átomos.
    """
    def __init__(self, atom1_idx: int, atom2_idx: int, order: int = 1) -> None:
        self.atom1_idx = atom1_idx
        self.atom2_idx = atom2_idx
        self.order = order

class Molecule:
    """
    Clase principal para manejar la estructura química.
    """
    def __init__(self) -> None:
        self.atoms: List[Atom] = []
        self.bonds: List[Bond] = []

    def add_atom(self, symbol: str, x: float, y: float) -> int:
        """Añade un átomo y retorna su índice."""
        atom = Atom(symbol, x, y)
        self.atoms.append(atom)
        return len(self.atoms) - 1

    def add_bond(self, atom1_idx: int, atom2_idx: int, order: int = 1) -> None:
        """Añade un enlace entre dos átomos."""
        bond = Bond(atom1_idx, atom2_idx, order)
        self.bonds.append(bond)

    def clear(self) -> None:
        """Limpia la molécula."""
        self.atoms = []
        self.bonds = []
