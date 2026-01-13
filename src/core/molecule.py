from core.model import Atom, Bond, MolGraph


class Molecule(MolGraph):
    """
    Deprecated compatibility shim.
    Prefer core.model.MolGraph for new code.
    """

    def add_atom(self, symbol: str, x: float, y: float, is_explicit: bool = False) -> int:
        atom = super().add_atom(symbol, x, y, is_explicit=is_explicit)
        return atom.id

    def add_bond(self, atom1_idx: int, atom2_idx: int, order: int = 1) -> None:
        super().add_bond(atom1_idx, atom2_idx, order)


__all__ = ["Atom", "Bond", "Molecule"]
