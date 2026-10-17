"""Renderer-independent geometry contracts for special molecular projections.

Projection geometry describes coordinates and projection-specific bond semantics;
global stroke widths, caps, joins, fonts, and colors remain DrawingStyle concerns.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

from chemuson.core.model import BondStyle, MolGraph


@dataclass(frozen=True, slots=True)
class ProjectionAtom:
    """One atom/label at a position in a two-dimensional projection."""

    element: str
    x: float
    y: float
    is_explicit: bool


@dataclass(frozen=True, slots=True)
class ProjectionBond:
    """Connectivity and projection-local emphasis between geometry atoms."""

    start: int
    end: int
    order: int = 1
    style: BondStyle = BondStyle.PLAIN


@dataclass(frozen=True, slots=True)
class ProjectionGeometry:
    """Immutable, deterministic geometry ready to adapt to the molecular model."""

    family: str
    bond_length: float
    atoms: tuple[ProjectionAtom, ...]
    bonds: tuple[ProjectionBond, ...]

    def __post_init__(self) -> None:
        if not self.family.strip():
            raise ValueError("projection family must not be empty")
        if not math.isfinite(self.bond_length) or self.bond_length <= 0.0:
            raise ValueError("projection bond_length must be finite and positive")
        if not self.atoms:
            raise ValueError("projection geometry must contain at least one atom")
        for atom in self.atoms:
            if not atom.element.strip():
                raise ValueError("projection atom labels must not be empty")
            if not math.isfinite(atom.x) or not math.isfinite(atom.y):
                raise ValueError("projection coordinates must be finite")
        for bond in self.bonds:
            if not (0 <= bond.start < len(self.atoms)):
                raise ValueError(f"projection bond start index out of range: {bond.start}")
            if not (0 <= bond.end < len(self.atoms)):
                raise ValueError(f"projection bond end index out of range: {bond.end}")
            if bond.start == bond.end:
                raise ValueError("projection bond endpoints must differ")
            if (
                isinstance(bond.order, bool)
                or not isinstance(bond.order, int)
                or bond.order < 1
            ):
                raise ValueError("projection bond order must be a positive integer")
            if not isinstance(bond.style, BondStyle):
                raise ValueError("projection bond style must be a BondStyle")

    def to_molgraph(self) -> MolGraph:
        """Create a fresh molecular graph without adding stereo annotations."""
        graph = MolGraph()
        atom_ids: list[int] = []
        for atom in self.atoms:
            created = graph.add_atom(
                atom.element,
                atom.x,
                atom.y,
                is_explicit=atom.is_explicit,
            )
            atom_ids.append(created.id)
        for bond in self.bonds:
            graph.add_bond(
                atom_ids[bond.start],
                atom_ids[bond.end],
                order=bond.order,
                style=bond.style,
            )
        return graph


def projection_atom(element: str, x: float, y: float, *, explicit: bool | None = None) -> ProjectionAtom:
    """Build a projection atom using the editor's implicit-carbon convention."""
    if explicit is None:
        explicit = element != "C"
    return ProjectionAtom(element, float(x), float(y), bool(explicit))


def projection_bond(
    start: int,
    end: int,
    *,
    order: int = 1,
    style: BondStyle = BondStyle.PLAIN,
) -> ProjectionBond:
    """Build a projection bond with validated local drawing semantics."""
    return ProjectionBond(start, end, order, style)
