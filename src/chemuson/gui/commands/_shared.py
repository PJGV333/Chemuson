"""Comandos de edición con soporte de deshacer/rehacer (QUndoCommand).

Este módulo centraliza operaciones de modificación del modelo/GUI
para integrarlas con la pila de undo/redo de Qt.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Dict, Iterable, List, Optional, Tuple

from PyQt6.QtCore import QPointF, QRectF
from PyQt6.QtGui import QFont, QUndoCommand, QColor
from PyQt6.QtWidgets import QGraphicsTextItem

from chemuson.core.model import BondStyle, BondStereo, MolGraph, bond_is_structural
from chemuson.gui.diagram_models import SemanticDiagram
from chemuson.gui.geom import angle_deg, angle_distance_deg, endpoint_from_angle_len

_IMPLICIT_ELEMENTS = {"C"}
_ANCHOR_UNSET = object()
_DONOR_UNSET = object()
_FLEX_CURVE_UNSET = object()
_SPHERE_STYLE_UNSET = object()
_BOND_LENGTH_UNSET = object()
_BOND_STROKE_UNSET = object()
_LABEL_SCALE_UNSET = object()
_OPACITY_UNSET = object()


def _default_is_explicit(element: str) -> bool:
    """Determina si un elemento debe mostrarse explícitamente por defecto."""
    return element not in _IMPLICIT_ELEMENTS


def _atom_degree(model: MolGraph, atom_id: int) -> int:
    """Calcula el grado estructural (sin enlaces intermoleculares) de un átomo."""
    degree_getter = getattr(model, "atom_degree", None)
    if callable(degree_getter):
        try:
            return int(degree_getter(atom_id))
        except Exception:
            pass
    return sum(
        1
        for bond in model.bonds.values()
        if bond_is_structural(bond) and (bond.a1_id == atom_id or bond.a2_id == atom_id)
    )


def _is_hidden_carbon_placeholder(model: MolGraph, atom_id: int) -> bool:
    """Indica si el átomo es un carbono implícito auto-generado sin semántica propia."""
    checker = getattr(model, "is_hidden_carbon_placeholder", None)
    if callable(checker):
        return bool(checker(atom_id))
    atom = model.atoms.get(atom_id)
    if atom is None:
        return False
    return atom.element == "C" and not atom.is_explicit


def _neighbor_angles_deg(model: MolGraph, atom_id: int) -> list[float]:
    """Devuelve los ángulos hacia los vecinos de un átomo."""
    anchor = model.get_atom(atom_id)
    origin = QPointF(anchor.x, anchor.y)
    angles: list[float] = []
    for bond in model.bonds.values():
        if not bond_is_structural(bond):
            continue
        if bond.a1_id == atom_id:
            other = model.get_atom(bond.a2_id)
        elif bond.a2_id == atom_id:
            other = model.get_atom(bond.a1_id)
        else:
            continue
        angles.append(angle_deg(origin, QPointF(other.x, other.y)))
    return angles


def _select_hydrogen_angles(existing_angles_deg: list[float], count: int) -> list[float]:
    """Selecciona ángulos base para colocar hidrógenos implícitos."""
    base_angles = [0.0, 120.0, 240.0]
    if count <= 0:
        return []
    if not existing_angles_deg:
        return base_angles[:count]

    def min_distance(angle: float) -> float:
        """Calcula la distancia mínima a los ángulos existentes."""
        return min(angle_distance_deg(angle, existing) for existing in existing_angles_deg)

    ordered = sorted(base_angles, key=min_distance, reverse=True)
    return ordered[:count]


def _bond_length_from_view(view) -> float:
    """Obtiene la longitud de enlace configurada desde la vista."""
    return getattr(getattr(view, "state", None), "bond_length", 40.0)


def _resolve_atom_label_spec(view, label: str) -> dict:
    """Normaliza etiquetas de UI a la especificación química persistente."""
    resolver = getattr(view, "resolve_atom_label_spec", None)
    if callable(resolver):
        try:
            spec = resolver(label)
            if isinstance(spec, dict):
                return spec
        except Exception:
            pass
    return {
        "element": label,
        "group_h_cap": None,
        "explicit_h": None,
        "no_implicit": False,
    }


def _validate_view_structure(view) -> None:
    """Dispara la validación visual si la vista la soporta."""
    requester = getattr(view, "request_structure_validation", None)
    if callable(requester):
        requester()
        return
    validator = getattr(view, "validate_structure", None)
    if callable(validator):
        validator()


def _remove_hydrogen_specs(model: MolGraph, view, specs: list[tuple[int, float, float, int]]) -> None:
    """Elimina hidrógenos y enlaces descritos por `specs` del modelo y la vista."""
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
    """Recrea hidrógenos y enlaces a partir de `specs`."""
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
    """Crea hidrógenos explícitos alrededor de un átomo ancla.

    Returns:
        Lista de especificaciones `(atom_id, x, y, bond_id)` creadas.
    """
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
    """Recolecta hidrógenos unidos a un átomo para poder restaurarlos."""
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


__all__ = [name for name in globals() if not name.startswith("__")]

