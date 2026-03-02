"""Cálculo de hidrógenos implícitos según valencias típicas."""

from __future__ import annotations

from typing import Dict

from chemuson.chemname.molview import MolView

# Valencias típicas usadas para inferir H implícitos en cálculos sencillos.
TYPICAL_VALENCE: Dict[str, int] = {
    "H": 1,
    "C": 4,
    "N": 3,
    "O": 2,
    "F": 1,
    "Cl": 1,
    "Br": 1,
    "I": 1,
    "S": 2,
    "P": 3,
}


def implicit_h_count(view: MolView, atom_id: int) -> int:
    """Calcula los hidrógenos implícitos para un átomo.

    Args:
        view: Vista del grafo molecular con operaciones de consulta.
        atom_id: Identificador del átomo a evaluar.

    Returns:
        Número de H implícitos estimados (>= 0).

    Side Effects:
        No tiene efectos laterales.
    """
    graph = getattr(view, "graph", None)
    if graph is not None and hasattr(graph, "implicit_h_count"):
        try:
            return int(max(0, graph.implicit_h_count(atom_id)))
        except Exception:
            pass

    element = view.element(atom_id)
    typical = TYPICAL_VALENCE.get(element)
    if typical is None:
        return 0

    atom = None
    try:
        atom = view._get_atom(atom_id)  # type: ignore[attr-defined]
    except Exception:
        atom = None
    if atom is not None:
        no_implicit = False
        if isinstance(atom, dict):
            no_implicit = bool(atom.get("no_implicit", False))
        else:
            no_implicit = bool(getattr(atom, "no_implicit", False))
        if no_implicit:
            return 0

    bond_order_sum = 0
    for nbr in view.neighbors(atom_id):
        order = view.bond_order_between(atom_id, nbr)
        if view.bond_is_aromatic(atom_id, nbr):
            order = 1.5
        bond_order_sum += order

    charge = view.atom_charge(atom_id)
    # Ajustes simples de especies iónicas frecuentes para evitar falsos positivos.
    if element == "N" and charge > 0:
        typical = 4
    elif element == "O" and charge < 0:
        typical = 1
    explicit_h = view.explicit_h(atom_id)
    implicit = typical - bond_order_sum - explicit_h - max(charge, 0)
    if implicit < 0:
        return 0
    rounded = int(round(implicit))
    if rounded < 0:
        return 0
    return rounded
