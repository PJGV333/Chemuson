"""Cálculo de hidrógenos implícitos según valencias típicas."""

from __future__ import annotations

from typing import Dict

from chemname.molview import MolView

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
    element = view.element(atom_id)
    typical = TYPICAL_VALENCE.get(element)
    if typical is None:
        return 0

    bond_order_sum = 0
    for nbr in view.neighbors(atom_id):
        bond_order_sum += view.bond_order_between(atom_id, nbr)

    charge = view.atom_charge(atom_id)
    implicit = typical - bond_order_sum - charge
    if implicit < 0:
        return 0
    return int(implicit)
