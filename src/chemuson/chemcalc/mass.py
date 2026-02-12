"""Cálculo de masas moleculares a partir de fórmulas."""

from __future__ import annotations

from typing import Dict

# Pesos atómicos promedio (u) para cálculo aproximado de masa molecular.
ATOMIC_WEIGHTS: Dict[str, float] = {
    "H": 1.00794,
    "C": 12.0107,
    "N": 14.0067,
    "O": 15.9994,
    "F": 18.998403,
    "Cl": 35.453,
    "Br": 79.904,
    "I": 126.90447,
}


def molecular_weight(formula_dict: Dict[str, int]) -> float:
    """Calcula el peso molecular a partir de una fórmula.

    Args:
        formula_dict: Diccionario de elemento -> conteo.

    Returns:
        Masa molecular aproximada en unidades atómicas (u).

    Raises:
        ValueError: Si el peso atómico de un elemento no está disponible.

    Side Effects:
        No tiene efectos laterales.
    """
    total = 0.0
    for element, count in formula_dict.items():
        weight = ATOMIC_WEIGHTS.get(element)
        if weight is None:
            raise ValueError(f"Atomic weight not available for {element}")
        total += weight * count
    return total
