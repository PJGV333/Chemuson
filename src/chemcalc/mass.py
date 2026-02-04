from __future__ import annotations

from typing import Dict

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
    """Return molecular weight from a formula dict."""
    total = 0.0
    for element, count in formula_dict.items():
        weight = ATOMIC_WEIGHTS.get(element)
        if weight is None:
            raise ValueError(f"Atomic weight not available for {element}")
        total += weight * count
    return total
