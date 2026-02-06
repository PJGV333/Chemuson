"""Cálculo y formateo de fórmulas moleculares.

Este módulo agrega utilidades para contar elementos a partir de un grafo
molecular y formatear la fórmula siguiendo el orden de Hill.
"""

from __future__ import annotations

from typing import Dict

from chemname.molview import MolView
from .valence import implicit_h_count


def molecular_formula(graph) -> Dict[str, int]:
    """Calcula la fórmula molecular como diccionario de elemento -> conteo.

    Args:
        graph: Grafo molecular compatible con `MolView`.

    Returns:
        Diccionario con símbolos atómicos y sus cantidades totales.

    Side Effects:
        No tiene efectos laterales; solo calcula y devuelve datos.
    """
    view = MolView(graph)
    counts: Dict[str, int] = {}

    for atom_id in view.atoms():
        element = view.element(atom_id)
        counts[element] = counts.get(element, 0) + 1
        explicit_h = view.explicit_h(atom_id)
        if explicit_h:
            counts["H"] = counts.get("H", 0) + int(explicit_h)

    for atom_id in view.atoms():
        if view.element(atom_id) == "H":
            continue
        implicit = implicit_h_count(view, atom_id)
        if implicit:
            counts["H"] = counts.get("H", 0) + int(implicit)

    return {element: count for element, count in counts.items() if count > 0}


def format_formula(formula_dict: Dict[str, int]) -> str:
    """Formatea una fórmula usando el orden de Hill (C, H, luego alfabético).

    Args:
        formula_dict: Diccionario con símbolos de elementos y cantidades.

    Returns:
        Cadena con la fórmula formateada (p. ej., "C6H6O").

    Side Effects:
        No tiene efectos laterales.
    """
    if not formula_dict:
        return ""
    order = []
    if "C" in formula_dict:
        order.append("C")
    if "H" in formula_dict:
        order.append("H")
    for element in sorted(e for e in formula_dict.keys() if e not in {"C", "H"}):
        order.append(element)

    parts = []
    for element in order:
        count = formula_dict.get(element, 0)
        if count <= 0:
            continue
        parts.append(element if count == 1 else f"{element}{count}")
    return "".join(parts)
