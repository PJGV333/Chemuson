from __future__ import annotations

from typing import Dict

from chemname.molview import MolView
from .valence import implicit_h_count


def molecular_formula(graph) -> Dict[str, int]:
    """Return a molecular formula as a dict of element -> count."""
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
    """Format a formula dict using Hill order (C, H, then alphabetical)."""
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
