"""Detectores de grupos funcionales complementarios para IUPAC-lite.

Este módulo centraliza patrones simples usados por `locants`, `ring_naming`
y el motor principal para mantener la lógica modular.
"""

from __future__ import annotations

from chemuson.chemcalc.valence import implicit_h_count

from .molview import MolView

_ALKYL_PREFIX = {
    1: "methyl",
    2: "ethyl",
    3: "propyl",
    4: "butyl",
    5: "pentyl",
    6: "hexyl",
}


def isotope_substituent_name(view: MolView, atom_id: int) -> str | None:
    """Devuelve prefijo isotópico para átomos sustituyentes simples."""
    isotope = view.isotope(atom_id)
    if isotope is None:
        return None
    elem = view.element(atom_id)
    if elem == "H":
        if isotope == 2:
            return "deuterio"
        if isotope == 3:
            return "tritio"
    return f"({isotope}{elem})"


def radical_substituent_name(
    view: MolView,
    start_atom: int,
    parent_set: set[int],
) -> str | None:
    """Detecta sustituyente radical simple terminal (oxyl)."""
    if not view.has_radical(start_atom):
        return None
    if view.element(start_atom) != "O":
        return None
    heavy_neighbors = [nbr for nbr in view.neighbors(start_atom) if view.element(nbr) != "H"]
    parent_neighbors = [nbr for nbr in heavy_neighbors if nbr in parent_set]
    if len(parent_neighbors) != 1:
        return None
    if len(heavy_neighbors) != 1:
        return None
    return "oxyl"


def thiol_substituent_name(
    view: MolView, start_atom: int, parent_set: set[int]
) -> str | None:
    """Detecta tiol terminal como prefijo mercapto."""
    if view.element(start_atom) != "S":
        return None
    heavy_neighbors = [nbr for nbr in view.neighbors(start_atom) if view.element(nbr) != "H"]
    parent_neighbors = [nbr for nbr in heavy_neighbors if nbr in parent_set]
    if len(parent_neighbors) != 1 or len(heavy_neighbors) != 1:
        return None
    if view.bond_order_between(start_atom, parent_neighbors[0]) != 1:
        return None
    h_total = implicit_h_count(view, start_atom) + view.explicit_h(start_atom)
    if h_total < 1:
        return None
    return "mercapto"


def azido_substituent_name(
    view: MolView, start_atom: int, parent_set: set[int]
) -> str | None:
    """Detecta sustituyente azido (-N3)."""
    if view.element(start_atom) != "N":
        return None
    heavy_neighbors = [nbr for nbr in view.neighbors(start_atom) if view.element(nbr) != "H"]
    parent_neighbors = [nbr for nbr in heavy_neighbors if nbr in parent_set]
    if len(parent_neighbors) != 1:
        return None
    n2_candidates = [nbr for nbr in heavy_neighbors if nbr not in parent_set and view.element(nbr) == "N"]
    if len(n2_candidates) != 1:
        return None
    n2 = n2_candidates[0]
    n2_heavy = [nbr for nbr in view.neighbors(n2) if view.element(nbr) != "H"]
    if start_atom not in n2_heavy:
        return None
    tail = [nbr for nbr in n2_heavy if nbr != start_atom and view.element(nbr) == "N"]
    if len(tail) != 1:
        return None
    n3 = tail[0]
    n3_heavy = [nbr for nbr in view.neighbors(n3) if view.element(nbr) != "H"]
    if n3_heavy != [n2]:
        return None
    order_12 = view.bond_order_between(start_atom, n2)
    order_23 = view.bond_order_between(n2, n3)
    if order_12 not in {1, 2} or order_23 not in {1, 2}:
        return None
    return "azido"


def peroxy_substituent_name(
    view: MolView, start_atom: int, parent_set: set[int]
) -> str | None:
    """Detecta sustituyente peroxy simple (-O-O-)."""
    if view.element(start_atom) != "O":
        return None
    heavy_neighbors = [nbr for nbr in view.neighbors(start_atom) if view.element(nbr) != "H"]
    parent_neighbors = [nbr for nbr in heavy_neighbors if nbr in parent_set]
    if len(parent_neighbors) != 1:
        return None
    oxygen_neighbors = [
        nbr
        for nbr in heavy_neighbors
        if nbr not in parent_set and view.element(nbr) == "O"
    ]
    if len(oxygen_neighbors) != 1:
        return None
    if view.bond_order_between(start_atom, oxygen_neighbors[0]) != 1:
        return None
    return "peroxy"


def sulfinyl_substituent_name(
    view: MolView, start_atom: int, parent_set: set[int]
) -> str | None:
    """Detecta sustituyentes sulfinilo tipo R-S(=O)-."""
    if view.element(start_atom) != "S":
        return None
    heavy_neighbors = [nbr for nbr in view.neighbors(start_atom) if view.element(nbr) != "H"]
    parent_neighbors = [nbr for nbr in heavy_neighbors if nbr in parent_set]
    if len(parent_neighbors) != 1:
        return None
    o_double = [
        nbr
        for nbr in heavy_neighbors
        if nbr not in parent_set
        and view.element(nbr) == "O"
        and view.bond_order_between(start_atom, nbr) == 2
    ]
    if len(o_double) != 1:
        return None
    rest = [nbr for nbr in heavy_neighbors if nbr not in set(parent_neighbors) | set(o_double)]
    if len(rest) == 1 and view.element(rest[0]) == "C":
        alkyl = _simple_alkyl_prefix(view, rest[0], parent_set | {start_atom} | set(o_double))
        if alkyl:
            return f"{alkyl}sulfinyl"
    return "sulfinyl"


def sulfonyl_substituent_name(
    view: MolView, start_atom: int, parent_set: set[int]
) -> str | None:
    """Detecta sustituyentes sulfonilo tipo R-S(=O)2-."""
    if view.element(start_atom) != "S":
        return None
    heavy_neighbors = [nbr for nbr in view.neighbors(start_atom) if view.element(nbr) != "H"]
    parent_neighbors = [nbr for nbr in heavy_neighbors if nbr in parent_set]
    if len(parent_neighbors) != 1:
        return None
    o_double = [
        nbr
        for nbr in heavy_neighbors
        if nbr not in parent_set
        and view.element(nbr) == "O"
        and view.bond_order_between(start_atom, nbr) == 2
    ]
    if len(o_double) != 2:
        return None
    rest = [nbr for nbr in heavy_neighbors if nbr not in set(parent_neighbors) | set(o_double)]
    if len(rest) == 1 and view.element(rest[0]) == "C":
        alkyl = _simple_alkyl_prefix(view, rest[0], parent_set | {start_atom} | set(o_double))
        if alkyl:
            return f"{alkyl}sulfonyl"
    return "sulfonyl"


def sulfonamido_substituent_name(
    view: MolView, start_atom: int, parent_set: set[int]
) -> str | None:
    """Detecta sustituyentes sulfonamido (-SO2NR-)."""
    if view.element(start_atom) != "S":
        return None
    heavy_neighbors = [nbr for nbr in view.neighbors(start_atom) if view.element(nbr) != "H"]
    parent_neighbors = [nbr for nbr in heavy_neighbors if nbr in parent_set]
    if len(parent_neighbors) != 1:
        return None
    o_double = [
        nbr
        for nbr in heavy_neighbors
        if view.element(nbr) == "O" and view.bond_order_between(start_atom, nbr) == 2
    ]
    n_single = [
        nbr
        for nbr in heavy_neighbors
        if view.element(nbr) == "N" and view.bond_order_between(start_atom, nbr) == 1
    ]
    if len(o_double) == 2 and len(n_single) == 1:
        return "sulfonamido"
    return None


def detect_sulfonic_attachment(
    view: MolView, sulfur_atom: int, parent_set: set[int]
) -> dict | None:
    """Detecta patrón de ácido sulfónico/sulfonato unido al esqueleto."""
    if view.element(sulfur_atom) != "S":
        return None
    heavy_neighbors = [nbr for nbr in view.neighbors(sulfur_atom) if view.element(nbr) != "H"]
    parent_neighbors = [nbr for nbr in heavy_neighbors if nbr in parent_set]
    if len(parent_neighbors) != 1:
        return None
    o_double = [
        nbr
        for nbr in heavy_neighbors
        if view.element(nbr) == "O" and view.bond_order_between(sulfur_atom, nbr) == 2
    ]
    o_single = [
        nbr
        for nbr in heavy_neighbors
        if view.element(nbr) == "O" and view.bond_order_between(sulfur_atom, nbr) == 1
    ]
    if len(o_double) != 2 or len(o_single) != 1:
        return None
    oxy = o_single[0]
    h_total = implicit_h_count(view, oxy) + view.explicit_h(oxy)
    oxy_charge = view.formal_charge(oxy)
    aux = set(o_double) | {oxy}
    if h_total >= 1:
        return {
            "kind": "sulfonic_acid",
            "parent_atom": parent_neighbors[0],
            "aux_atoms": aux,
        }
    if oxy_charge <= -1:
        return {
            "kind": "sulfonate",
            "parent_atom": parent_neighbors[0],
            "aux_atoms": aux,
        }
    return None


def _simple_alkyl_prefix(
    view: MolView,
    start_atom: int,
    blocked_atoms: set[int],
    max_len: int = 6,
) -> str | None:
    """Obtiene prefijo alquilo lineal para sustituyentes sobre heteroátomos."""
    if view.element(start_atom) != "C":
        return None
    prev = None
    current = start_atom
    length = 0
    visited: set[int] = set()
    while True:
        if current in visited:
            return None
        visited.add(current)
        if view.element(current) != "C":
            return None
        length += 1
        if length > max_len:
            return None
        next_candidates = []
        for nbr in view.neighbors(current):
            if nbr in blocked_atoms or nbr == prev or view.element(nbr) == "H":
                continue
            if view.element(nbr) != "C":
                return None
            if view.bond_order_between(current, nbr) != 1:
                return None
            next_candidates.append(nbr)
        if len(next_candidates) > 1:
            return None
        if not next_candidates:
            break
        prev, current = current, next_candidates[0]
    return _ALKYL_PREFIX.get(length)
