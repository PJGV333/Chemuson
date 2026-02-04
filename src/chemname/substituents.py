from __future__ import annotations

from typing import Dict, Set

from .errors import ChemNameNotSupported
from .molview import MolView

HALO_MAP: Dict[str, str] = {
    "F": "fluoro",
    "Cl": "chloro",
    "Br": "bromo",
    "I": "iodo",
}

ALKANE_PARENT: Dict[int, str] = {
    1: "methane",
    2: "ethane",
    3: "propane",
    4: "butane",
    5: "pentane",
    6: "hexane",
    7: "heptane",
    8: "octane",
    9: "nonane",
    10: "decane",
    11: "undecane",
    12: "dodecane",
    13: "tridecane",
    14: "tetradecane",
    15: "pentadecane",
    16: "hexadecane",
    17: "heptadecane",
    18: "octadecane",
    19: "nonadecane",
    20: "eicosane",
}

ALKYL: Dict[int, str] = {
    1: "methyl",
    2: "ethyl",
    3: "propyl",
    4: "butyl",
    5: "pentyl",
    6: "hexyl",
    7: "heptyl",
    8: "octyl",
    9: "nonyl",
    10: "decyl",
    11: "undecyl",
    12: "dodecyl",
    13: "tridecyl",
    14: "tetradecyl",
    15: "pentadecyl",
    16: "hexadecyl",
    17: "heptadecyl",
    18: "octadecyl",
    19: "nonadecyl",
    20: "eicosyl",
}


def alkane_root(parent: str) -> str:
    if parent.endswith("ane"):
        return parent[:-3]
    return parent


def parent_name(
    length: int,
    unsat_order: int | None = None,
    unsat_locant: int | None = None,
    suffix: str | None = None,
    suffix_locant: int | None = None,
) -> str:
    parent = ALKANE_PARENT.get(length)
    if parent is None:
        raise ChemNameNotSupported("Unsupported parent length")

    if suffix is not None:
        if suffix_locant is None:
            raise ChemNameNotSupported("Missing suffix locant")
        if unsat_order is not None:
            if unsat_locant is None:
                raise ChemNameNotSupported("Missing unsaturation locant")
            root = alkane_root(parent)
            infix = "en" if unsat_order == 2 else "yn"
            base = f"{root}-{unsat_locant}-{infix}"
        else:
            base = parent[:-1] if parent.endswith("e") else parent
        return f"{base}-{suffix_locant}-{suffix}"

    if unsat_order is not None:
        if unsat_locant is None:
            raise ChemNameNotSupported("Missing unsaturation locant")
        root = alkane_root(parent)
        infix = "ene" if unsat_order == 2 else "yne"
        return f"{root}-{unsat_locant}-{infix}"

    return parent


def alkyl_length_linear(view: MolView, start_atom: int, chain_set: Set[int]) -> int:
    """Return the length of a linear alkyl substituent.

    Raises ChemNameNotSupported if the branch is not a simple linear alkyl.
    """
    if start_atom in chain_set:
        raise ChemNameNotSupported("Branch atom is part of the chain")

    visited: Set[int] = set()
    current = start_atom
    prev = None
    length = 0

    while True:
        if current in visited:
            raise ChemNameNotSupported("Branch cycle detected")
        visited.add(current)

        if view.element(current) != "C":
            raise ChemNameNotSupported("Non-carbon in alkyl branch")

        length += 1
        neighbors = [nbr for nbr in view.neighbors(current) if nbr not in chain_set]
        next_candidates = []
        for nbr in neighbors:
            if nbr == prev:
                continue
            elem = view.element(nbr)
            if elem == "H":
                continue
            if elem != "C":
                raise ChemNameNotSupported("Non-carbon in alkyl branch")
            next_candidates.append(nbr)

        if len(next_candidates) > 1:
            raise ChemNameNotSupported("Alkyl branch is branched")
        if not next_candidates:
            break

        prev, current = current, next_candidates[0]

    return length
