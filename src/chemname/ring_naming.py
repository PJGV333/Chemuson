from __future__ import annotations

from typing import Iterable, List, Sequence

from chemcalc.valence import implicit_h_count

from .errors import ChemNameNotSupported
from .locants import Sub, orientation_key
from .molview import MolView
from .substituents import ALKYL, HALO_MAP, alkyl_length_linear


def enumerate_ring_numberings(ring_order: Sequence[int]) -> List[List[int]]:
    n = len(ring_order)
    numberings: List[List[int]] = []
    for i in range(n):
        forward = [ring_order[(i + j) % n] for j in range(n)]
        backward = [ring_order[(i - j) % n] for j in range(n)]
        numberings.append(forward)
        numberings.append(backward)
    return numberings


def choose_ring_orientation(
    view: MolView,
    ring_order: Sequence[int],
    opts,
    allow_hydroxy: bool = False,
    allow_nitro: bool = False,
) -> List[int]:
    best = None
    best_key = None
    for numbering in enumerate_ring_numberings(ring_order):
        subs = ring_substituents(
            view,
            numbering,
            allow_hydroxy=allow_hydroxy,
            allow_nitro=allow_nitro,
        )
        locants = sorted(sub.locant for sub in subs)
        key = orientation_key(subs, opts, primary_locants=locants)
        if best_key is None or key < best_key:
            best_key = key
            best = numbering
    return list(best) if best is not None else list(ring_order)


def ring_substituents(
    view: MolView,
    ring_order: Sequence[int],
    allow_hydroxy: bool = False,
    allow_nitro: bool = False,
) -> List[Sub]:
    ring_set = set(ring_order)
    index_map = {atom_id: idx + 1 for idx, atom_id in enumerate(ring_order)}
    subs: List[Sub] = []
    for atom_id in ring_order:
        locant = index_map[atom_id]
        for nbr in view.neighbors(atom_id):
            if nbr in ring_set:
                continue
            if view.element(nbr) == "H":
                continue
            name = _substituent_name_for_neighbor(
                view,
                atom_id,
                nbr,
                ring_set,
                allow_hydroxy=allow_hydroxy,
                allow_nitro=allow_nitro,
            )
            subs.append(Sub(name, locant))
    return subs


def _substituent_name_for_neighbor(
    view: MolView,
    ring_atom: int,
    nbr: int,
    ring_set: Iterable[int],
    allow_hydroxy: bool,
    allow_nitro: bool,
) -> str:
    elem = view.element(nbr)
    if elem in HALO_MAP:
        return HALO_MAP[elem]
    if elem == "C":
        length = alkyl_length_linear(view, nbr, set(ring_set))
        name = ALKYL.get(length)
        if name is None:
            raise ChemNameNotSupported("Unsupported alkyl length")
        return name
    if allow_hydroxy and elem == "O":
        if view.bond_order_between(ring_atom, nbr) != 1:
            raise ChemNameNotSupported("Unsupported hydroxy bond")
        heavy_neighbors = [n for n in view.neighbors(nbr) if view.element(n) != "H"]
        if len(heavy_neighbors) != 1:
            raise ChemNameNotSupported("Unsupported hydroxy substituent")
        if implicit_h_count(view, nbr) + view.explicit_h(nbr) < 1:
            raise ChemNameNotSupported("Unsupported hydroxy substituent")
        return "hydroxy"
    if allow_nitro and elem == "N":
        if implicit_h_count(view, nbr) + view.explicit_h(nbr) != 0:
            raise ChemNameNotSupported("Unsupported nitro substituent")
        heavy_neighbors = [n for n in view.neighbors(nbr) if view.element(n) != "H"]
        if len(heavy_neighbors) != 3:
            raise ChemNameNotSupported("Unsupported nitro substituent")
        oxygen_neighbors = [n for n in heavy_neighbors if view.element(n) == "O"]
        if len(oxygen_neighbors) != 2:
            raise ChemNameNotSupported("Unsupported nitro substituent")
        for o_id in oxygen_neighbors:
            order = view.bond_order_between(nbr, o_id)
            if order not in {1, 2}:
                raise ChemNameNotSupported("Unsupported nitro substituent")
        return "nitro"

    raise ChemNameNotSupported("Unsupported substituent")
