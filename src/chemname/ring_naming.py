from __future__ import annotations

from typing import Iterable, List, Sequence

from chemcalc.valence import implicit_h_count

from .errors import ChemNameNotSupported
from .locants import Sub, orientation_key
from .molview import MolView
from .substituents import (
    HALO_MAP,
    alkoxy_substituent_name,
    amino_substituent_name,
    alkyl_substituent_name,
    halomethyl_substituent_name,
    nitro_substituent_name,
    ring_substituent_name,
)
from .rings import RingContext


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
    allow_amino: bool = False,
    allow_alkoxy: bool = False,
    forbid_hetero_substituents: bool = False,
    ring_ctx: RingContext | None = None,
) -> List[int]:
    best = None
    best_key = None
    for numbering in enumerate_ring_numberings(ring_order):
        subs = ring_substituents(
            view,
            numbering,
            allow_hydroxy=allow_hydroxy,
            allow_nitro=allow_nitro,
            allow_amino=allow_amino,
            allow_alkoxy=allow_alkoxy,
            forbid_hetero_substituents=forbid_hetero_substituents,
            ring_ctx=ring_ctx,
        )
        locants = sorted(sub.locant for sub in subs)
        key = orientation_key(subs, opts, primary_locants=locants)
        if best_key is None or key < best_key:
            best_key = key
            best = numbering
    return list(best) if best is not None else list(ring_order)


def choose_hetero_ring_orientation(
    view: MolView,
    ring_order: Sequence[int],
    hetero_atoms: Iterable[int],
    opts,
    allow_hydroxy: bool = False,
    allow_nitro: bool = False,
    allow_amino: bool = False,
    allow_alkoxy: bool = False,
    preferred_start_atoms: Iterable[int] | None = None,
    forbid_hetero_substituents: bool = False,
    ring_ctx: RingContext | None = None,
) -> List[int]:
    hetero_set = set(hetero_atoms)
    preferred_set = set(preferred_start_atoms or [])
    best = None
    best_key = None
    for numbering in enumerate_ring_numberings(ring_order):
        if preferred_set:
            if numbering[0] not in preferred_set:
                continue
        else:
            if numbering[0] not in hetero_set:
                continue
        subs = ring_substituents(
            view,
            numbering,
            allow_hydroxy=allow_hydroxy,
            allow_nitro=allow_nitro,
            allow_amino=allow_amino,
            allow_alkoxy=allow_alkoxy,
            forbid_hetero_substituents=forbid_hetero_substituents,
            ring_ctx=ring_ctx,
        )
        hetero_locants = [idx + 1 for idx, atom_id in enumerate(numbering) if atom_id in hetero_set]
        key = orientation_key(subs, opts, primary_locants=hetero_locants)
        if best_key is None or key < best_key:
            best_key = key
            best = numbering

    if best is None and preferred_set:
        return choose_hetero_ring_orientation(
            view,
            ring_order,
            hetero_atoms,
            opts,
            allow_hydroxy=allow_hydroxy,
            allow_nitro=allow_nitro,
            allow_amino=allow_amino,
            allow_alkoxy=allow_alkoxy,
            preferred_start_atoms=None,
            forbid_hetero_substituents=forbid_hetero_substituents,
            ring_ctx=ring_ctx,
        )
    return list(best) if best is not None else list(ring_order)


def ring_substituents(
    view: MolView,
    ring_order: Sequence[int],
    allow_hydroxy: bool = False,
    allow_nitro: bool = False,
    allow_amino: bool = False,
    allow_alkoxy: bool = False,
    forbid_hetero_substituents: bool = False,
    ring_ctx: RingContext | None = None,
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
            if forbid_hetero_substituents and view.element(atom_id) != "C":
                raise ChemNameNotSupported("Unsupported hetero substituent")
            name = _substituent_name_for_neighbor(
                view,
                atom_id,
                nbr,
                ring_set,
                allow_hydroxy=allow_hydroxy,
                allow_nitro=allow_nitro,
                allow_amino=allow_amino,
                allow_alkoxy=allow_alkoxy,
                ring_ctx=ring_ctx,
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
    allow_amino: bool,
    allow_alkoxy: bool,
    ring_ctx: RingContext | None = None,
) -> str:
    elem = view.element(nbr)
    if elem in HALO_MAP:
        return HALO_MAP[elem]
    if elem == "C":
        halo_name = halomethyl_substituent_name(view, nbr, set(ring_set))
        if halo_name is not None:
            return halo_name
        ring_name = ring_substituent_name(view, nbr, set(ring_set), ring_ctx)
        if ring_name is not None:
            return ring_name
        return alkyl_substituent_name(view, nbr, set(ring_set))
    if elem == "O":
        if allow_alkoxy:
            name = alkoxy_substituent_name(view, nbr, set(ring_set))
            if name is not None:
                return name
        if not allow_hydroxy:
            raise ChemNameNotSupported("Unsupported substituent")
        if view.bond_order_between(ring_atom, nbr) != 1:
            raise ChemNameNotSupported("Unsupported hydroxy bond")
        heavy_neighbors = [n for n in view.neighbors(nbr) if view.element(n) != "H"]
        if len(heavy_neighbors) != 1:
            raise ChemNameNotSupported("Unsupported hydroxy substituent")
        if implicit_h_count(view, nbr) + view.explicit_h(nbr) < 1:
            raise ChemNameNotSupported("Unsupported hydroxy substituent")
        return "hydroxy"
    if elem == "N":
        if allow_nitro:
            name = nitro_substituent_name(view, nbr, set(ring_set))
            if name is not None:
                return name
        if allow_amino:
            name = amino_substituent_name(view, nbr, set(ring_set))
            if name is not None:
                return name

    raise ChemNameNotSupported("Unsupported substituent")
