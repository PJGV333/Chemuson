from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence, Set, Tuple

from .errors import ChemNameNotSupported
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


@dataclass(frozen=True)
class Sub:
    name: str
    locant: int


def substituents_on_chain(
    view: MolView,
    chain: Sequence[int],
    ignore_atoms: Optional[Set[int]] = None,
    ring_ctx: RingContext | None = None,
) -> List[Sub]:
    """Return substituents attached to the given chain."""
    chain_set = set(chain)
    ignored = ignore_atoms or set()
    index_map = {atom_id: idx + 1 for idx, atom_id in enumerate(chain)}
    substituents: List[Sub] = []

    for atom_id in chain:
        locant = index_map[atom_id]
        for nbr in view.neighbors(atom_id):
            if nbr in ignored:
                continue
            if nbr in chain_set:
                continue
            elem = view.element(nbr)
            if elem == "H":
                continue
            if elem in HALO_MAP:
                substituents.append(Sub(HALO_MAP[elem], locant))
                continue
            if elem == "C":
                halo_name = halomethyl_substituent_name(view, nbr, chain_set)
                if halo_name is not None:
                    substituents.append(Sub(halo_name, locant))
                    continue
                ring_name = ring_substituent_name(view, nbr, chain_set, ring_ctx)
                if ring_name is not None:
                    substituents.append(Sub(ring_name, locant))
                    continue
                name = alkyl_substituent_name(view, nbr, chain_set)
                substituents.append(Sub(name, locant))
                continue
            if elem == "O":
                name = alkoxy_substituent_name(view, nbr, chain_set)
                if name is not None:
                    substituents.append(Sub(name, locant))
                    continue
            if elem == "N":
                name = nitro_substituent_name(view, nbr, chain_set)
                if name is not None:
                    substituents.append(Sub(name, locant))
                    continue
                name = amino_substituent_name(view, nbr, chain_set)
                if name is not None:
                    substituents.append(Sub(name, locant))
                    continue
            raise ChemNameNotSupported("Unsupported substituent")

    return substituents


def choose_orientation(view: MolView, chain: Sequence[int], opts) -> List[int]:
    """Choose chain orientation based on substituent locants."""
    forward = list(chain)
    reverse = list(reversed(chain))

    forward_subs = substituents_on_chain(view, forward)
    reverse_subs = substituents_on_chain(view, reverse)

    if opts.prefer_alphabetical_tiebreak:
        f_key = _alphabetical_key(forward_subs)
        r_key = _alphabetical_key(reverse_subs)
    else:
        f_key = _locant_key(forward_subs)
        r_key = _locant_key(reverse_subs)

    if r_key < f_key:
        return reverse
    return forward


def _alphabetical_key(subs: Iterable[Sub]) -> Tuple[Tuple[str, int], ...]:
    return tuple(sorted(((sub.name, sub.locant) for sub in subs)))


def _locant_key(subs: Iterable[Sub]) -> Tuple[int, ...]:
    return tuple(sorted(sub.locant for sub in subs))


def orientation_key(
    subs: Iterable[Sub],
    opts,
    primary_locants: Iterable[int] = (),
    secondary_locants: Iterable[int] = (),
) -> Tuple[Tuple[int, ...], Tuple[int, ...], Tuple]:
    primary = tuple(sorted(primary_locants))
    secondary = tuple(sorted(secondary_locants))
    if opts.prefer_alphabetical_tiebreak:
        sub_key = _alphabetical_key(subs)
    else:
        sub_key = _locant_key(subs)
    return primary, secondary, sub_key
