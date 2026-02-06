"""Cálculo de locantes y sustituyentes para la nomenclatura.

Este módulo identifica sustituyentes en una cadena principal y define
criterios de orientación para minimizar locantes.
"""

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
    """Representa un sustituyente con su nombre y locante."""
    name: str
    locant: int


def substituents_on_chain(
    view: MolView,
    chain: Sequence[int],
    ignore_atoms: Optional[Set[int]] = None,
    ring_ctx: RingContext | None = None,
) -> List[Sub]:
    """Devuelve los sustituyentes conectados a la cadena indicada.

    Args:
        view: Vista del grafo molecular.
        chain: Secuencia de IDs que forman la cadena principal.
        ignore_atoms: Conjunto de átomos a ignorar (p. ej., grupo funcional).
        ring_ctx: Contexto de anillos para detectar sustituyentes aromáticos.

    Returns:
        Lista de sustituyentes con su locante en la cadena.

    Raises:
        ChemNameNotSupported: Si se encuentra un sustituyente no soportado.
    """
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
                # Priorizamos sustituyentes especiales antes del alquilo genérico.
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
                # Nitro tiene prioridad sobre amino cuando aplica.
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
    """Elige la orientación de la cadena según reglas de locantes.

    Args:
        view: Vista del grafo molecular.
        chain: Cadena principal candidata.
        opts: Opciones del motor (desempate alfabético o por locantes).

    Returns:
        Cadena orientada (lista de IDs) que minimiza locantes.
    """
    forward = list(chain)
    reverse = list(reversed(chain))

    forward_subs = substituents_on_chain(view, forward)
    reverse_subs = substituents_on_chain(view, reverse)

    # Desempate: alfabético por sustituyente o por locantes según opciones.
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
    """Clave de orden basada en (nombre, locante) para desempate."""
    return tuple(sorted(((sub.name, sub.locant) for sub in subs)))


def _locant_key(subs: Iterable[Sub]) -> Tuple[int, ...]:
    """Clave de orden basada solo en locantes."""
    return tuple(sorted(sub.locant for sub in subs))


def orientation_key(
    subs: Iterable[Sub],
    opts,
    primary_locants: Iterable[int] = (),
    secondary_locants: Iterable[int] = (),
) -> Tuple[Tuple[int, ...], Tuple[int, ...], Tuple]:
    """Crea una clave compuesta de orientación para comparar alternativas.

    Args:
        subs: Sustituyentes detectados.
        opts: Opciones de desempate.
        primary_locants: Locantes prioritarios (p. ej., grupo funcional).
        secondary_locants: Locantes secundarios (p. ej., insaturaciones).

    Returns:
        Tupla comparable que minimiza locantes y aplica desempates.
    """
    primary = tuple(sorted(primary_locants))
    secondary = tuple(sorted(secondary_locants))
    if opts.prefer_alphabetical_tiebreak:
        sub_key = _alphabetical_key(subs)
    else:
        sub_key = _locant_key(subs)
    return primary, secondary, sub_key
