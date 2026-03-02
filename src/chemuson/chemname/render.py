"""Renderizado final de nombres IUPAC-lite."""

from __future__ import annotations

from collections import defaultdict
import re
from typing import Dict, Iterable, List

from .errors import ChemNameNotSupported
from .locants import Sub

# Prefijos multiplicativos para sustituyentes idénticos.
MULTIPLIER = {
    2: "di",
    3: "tri",
    4: "tetra",
    5: "penta",
    6: "hexa",
    7: "hepta",
    8: "octa",
    9: "nona",
    10: "deca",
}


def render_name(
    substituents: Iterable[Sub],
    parent: str,
    always_include_locant: bool = True,
    stereo_descriptors: Iterable[str] | None = None,
) -> str:
    """Renderiza el nombre combinando sustituyentes y padre.

    Args:
        substituents: Iterable de sustituyentes con locantes.
        parent: Nombre del padre (cadena o anillo principal).
        always_include_locant: Forzar locante incluso si es 1 único.
        stereo_descriptors: Lista de designadores estereoquímicos (p. ej., 2R, 3S, E).

    Returns:
        Nombre final en formato IUPAC-lite.

    Raises:
        ChemNameNotSupported: Si hay demasiados sustituyentes idénticos.
    """
    stereo_tokens = _sort_stereo_tokens(stereo_descriptors or [])
    stereo_prefix = f"({','.join(stereo_tokens)})-" if stereo_tokens else ""

    groups: Dict[str, List[int]] = defaultdict(list)
    for sub in substituents:
        groups[sub.name].append(sub.locant)

    if not groups:
        return f"{stereo_prefix}{parent}"

    blocks: List[str] = []
    for name in sorted(groups.keys()):
        locants = sorted(groups[name])
        if len(locants) == 1:
            if not always_include_locant and locants[0] == 1 and len(groups) == 1:
                blocks.append(f"{name}")
            else:
                blocks.append(f"{locants[0]}-{name}")
        else:
            prefix = MULTIPLIER.get(len(locants))
            if prefix is None:
                raise ChemNameNotSupported("Too many identical substituents")
            locant_str = ",".join(str(loc) for loc in locants)
            blocks.append(f"{locant_str}-{prefix}{name}")

    return f"{stereo_prefix}{'-'.join(blocks)}{parent}"


def _sort_stereo_tokens(tokens: Iterable[str]) -> list[str]:
    """Ordena descriptores: helicidad > axialidad > R/S/E/Z > endo/exo/si/re."""
    clean = [str(token).strip() for token in tokens if str(token).strip()]
    clean.sort(key=_stereo_sort_key)
    return clean


def _stereo_sort_key(token: str) -> tuple[int, int, int, str]:
    """Devuelve clave estable para ordenar descriptores estereo."""
    normalized = str(token).strip()
    if normalized in {"M", "P"}:
        order = {"M": 0, "P": 1}
        return (0, 0, order.get(normalized, 9), normalized)

    match = re.fullmatch(
        r"(?P<loc>\d+)(?P<label>R_a|S_a|R|S|E|Z|endo|exo|si|re)",
        normalized,
    )
    if match:
        loc = int(match.group("loc"))
        label = match.group("label")
        if label in {"R_a", "S_a"}:
            order = {"R_a": 0, "S_a": 1}
            return (1, loc, order.get(label, 9), normalized)
        if label in {"R", "S", "E", "Z"}:
            order = {"E": 0, "Z": 1, "R": 2, "S": 3}
            return (2, loc, order.get(label, 9), normalized)
        order = {"endo": 0, "exo": 1, "si": 2, "re": 3}
        return (3, loc, order.get(label, 9), normalized)

    if normalized in {"R_a", "S_a"}:
        order = {"R_a": 0, "S_a": 1}
        return (1, 0, order.get(normalized, 9), normalized)
    if normalized in {"R", "S", "E", "Z"}:
        order = {"E": 0, "Z": 1, "R": 2, "S": 3}
        return (2, 0, order.get(normalized, 9), normalized)
    if normalized in {"endo", "exo", "si", "re"}:
        order = {"endo": 0, "exo": 1, "si": 2, "re": 3}
        return (3, 0, order.get(normalized, 9), normalized)

    return (4, 0, 0, normalized)
