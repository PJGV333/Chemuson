"""Renderizado final de nombres IUPAC-lite."""

from __future__ import annotations

from collections import defaultdict
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
    stereo_tokens = [str(token).strip() for token in (stereo_descriptors or []) if str(token).strip()]
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
