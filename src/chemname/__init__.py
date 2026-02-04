"""IUPAC-lite naming for acyclic hydrocarbons.

Supported (v0):
- Acyclic carbon chains (alkanes), halogens, linear alkyl substituents.
- Single C=C or C#C within the parent chain.
- Single -ol or -amine functional group.
- Simple cycloalkanes (C3–C8) with simple substituents.
- Benzene with simple substituents (halogens, alkyl, hydroxy, nitro).
- Heteroaromatics: furan, thiophene, pyrrole, pyridine, diazines, imidazole, oxazole, thiazole.
- Naphthalene and mono-substituted derivatives (simple substituents).
Outside scope returns "N/D".
"""

from .engine import iupac_name
from .options import NameOptions

__all__ = ["iupac_name", "NameOptions"]
