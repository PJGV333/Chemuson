"""IUPAC-lite naming for selected hydrocarbons and simple aromatics.

Supported (v0):
- Acyclic carbon chains (alkanes), halogens, linear alkyl substituents.
- Single C=C or C#C within the parent chain.
- Single -ol or -amine functional group.
- Simple cycloalkanes (C3–C8) with simple substituents.
- Benzene with simple substituents (halogens, alkyl, hydroxy, nitro).
- Heteroaromatics: furan, thiophene, pyrrole, pyridine, diazines, imidazole, oxazole, thiazole.
- Naphthalene with mono/di-substituted derivatives (simple substituents).
- Branched alkyl substituents: isopropyl, sec-butyl, tert-butyl, neopentyl.
- Ring substituents: phenyl, benzyl, cyclohexyl.
- Fused aromatics: anthracene, phenanthrene (mono-substituted).
- Fused heteroaromatics: quinoline, isoquinoline, indole (mono-substituted).
Outside scope returns "N/D".
"""

from .engine import iupac_name
from .options import NameOptions

__all__ = ["iupac_name", "NameOptions"]
