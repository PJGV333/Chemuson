"""Paquete principal de Chemuson.

Agrupa los submódulos de GUI, lógica química y utilidades.
"""

from chemuson._version import __version__
from chemuson.version import get_app_version

__all__ = [
    "__version__",
    "get_app_version",
]
