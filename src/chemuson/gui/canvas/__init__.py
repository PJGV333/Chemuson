"""API publica compatible para el paquete canvas de Chemuson."""
from __future__ import annotations

from . import _shared as _shared_module
from .canvas_view import ChemusonCanvas

_SHARED_EXPORTS = tuple(
    name
    for name in dir(_shared_module)
    if not name.startswith("_")
)

globals().update({name: getattr(_shared_module, name) for name in _SHARED_EXPORTS})

__all__ = [*_SHARED_EXPORTS, "ChemusonCanvas"]
