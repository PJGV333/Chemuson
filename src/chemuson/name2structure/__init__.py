"""Conversión Name→Structure con conectores desacoplados."""

from .service import (
    NameToStructureResult,
    PubChemNameConnector,
    StaticNameConnector,
    resolve_name_to_structure,
)

__all__ = [
    "NameToStructureResult",
    "PubChemNameConnector",
    "StaticNameConnector",
    "resolve_name_to_structure",
]
