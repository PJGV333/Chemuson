"""Conversión Name→Structure con conectores desacoplados."""

from .identity import (
    MolecularIdentityStatus,
    MolecularIdentityVerification,
    extract_requested_molecule_name,
    verify_molecular_identity,
)
from .service import (
    NameToStructureResult,
    PubChemNameConnector,
    StaticNameConnector,
    resolve_name_to_structure,
)

__all__ = [
    "MolecularIdentityStatus",
    "MolecularIdentityVerification",
    "extract_requested_molecule_name",
    "verify_molecular_identity",
    "NameToStructureResult",
    "PubChemNameConnector",
    "StaticNameConnector",
    "resolve_name_to_structure",
]
