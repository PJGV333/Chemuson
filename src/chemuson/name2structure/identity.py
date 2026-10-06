"""Conservative identity verification for explicit molecule-name requests."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import re
from typing import Protocol

from chemuson.core.model import MolGraph
from chemuson.name2structure.service import (
    NameToStructureResult,
    resolve_name_to_structure,
)


class MolecularIdentityStatus(str, Enum):
    NOT_APPLICABLE = "not_applicable"
    UNVERIFIED = "unverified"
    VERIFIED = "verified"
    MISMATCH = "mismatch"
    REFERENCE_ERROR = "reference_error"


@dataclass(frozen=True, slots=True)
class MolecularIdentityVerification:
    status: MolecularIdentityStatus
    requested_name: str | None = None
    reference_identifier: str | None = None
    reason_code: str | None = None


class IdentityResolver(Protocol):
    """Name resolver contract with explicit, per-call network permission."""

    def __call__(
        self,
        name: str,
        *,
        allow_network: bool,
    ) -> NameToStructureResult: ...

_NAME_REQUEST = re.compile(
    r"^\s*(?:draw|generate|dibuja|genera)(?:\s+(?:the|la|el))?\s+(.+?)\s*[.!?]?\s*$",
    re.IGNORECASE,
)
_OPEN_ENDED_PREFIXES = (
    "a molecule ",
    "a compound ",
    "a structure ",
    "molecule ",
    "compound ",
    "structure ",
    "una molecula ",
    "una molécula ",
    "un compuesto ",
    "una estructura ",
)
_OPEN_ENDED_CUES = (
    " with ",
    " having ",
    " containing ",
    " con ",
    " que tenga ",
    " que contiene ",
)


def extract_requested_molecule_name(request: str) -> str | None:
    """Extract a subject only from a small, explicit name-request grammar."""
    if not isinstance(request, str) or len(request) > 512:
        return None
    match = _NAME_REQUEST.fullmatch(request)
    if match is None:
        return None
    name = match.group(1).strip().strip(".?! ")
    folded = f" {name.casefold()} "
    if not name or len(name) > 160 or "\n" in name or "\r" in name:
        return None
    if any(name.casefold().startswith(prefix) for prefix in _OPEN_ENDED_PREFIXES):
        return None
    if any(cue in folded for cue in _OPEN_ENDED_CUES):
        return None
    return name


def verify_molecular_identity(
    request: str,
    proposed_graph: MolGraph,
    *,
    resolver: IdentityResolver | None = None,
    allow_network: bool = False,
    enabled: bool = True,
) -> MolecularIdentityVerification:
    """Compare isolated InChI identity using an explicit network policy.

    Offline lookup is the default. Injected resolvers receive ``allow_network``
    as a keyword argument just like the built-in Name→Structure resolver.
    """
    if enabled is not True:
        return MolecularIdentityVerification(
            MolecularIdentityStatus.NOT_APPLICABLE,
            reason_code="verification_disabled",
        )
    allow_network = allow_network is True
    requested_name = extract_requested_molecule_name(request)
    if requested_name is None:
        return MolecularIdentityVerification(MolecularIdentityStatus.NOT_APPLICABLE)
    if not isinstance(proposed_graph, MolGraph) or not proposed_graph.atoms:
        return MolecularIdentityVerification(
            MolecularIdentityStatus.REFERENCE_ERROR,
            requested_name=requested_name,
            reason_code="proposal_unavailable",
        )

    try:
        reference = (
            resolver(requested_name, allow_network=allow_network)
            if resolver is not None
            else resolve_name_to_structure(
                requested_name,
                allow_network=allow_network,
                timeout_s=8.0,
            )
        )
    except Exception:
        return MolecularIdentityVerification(
            MolecularIdentityStatus.REFERENCE_ERROR,
            requested_name=requested_name,
            reason_code="resolver_error",
        )

    if not isinstance(reference, NameToStructureResult) or not reference.ok or reference.graph is None:
        message = str(getattr(reference, "message", "") or "").casefold()
        failed = any(token in message for token in ("timeout", "urlerror", "httperror", "network", "conversion_failed"))
        status = (
            MolecularIdentityStatus.REFERENCE_ERROR
            if failed
            else MolecularIdentityStatus.UNVERIFIED
        )
        return MolecularIdentityVerification(
            status,
            requested_name=requested_name,
            reason_code=(
                "reference_unavailable"
                if failed
                else "reference_not_found"
                if allow_network
                else "reference_not_found_offline"
            ),
        )
    if reference.confidence < 0.7:
        return MolecularIdentityVerification(
            MolecularIdentityStatus.UNVERIFIED,
            requested_name=requested_name,
            reason_code="reference_confidence_low",
        )

    try:
        from chemuson.chemio.rdkit_safe import molgraph_to_inchi_isolated

        proposed_inchi, proposed_error = molgraph_to_inchi_isolated(
            proposed_graph,
            timeout_s=5.0,
        )
        reference_inchi, reference_error = molgraph_to_inchi_isolated(
            reference.graph,
            timeout_s=5.0,
        )
    except Exception:
        proposed_inchi = reference_inchi = None
        proposed_error = reference_error = "worker_error"
    if (
        proposed_error
        or reference_error
        or not proposed_inchi
        or not reference_inchi
    ):
        return MolecularIdentityVerification(
            MolecularIdentityStatus.REFERENCE_ERROR,
            requested_name=requested_name,
            reason_code="canonicalization_error",
        )

    resolved = reference.resolved_name or reference.query or requested_name
    identifier = f"{reference.source}:{resolved}"
    return MolecularIdentityVerification(
        MolecularIdentityStatus.VERIFIED
        if proposed_inchi == reference_inchi
        else MolecularIdentityStatus.MISMATCH,
        requested_name=requested_name,
        reference_identifier=identifier,
    )
