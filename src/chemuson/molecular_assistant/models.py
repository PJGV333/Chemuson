"""Provider-neutral request and result contracts for molecular generation."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from chemuson.core.model import MolGraph


class MolecularAssistantStatus(str, Enum):
    """Stable outcomes exposed by the molecular-assistant service."""

    SUCCESS = "success"
    INVALID_REQUEST = "invalid_request"
    INVALID_STRUCTURE = "invalid_structure"
    VALIDATION_ERROR = "validation_error"
    PROVIDER_ERROR = "provider_error"
    MALFORMED_RESPONSE = "malformed_response"
    CANCELLED = "cancelled"


@dataclass(frozen=True)
class MolecularAssistantRequest:
    """A single natural-language description; no document or canvas context."""

    description: str


@dataclass(frozen=True)
class ProviderResponse:
    """Raw structured-response content and optional transport model identity."""

    content: str
    model_id: str | None = None


@dataclass(frozen=True)
class MolecularAssistantResult:
    """Structured outcome; failed operations can never carry a usable graph."""

    status: MolecularAssistantStatus
    provider_id: str
    model_id: str | None = None
    proposed_smiles: str | None = None
    graph: MolGraph | None = None
    validation_passed: bool | None = None
    reason_code: str | None = None

    def __post_init__(self) -> None:
        if self.status == MolecularAssistantStatus.SUCCESS:
            if not isinstance(self.graph, MolGraph) or not self.graph.atoms:
                raise ValueError("success requires a non-empty MolGraph")
            if self.validation_passed is not True or self.reason_code is not None:
                raise ValueError("success requires validation_passed=true and no reason")
            return

        if self.graph is not None:
            raise ValueError("failed results must not contain a MolGraph")
        if self.reason_code is None:
            raise ValueError("failed results require a stable reason_code")
        if self.status == MolecularAssistantStatus.INVALID_STRUCTURE:
            if self.validation_passed is not False:
                raise ValueError("invalid_structure requires validation_passed=false")
        elif self.validation_passed is not None:
            raise ValueError("non-chemical failures require validation_passed=null")
