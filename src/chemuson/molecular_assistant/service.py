"""Strict response decoding and isolated ChemIO validation orchestration."""

from __future__ import annotations

import json
from collections.abc import Sequence

from chemuson.chemio.rdkit_safe import smiles_to_molgraph_isolated
from chemuson.core.model import MolGraph
from chemuson.molecular_assistant.limits import (
    CHEMIO_VALIDATION_TIMEOUT_S,
    MAX_MODEL_CONTENT_BYTES,
    MAX_PROMPT_BYTES,
    MAX_SMILES_BYTES,
)
from chemuson.molecular_assistant.models import (
    MolecularAssistantRequest,
    MolecularAssistantResult,
    MolecularAssistantStatus,
    MolecularTransformationRequest,
    ProviderResponse,
)
from chemuson.molecular_assistant.provider import (
    MolecularStructureProvider,
    ProviderCancelled,
    ProviderError,
    ProviderErrorCode,
)


class _DuplicateJsonKey(ValueError):
    pass


def _unique_object(pairs: Sequence[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise _DuplicateJsonKey
        result[key] = value
    return result


def _reject_json_constant(_value: str) -> None:
    raise ValueError("non-standard JSON constant")


class MolecularAssistant:
    """Propose one structure, decode exactly one SMILES field, validate via M01."""

    def __init__(self, provider: MolecularStructureProvider) -> None:
        self.provider = provider

    def transform(
        self,
        request: MolecularTransformationRequest,
    ) -> MolecularAssistantResult:
        """Transform one source SMILES through the same strict generation/validation path."""
        if (
            not isinstance(request, MolecularTransformationRequest)
            or not isinstance(request.source_smiles, str)
            or not request.source_smiles.strip()
            or not isinstance(request.instruction, str)
            or not request.instruction.strip()
        ):
            return self._failure(
                MolecularAssistantStatus.INVALID_REQUEST,
                "invalid_transformation",
                provider_id=self._provider_id(),
                model_id=self._configured_model_id(),
            )
        description = (
            "Transform the complete molecule represented by this source SMILES. "
            "Return one complete replacement molecule as a SMILES proposal.\n\n"
            f"Source SMILES: {request.source_smiles.strip()}\n"
            f"Requested transformation: {request.instruction.strip()}"
        )
        return self.generate(MolecularAssistantRequest(description))

    def generate(self, request: MolecularAssistantRequest) -> MolecularAssistantResult:
        """Generate a proposal without logging input or mutating application state."""
        provider_id = self._provider_id()
        configured_model_id = self._configured_model_id()

        if not isinstance(request, MolecularAssistantRequest) or not isinstance(request.description, str):
            return self._failure(
                MolecularAssistantStatus.INVALID_REQUEST,
                "invalid_prompt",
                provider_id=provider_id,
                model_id=configured_model_id,
            )
        try:
            prompt_bytes = request.description.encode("utf-8")
        except UnicodeEncodeError:
            return self._failure(
                MolecularAssistantStatus.INVALID_REQUEST,
                "invalid_prompt",
                provider_id=provider_id,
                model_id=configured_model_id,
            )
        if not request.description.strip():
            return self._failure(
                MolecularAssistantStatus.INVALID_REQUEST,
                "empty_prompt",
                provider_id=provider_id,
                model_id=configured_model_id,
            )
        if len(prompt_bytes) > MAX_PROMPT_BYTES:
            return self._failure(
                MolecularAssistantStatus.INVALID_REQUEST,
                "request_too_large",
                provider_id=provider_id,
                model_id=configured_model_id,
            )

        try:
            response = self.provider.generate(request)
        except ProviderCancelled:
            return self._failure(
                MolecularAssistantStatus.CANCELLED,
                "cancelled",
                provider_id=provider_id,
                model_id=configured_model_id,
            )
        except ProviderError as exc:
            if exc.code == ProviderErrorCode.RESPONSE_TOO_LARGE:
                return self._failure(
                    MolecularAssistantStatus.MALFORMED_RESPONSE,
                    "response_too_large",
                    provider_id=provider_id,
                    model_id=configured_model_id,
                )
            reason_code = {
                ProviderErrorCode.TIMEOUT: "timeout",
                ProviderErrorCode.NETWORK_ERROR: "network_error",
                ProviderErrorCode.HTTP_ERROR: "http_error",
            }.get(exc.code, "provider_error")
            return self._failure(
                MolecularAssistantStatus.PROVIDER_ERROR,
                reason_code,
                provider_id=provider_id,
                model_id=configured_model_id,
            )
        except Exception:
            return self._failure(
                MolecularAssistantStatus.PROVIDER_ERROR,
                "provider_error",
                provider_id=provider_id,
                model_id=configured_model_id,
            )

        if (
            not isinstance(response, ProviderResponse)
            or not isinstance(response.content, str)
            or (response.model_id is not None and not isinstance(response.model_id, str))
        ):
            return self._failure(
                MolecularAssistantStatus.PROVIDER_ERROR,
                "provider_error",
                provider_id=provider_id,
                model_id=configured_model_id,
            )
        model_id = response.model_id.strip() if response.model_id else configured_model_id

        smiles, response_reason = self._decode_smiles(response.content)
        if response_reason is not None:
            return self._failure(
                MolecularAssistantStatus.MALFORMED_RESPONSE,
                response_reason,
                provider_id=provider_id,
                model_id=model_id,
            )

        try:
            graph, parser_error = smiles_to_molgraph_isolated(
                smiles,
                timeout_s=CHEMIO_VALIDATION_TIMEOUT_S,
            )
        except Exception:
            graph, parser_error = None, None

        if parser_error == "invalid_input":
            return self._failure(
                MolecularAssistantStatus.INVALID_STRUCTURE,
                "invalid_smiles",
                provider_id=provider_id,
                model_id=model_id,
                proposed_smiles=smiles,
                validation_passed=False,
            )
        if parser_error == "timeout":
            return self._failure(
                MolecularAssistantStatus.VALIDATION_ERROR,
                "parser_timeout",
                provider_id=provider_id,
                model_id=model_id,
                proposed_smiles=smiles,
            )
        if parser_error == "rdkit_unavailable":
            return self._failure(
                MolecularAssistantStatus.VALIDATION_ERROR,
                "parser_unavailable",
                provider_id=provider_id,
                model_id=model_id,
                proposed_smiles=smiles,
            )
        if parser_error is not None or not isinstance(graph, MolGraph) or not graph.atoms:
            return self._failure(
                MolecularAssistantStatus.VALIDATION_ERROR,
                "parser_error",
                provider_id=provider_id,
                model_id=model_id,
                proposed_smiles=smiles,
            )

        return MolecularAssistantResult(
            status=MolecularAssistantStatus.SUCCESS,
            provider_id=provider_id,
            model_id=model_id,
            proposed_smiles=smiles,
            graph=graph,
            validation_passed=True,
        )

    def _provider_id(self) -> str:
        value = getattr(self.provider, "provider_id", None)
        return value if isinstance(value, str) and value.strip() else "custom"

    def _configured_model_id(self) -> str | None:
        value = getattr(self.provider, "model_id", None)
        return value if isinstance(value, str) and value.strip() else None

    @staticmethod
    def _decode_smiles(content: str) -> tuple[str | None, str | None]:
        try:
            content_size = len(content.encode("utf-8"))
        except UnicodeEncodeError:
            return None, "invalid_json"
        if content_size > MAX_MODEL_CONTENT_BYTES:
            return None, "response_too_large"

        try:
            payload = json.loads(
                content,
                object_pairs_hook=_unique_object,
                parse_constant=_reject_json_constant,
            )
        except (json.JSONDecodeError, ValueError, RecursionError):
            return None, "invalid_json"
        if not isinstance(payload, dict):
            return None, "invalid_json"
        if "smiles" not in payload:
            return None, "missing_smiles"
        if set(payload) != {"smiles"}:
            return None, "unexpected_fields"
        smiles = payload["smiles"]
        if not isinstance(smiles, str):
            return None, "invalid_json"
        try:
            smiles_size = len(smiles.encode("utf-8"))
        except UnicodeEncodeError:
            return None, "invalid_json"
        if smiles_size > MAX_SMILES_BYTES:
            return None, "response_too_large"
        if not smiles.strip():
            return None, "empty_smiles"
        return smiles, None

    @staticmethod
    def _failure(
        status: MolecularAssistantStatus,
        reason_code: str,
        *,
        provider_id: str,
        model_id: str | None,
        proposed_smiles: str | None = None,
        validation_passed: bool | None = None,
    ) -> MolecularAssistantResult:
        return MolecularAssistantResult(
            status=status,
            provider_id=provider_id,
            model_id=model_id,
            proposed_smiles=proposed_smiles,
            validation_passed=validation_passed,
            reason_code=reason_code,
        )
