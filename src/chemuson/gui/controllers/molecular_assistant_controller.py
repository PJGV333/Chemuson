"""Asynchronous GUI adapter for the provider-neutral molecular assistant."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, replace
from enum import Enum
import math
from typing import Any

from PyQt6.QtCore import QObject, QThread, pyqtSignal, pyqtSlot

from chemuson.chemio.rdkit_io import molgraph_to_smiles_isolated_or_error
from chemuson.chemio.rdkit_safe import smiles_to_molgraph_isolated
from chemuson.core.model import MolGraph
from chemuson.name2structure import (
    MolecularIdentityStatus,
    MolecularIdentityVerification,
    NameToStructureResult,
    extract_requested_molecule_name,
    resolve_name_to_structure,
    verify_molecular_identity,
)
from chemuson.molecular_assistant import (
    MolecularAssistant,
    MolecularAssistantRequest,
    MolecularAssistantResult,
    MolecularTransformationRequest,
    MolecularAssistantStatus,
    OpenAICompatibleConfig,
    OPENAI_COMPATIBLE_PROFILES,
    OpenAICompatibleProvider,
    get_openai_compatible_profile,
)


class MolecularResolutionMethod(str, Enum):
    AI = "ai"
    AI_REFERENCE = "ai_reference"
    REFERENCE = "reference"


class StructureOrigin(str, Enum):
    AI = "ai"
    REFERENCE = "reference"
    AI_VERIFIED_BY_REFERENCE = "ai_verified_by_reference"
    AI_MISMATCH_REFERENCE = "ai_mismatch_reference"


AssistantRequest = MolecularAssistantRequest | MolecularTransformationRequest
ResultGenerator = Callable[[AssistantRequest, OpenAICompatibleConfig], MolecularAssistantResult]
ReferenceResolver = Callable[..., NameToStructureResult]
IdentityVerifier = Callable[
    [str, MolecularAssistantResult, NameToStructureResult], MolecularIdentityVerification
]


def _generate_structure(
    request: AssistantRequest,
    config: OpenAICompatibleConfig,
) -> MolecularAssistantResult:
    assistant = MolecularAssistant(OpenAICompatibleProvider(config))
    if isinstance(request, MolecularTransformationRequest):
        return assistant.transform(request)
    return assistant.generate(request)


def _provider_failure(config: OpenAICompatibleConfig) -> MolecularAssistantResult:
    return MolecularAssistantResult(
        status=MolecularAssistantStatus.PROVIDER_ERROR,
        provider_id=config.provider_id,
        model_id=config.model,
        reason_code="provider_error",
    )


def _resolve_reference(name: str, allow_network: bool) -> NameToStructureResult:
    return resolve_name_to_structure(name, allow_network=allow_network, timeout_s=8.0)


def _verify_reference_identity(
    request: str,
    result: MolecularAssistantResult,
    reference: NameToStructureResult,
) -> MolecularIdentityVerification:
    if result.graph is None:
        return MolecularIdentityVerification(
            MolecularIdentityStatus.REFERENCE_ERROR,
            requested_name=extract_requested_molecule_name(request),
            reason_code="proposal_unavailable",
        )
    return verify_molecular_identity(
        request,
        result.graph,
        allow_network=False,
        resolver=lambda _name, *, allow_network: reference,
    )


@dataclass(frozen=True)
class MolecularAssistantResolution:
    """Transient AI/reference candidates and their controlled provenance."""

    method: MolecularResolutionMethod
    ai_result: MolecularAssistantResult | None
    reference_result: NameToStructureResult | None
    identity: MolecularIdentityVerification
    default_candidate: str | None
    origin: StructureOrigin | None
    failure_reason: str | None = None

    def graph_for_candidate(self, candidate: str | None) -> MolGraph | None:
        if candidate == "ai" and self.ai_result is not None:
            return self.ai_result.graph
        if candidate == "reference" and self.reference_result is not None:
            return self.reference_result.graph
        return None

    def smiles_for_candidate(self, candidate: str | None) -> str | None:
        if candidate == "ai" and self.ai_result is not None:
            return self.ai_result.proposed_smiles
        if candidate == "reference" and self.reference_result is not None:
            return self.reference_result.smiles or None
        return None

    def origin_for_candidate(self, candidate: str | None) -> StructureOrigin | None:
        if candidate == "reference":
            return StructureOrigin.REFERENCE
        if candidate != "ai":
            return None
        if self.identity.status is MolecularIdentityStatus.MISMATCH:
            return StructureOrigin.AI_MISMATCH_REFERENCE
        if self.identity.status is MolecularIdentityStatus.VERIFIED:
            return StructureOrigin.AI_VERIFIED_BY_REFERENCE
        return StructureOrigin.AI

    @property
    def graph(self) -> MolGraph | None:
        return self.graph_for_candidate(self.default_candidate)

    @property
    def status(self) -> MolecularAssistantStatus:
        if self.graph is not None:
            return MolecularAssistantStatus.SUCCESS
        if self.ai_result is not None and self.ai_result.status is not MolecularAssistantStatus.SUCCESS:
            return self.ai_result.status
        return MolecularAssistantStatus.PROVIDER_ERROR

    @property
    def reason_code(self) -> str | None:
        if self.graph is not None:
            return None
        if self.ai_result is not None and self.ai_result.status is not MolecularAssistantStatus.SUCCESS:
            return self.ai_result.reason_code
        return self.failure_reason or "reference_unavailable"

    @property
    def validation_passed(self) -> bool | None:
        if self.graph is not None:
            return True
        return self.ai_result.validation_passed if self.ai_result is not None else None

    @property
    def proposed_smiles(self) -> str | None:
        return self.smiles_for_candidate(self.default_candidate)

    @property
    def provider_id(self) -> str:
        if self.default_candidate == "ai" and self.ai_result is not None:
            return self.ai_result.provider_id
        return "reference"

    @property
    def model_id(self) -> str | None:
        if self.default_candidate == "ai" and self.ai_result is not None:
            return self.ai_result.model_id
        return None

    @property
    def ai_failure_reason(self) -> str | None:
        if self.ai_result is None or self.ai_result.status is MolecularAssistantStatus.SUCCESS:
            return None
        return self.ai_result.reason_code


class _MolecularAssistantWorker(QObject):
    """Run the selected generation/reference pipeline on the existing QThread."""

    finished = pyqtSignal(int, object, str, object)

    def __init__(
        self,
        job_id: int,
        description: str,
        config: OpenAICompatibleConfig | None,
        generator: ResultGenerator,
        source_graph: MolGraph | None,
        reference_resolver: ReferenceResolver,
        identity_verifier: IdentityVerifier,
        resolution_method: MolecularResolutionMethod,
        allow_external_reference: bool,
    ) -> None:
        super().__init__()
        self._job_id = int(job_id)
        self._description = description
        self._config = config
        self._generator = generator
        self._source_graph = source_graph
        self._reference_resolver = reference_resolver
        self._identity_verifier = identity_verifier
        self._resolution_method = resolution_method
        self._allow_external_reference = allow_external_reference

    def _resolve_reference(self, name: str) -> NameToStructureResult:
        try:
            result = self._reference_resolver(
                name,
                allow_network=self._allow_external_reference,
            )
        except Exception:
            return NameToStructureResult(
                name, None, "none", 0.0, message="reference_error"
            )
        if not isinstance(result, NameToStructureResult):
            return NameToStructureResult(
                name, None, "none", 0.0, message="reference_error"
            )
        if result.graph is None:
            message = str(result.message or "").strip().casefold()
            reason = (
                "reference_not_found"
                if message in {"not_found", "empty_query"}
                else "reference_error"
            )
            return NameToStructureResult(
                name, None, "none", 0.0, message=reason
            )
        smiles = result.smiles if isinstance(result.smiles, str) else ""
        try:
            smiles_size = len(smiles.encode("utf-8"))
        except UnicodeEncodeError:
            smiles_size = 16 * 1024 + 1
        if not smiles.strip() or smiles_size > 16 * 1024:
            return NameToStructureResult(
                name, None, "none", 0.0, message="invalid_reference"
            )
        try:
            confidence = float(result.confidence)
        except (TypeError, ValueError, OverflowError):
            confidence = -1.0
        if not math.isfinite(confidence) or not 0.0 <= confidence <= 1.0:
            return NameToStructureResult(
                name, None, "none", 0.0, message="invalid_reference"
            )
        try:
            graph, error = smiles_to_molgraph_isolated(smiles, timeout_s=8.0)
        except Exception:
            graph, error = None, "conversion_failed"
        if error or not isinstance(graph, MolGraph) or not graph.atoms:
            return NameToStructureResult(
                name, None, "none", 0.0, message="invalid_reference"
            )
        return replace(result, graph=graph, confidence=confidence, message="")

    @pyqtSlot()
    def run(self) -> None:
        source_smiles = ""
        request: AssistantRequest = MolecularAssistantRequest(self._description)
        if self._source_graph is not None:
            try:
                source_smiles = molgraph_to_smiles_isolated_or_error(
                    self._source_graph,
                    timeout_s=8.0,
                )
            except Exception:
                source_smiles = ""
            if not isinstance(source_smiles, str) or not source_smiles.strip():
                failure = MolecularAssistantResult(
                    status=MolecularAssistantStatus.VALIDATION_ERROR,
                    provider_id=self._config.provider_id if self._config else "unknown",
                    model_id=self._config.model if self._config else None,
                    reason_code="source_export_failed",
                )
                outcome = MolecularAssistantResolution(
                    MolecularResolutionMethod.AI,
                    failure,
                    None,
                    MolecularIdentityVerification(MolecularIdentityStatus.NOT_APPLICABLE),
                    None,
                    None,
                )
                self.finished.emit(self._job_id, outcome, "", outcome.identity)
                return
            request = MolecularTransformationRequest(
                source_smiles=source_smiles.strip(),
                instruction=self._description,
            )

        ai_result: MolecularAssistantResult | None = None
        if self._resolution_method is not MolecularResolutionMethod.REFERENCE:
            if self._config is None:
                ai_result = None
            else:
                try:
                    ai_result = self._generator(request, self._config)
                except Exception:
                    ai_result = _provider_failure(self._config)
                if not isinstance(ai_result, MolecularAssistantResult):
                    ai_result = _provider_failure(self._config)

        requested_name = (
            None
            if self._source_graph is not None
            else extract_requested_molecule_name(self._description)
        )
        reference: NameToStructureResult | None = None
        if (
            self._source_graph is None
            and requested_name is not None
            and self._resolution_method
            in {MolecularResolutionMethod.AI_REFERENCE, MolecularResolutionMethod.REFERENCE}
            and not QThread.currentThread().isInterruptionRequested()
        ):
            reference = self._resolve_reference(requested_name)

        identity = MolecularIdentityVerification(MolecularIdentityStatus.NOT_APPLICABLE)
        if self._resolution_method is MolecularResolutionMethod.AI_REFERENCE:
            if requested_name is not None and ai_result is not None:
                if reference is not None and reference.graph is not None:
                    try:
                        identity = self._identity_verifier(
                            self._description,
                            ai_result,
                            reference,
                        )
                    except Exception:
                        identity = MolecularIdentityVerification(
                            MolecularIdentityStatus.REFERENCE_ERROR,
                            requested_name=requested_name,
                            reason_code="verification_error",
                        )
                elif ai_result.status is MolecularAssistantStatus.SUCCESS:
                    reference_error = reference is not None and reference.message in {
                        "invalid_reference", "reference_error"
                    }
                    identity = MolecularIdentityVerification(
                        MolecularIdentityStatus.REFERENCE_ERROR
                        if reference_error
                        else MolecularIdentityStatus.UNVERIFIED,
                        requested_name=requested_name,
                        reason_code=(
                            "reference_invalid"
                            if reference is not None and reference.message == "invalid_reference"
                            else "reference_unavailable"
                            if reference_error
                            else "reference_not_found_offline"
                            if not self._allow_external_reference
                            else "reference_not_found"
                        ),
                    )
                else:
                    identity = MolecularIdentityVerification(
                        MolecularIdentityStatus.NOT_APPLICABLE,
                        requested_name=requested_name,
                        reason_code="proposal_unavailable",
                    )
        elif self._resolution_method is MolecularResolutionMethod.AI:
            if requested_name is not None:
                identity = MolecularIdentityVerification(
                    MolecularIdentityStatus.UNVERIFIED,
                    requested_name=requested_name,
                    reason_code="reference_not_requested",
                )
        elif requested_name is None:
            identity = MolecularIdentityVerification(
                MolecularIdentityStatus.NOT_APPLICABLE,
                reason_code="reference_name_required",
            )
        elif reference is not None and reference.graph is not None:
            identity = MolecularIdentityVerification(
                MolecularIdentityStatus.NOT_APPLICABLE,
                requested_name=requested_name,
                reason_code="reference_selected",
            )

        ai_ok = (
            ai_result is not None
            and ai_result.status is MolecularAssistantStatus.SUCCESS
            and ai_result.graph is not None
        )
        reference_ok = reference is not None and reference.graph is not None
        default_candidate: str | None = None
        origin: StructureOrigin | None = None
        failure_reason: str | None = None
        if self._resolution_method is MolecularResolutionMethod.REFERENCE:
            if requested_name is None:
                failure_reason = "reference_name_required"
            elif reference_ok:
                default_candidate = "reference"
                origin = StructureOrigin.REFERENCE
            else:
                failure_reason = self._reference_failure_reason(reference)
        elif ai_ok:
            default_candidate = "ai"
            origin = StructureOrigin.AI
            if reference_ok and identity.status is MolecularIdentityStatus.VERIFIED:
                origin = StructureOrigin.AI_VERIFIED_BY_REFERENCE
            elif reference_ok and identity.status is MolecularIdentityStatus.MISMATCH:
                default_candidate = "reference"
                origin = StructureOrigin.AI_MISMATCH_REFERENCE
        elif reference_ok and self._resolution_method is MolecularResolutionMethod.AI_REFERENCE:
            default_candidate = "reference"
            origin = StructureOrigin.REFERENCE
        elif ai_result is None:
            failure_reason = "reference_unavailable"
        elif ai_result.status is MolecularAssistantStatus.SUCCESS:
            failure_reason = "invalid_structure"
        else:
            failure_reason = ai_result.reason_code or "provider_error"

        outcome = MolecularAssistantResolution(
            method=self._resolution_method,
            ai_result=ai_result,
            reference_result=reference,
            identity=identity,
            default_candidate=default_candidate,
            origin=origin,
            failure_reason=failure_reason,
        )
        self.finished.emit(self._job_id, outcome, source_smiles, identity)

    @staticmethod
    def _reference_failure_reason(reference: NameToStructureResult | None) -> str:
        if reference is None:
            return "reference_unavailable"
        if reference.message == "reference_name_required":
            return "reference_name_required"
        if reference.message == "invalid_reference":
            return "reference_invalid"
        if reference.message == "reference_not_found":
            return "reference_not_found"
        return "reference_unavailable"


class MolecularAssistantController(QObject):
    """Own worker threads and relay only typed, stable M23 results to the GUI."""

    job_started = pyqtSignal(int)
    job_finished = pyqtSignal(int, object)
    source_smiles_ready = pyqtSignal(int, str)
    identity_ready = pyqtSignal(int, object)

    def __init__(
        self,
        parent: QObject | None = None,
        *,
        generator: ResultGenerator | None = None,
        identity_verifier: IdentityVerifier | None = None,
        reference_resolver: ReferenceResolver | None = None,
    ) -> None:
        super().__init__(parent)
        self._generator = generator or _generate_structure
        self._identity_verifier = identity_verifier or _verify_reference_identity
        self._reference_resolver = reference_resolver or _resolve_reference
        self._next_job_id = 1
        self._jobs: dict[int, tuple[QThread, _MolecularAssistantWorker]] = {}
        self._pending_results: dict[int, MolecularAssistantResolution] = {}
        self._pending_source_smiles: dict[int, str] = {}
        self._pending_identities: dict[int, object] = {}
        self._abandoned_jobs: set[int] = set()
        self._shutting_down = False

    @property
    def provider_profiles(self):
        """Expose immutable profile metadata to the UI without a direct M23 import."""
        return OPENAI_COMPATIBLE_PROFILES

    def start_job(
        self,
        description: str,
        *,
        base_url: str = "",
        model: str = "",
        api_key: str | None = None,
        supports_json_output: bool = False,
        provider_id: str = "openai-compatible",
        timeout_s: float = 60.0,
        max_tokens: int = 4096,
        source_graph: MolGraph | None = None,
        resolution_method: MolecularResolutionMethod | str = MolecularResolutionMethod.AI,
        allow_external_reference: bool = False,
    ) -> int | None:
        """Validate settings and start the selected AI/reference route."""
        if self._shutting_down:
            return None
        if not isinstance(description, str) or not description.strip():
            return None
        try:
            method = MolecularResolutionMethod(resolution_method)
        except (TypeError, ValueError):
            return None
        if not isinstance(allow_external_reference, bool):
            return None
        if source_graph is not None and method is not MolecularResolutionMethod.AI:
            return None

        config: OpenAICompatibleConfig | None = None
        if method is not MolecularResolutionMethod.REFERENCE:
            profile = get_openai_compatible_profile(provider_id)
            if profile is None:
                return None
            if profile.api_key_required and (
                not isinstance(api_key, str) or not api_key.strip()
            ):
                return None
            try:
                config = OpenAICompatibleConfig(
                    base_url=base_url.strip() if isinstance(base_url, str) else base_url,
                    model=model.strip() if isinstance(model, str) else model,
                    api_key=api_key if api_key else None,
                    supports_json_output=supports_json_output,
                    provider_id=profile.profile_id,
                    timeout_s=timeout_s,
                    max_tokens=max_tokens,
                )
            except (TypeError, ValueError):
                return None

        job_id = self._next_job_id
        self._next_job_id += 1
        thread = QThread(self)
        thread.setObjectName(f"MolecularAssistant-{job_id}")
        worker = _MolecularAssistantWorker(
            job_id,
            description,
            config,
            self._generator,
            source_graph,
            self._reference_resolver,
            self._identity_verifier,
            method,
            allow_external_reference,
        )
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.finished.connect(self._record_worker_finished)
        worker.finished.connect(thread.quit)
        worker.finished.connect(worker.deleteLater)
        thread.finished.connect(lambda job_id=job_id: self._on_thread_finished(job_id))
        thread.finished.connect(thread.deleteLater)
        self._jobs[job_id] = (thread, worker)
        self.job_started.emit(job_id)
        thread.start()
        return job_id

    def active_jobs(self) -> tuple[int, ...]:
        """Return job identifiers whose QThreads have not completed."""
        return tuple(sorted(self._jobs))

    def has_active_jobs(self) -> bool:
        """Whether any owned QThread is still running."""
        return any(thread.isRunning() for thread, _worker in self._jobs.values())

    @property
    def shutdown_complete(self) -> bool:
        """Whether shutdown began and every owned QThread has stopped."""
        return self._shutting_down and not self.has_active_jobs()

    def begin_shutdown(self) -> None:
        """Suppress late results and request cooperative worker interruption."""
        if self._shutting_down:
            return
        self._shutting_down = True
        for job_id, (thread, _worker) in self._jobs.items():
            self._abandoned_jobs.add(job_id)
            thread.requestInterruption()

    def abandon_job(self, job_id: int) -> None:
        """Suppress late delivery and cooperatively skip follow-up work, not in-flight HTTP."""
        job_id = int(job_id)
        job = self._jobs.get(job_id)
        if job is not None:
            self._abandoned_jobs.add(job_id)
            job[0].requestInterruption()

    @pyqtSlot(int, object, str, object)
    def _record_worker_finished(
        self,
        job_id: int,
        result: Any,
        source_smiles: str,
        identity: object,
    ) -> None:
        if self._shutting_down:
            return
        if isinstance(result, MolecularAssistantResolution):
            self._pending_results[int(job_id)] = result
            if isinstance(source_smiles, str) and source_smiles:
                self._pending_source_smiles[int(job_id)] = source_smiles
            if identity is not None:
                self._pending_identities[int(job_id)] = identity

    def _on_thread_finished(self, job_id: int) -> None:
        job_id = int(job_id)
        self._jobs.pop(job_id, None)
        result = self._pending_results.pop(job_id, None)
        source_smiles = self._pending_source_smiles.pop(job_id, None)
        identity = self._pending_identities.pop(job_id, None)
        if self._shutting_down or job_id in self._abandoned_jobs:
            self._abandoned_jobs.discard(job_id)
            return
        if result is not None:
            if source_smiles is not None:
                self.source_smiles_ready.emit(job_id, source_smiles)
            if identity is not None:
                self.identity_ready.emit(job_id, identity)
            self.job_finished.emit(job_id, result)
