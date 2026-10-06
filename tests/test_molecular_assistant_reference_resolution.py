from __future__ import annotations

import threading
import time

import pytest
from PyQt6.QtWidgets import QApplication
from PyQt6.QtTest import QSignalSpy

from chemuson.chemio.rdkit_safe import smiles_to_molgraph_isolated
from chemuson.core.model import MolGraph
from chemuson.gui.controllers import (
    MolecularAssistantController,
    MolecularResolutionMethod,
    StructureOrigin,
)
from chemuson.gui.controllers.molecular_assistant_controller import (
    MolecularAssistantResolution,
)
from chemuson.molecular_assistant import (
    MolecularAssistantResult,
    MolecularAssistantStatus,
    MolecularTransformationRequest,
)
from chemuson.name2structure import (
    MolecularIdentityStatus,
    NameToStructureResult,
)


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def _graph(smiles: str) -> MolGraph:
    graph, error = smiles_to_molgraph_isolated(smiles, timeout_s=5.0)
    if graph is None:
        pytest.skip(f"ChemIO isolated worker unavailable: {error}")
    return graph


def _success(smiles: str, graph: MolGraph | None = None) -> MolecularAssistantResult:
    return MolecularAssistantResult(
        status=MolecularAssistantStatus.SUCCESS,
        provider_id="fake-model-provider",
        model_id="fake-model",
        proposed_smiles=smiles,
        graph=graph or _graph(smiles),
        validation_passed=True,
    )


def _failure(
    reason: str,
    status: MolecularAssistantStatus = MolecularAssistantStatus.PROVIDER_ERROR,
) -> MolecularAssistantResult:
    return MolecularAssistantResult(
        status=status,
        provider_id="fake-model-provider",
        model_id="fake-model",
        reason_code=reason,
        finish_reason="length" if reason == "generation_exhausted" else None,
        completion_tokens=4096 if reason == "generation_exhausted" else None,
        reasoning_tokens=4096 if reason == "generation_exhausted" else None,
    )


def _reference(name: str, smiles: str, *, source: str = "pubchem", from_cache: bool = False):
    return NameToStructureResult(
        query=name,
        graph=_graph(smiles),
        source=source,
        confidence=0.95,
        smiles=smiles,
        resolved_name=name.title(),
        from_cache=from_cache,
    )


def _run_job(
    *,
    prompt: str,
    generator,
    resolver,
    method: MolecularResolutionMethod = MolecularResolutionMethod.AI_REFERENCE,
    allow_external: bool = False,
    source_graph: MolGraph | None = None,
    provider_id: str = "llama-cpp",
) -> MolecularAssistantResolution:
    app = QApplication.instance()
    assert app is not None
    controller = MolecularAssistantController(
        generator=generator,
        reference_resolver=resolver,
    )
    spy = QSignalSpy(controller.job_finished)
    job_id = controller.start_job(
        prompt,
        base_url="http://127.0.0.1:8080/v1" if method is not MolecularResolutionMethod.REFERENCE else "",
        model="offline-test-model" if method is not MolecularResolutionMethod.REFERENCE else "",
        provider_id=provider_id,
        source_graph=source_graph,
        resolution_method=method,
        allow_external_reference=allow_external,
    )
    assert job_id is not None
    assert spy.wait(8_000)
    deadline = time.monotonic() + 3.0
    while controller.active_jobs() and time.monotonic() < deadline:
        app.processEvents()
        time.sleep(0.005)
    assert controller.active_jobs() == ()
    result = spy[0][1]
    assert isinstance(result, MolecularAssistantResolution)
    controller.deleteLater()
    app.processEvents()
    return result


def test_ai_success_with_equivalent_nonidentical_smiles_is_verified():
    proposal = _success("OCC", _graph("CCO"))
    calls = []

    def resolver(name, *, allow_network):
        calls.append((name, allow_network))
        return _reference(name, "CCO", from_cache=True)

    outcome = _run_job(
        prompt="Draw ethanol",
        generator=lambda *_args: proposal,
        resolver=resolver,
    )

    assert calls == [("ethanol", False)]
    assert outcome.status is MolecularAssistantStatus.SUCCESS
    assert outcome.identity.status is MolecularIdentityStatus.VERIFIED
    assert outcome.default_candidate == "ai"
    assert outcome.origin is StructureOrigin.AI_VERIFIED_BY_REFERENCE
    assert outcome.ai_result is proposal
    assert outcome.reference_result.from_cache is True


def test_ai_mismatch_offers_both_unchanged_candidates_and_defaults_to_reference():
    proposal = _success("CCN")
    outcome = _run_job(
        prompt="Dibuja tetrandrina",
        generator=lambda *_args: proposal,
        resolver=lambda name, *, allow_network: _reference(name, "CCO"),
    )

    assert outcome.identity.status is MolecularIdentityStatus.MISMATCH
    assert outcome.default_candidate == "reference"
    assert outcome.origin is StructureOrigin.AI_MISMATCH_REFERENCE
    assert outcome.smiles_for_candidate("ai") == "CCN"
    assert outcome.smiles_for_candidate("reference") == "CCO"
    assert outcome.graph_for_candidate("ai") is proposal.graph
    assert outcome.graph_for_candidate("reference") is outcome.reference_result.graph


@pytest.mark.parametrize(
    ("status", "reason"),
    [
        (MolecularAssistantStatus.MALFORMED_RESPONSE, "invalid_json"),
        (MolecularAssistantStatus.PROVIDER_ERROR, "timeout"),
        (MolecularAssistantStatus.MALFORMED_RESPONSE, "generation_exhausted"),
    ],
)
def test_ai_failure_falls_back_to_validated_reference_and_preserves_ai_diagnostic(status, reason):
    ai_failure = _failure(reason, status)
    outcome = _run_job(
        prompt="Dibuja tetrandrina",
        generator=lambda *_args: ai_failure,
        resolver=lambda name, *, allow_network: _reference(name, "CCO"),
    )

    assert outcome.status is MolecularAssistantStatus.SUCCESS
    assert outcome.ai_result is ai_failure
    assert outcome.ai_failure_reason == reason
    assert outcome.default_candidate == "reference"
    assert outcome.origin is StructureOrigin.REFERENCE
    assert outcome.graph is outcome.reference_result.graph
    assert outcome.reference_result.smiles == "CCO"
    if reason == "generation_exhausted":
        assert outcome.ai_result.finish_reason == "length"
        assert outcome.ai_result.completion_tokens == 4096
        assert outcome.ai_result.reasoning_tokens == 4096


def test_ai_failure_and_missing_reference_remains_a_controlled_failure():
    ai_failure = _failure("timeout")
    outcome = _run_job(
        prompt="Draw caffeine",
        generator=lambda *_args: ai_failure,
        resolver=lambda name, *, allow_network: NameToStructureResult(
            name, None, "pubchem", 0.0, message="not_found"
        ),
    )

    assert outcome.status is MolecularAssistantStatus.PROVIDER_ERROR
    assert outcome.reason_code == "timeout"
    assert outcome.default_candidate is None
    assert outcome.graph is None
    assert outcome.ai_result is ai_failure


def test_ai_success_without_reference_is_unverified_and_stays_ai_origin():
    calls = []
    outcome = _run_job(
        prompt="Draw caffeine",
        generator=lambda *_args: _success("CCO"),
        resolver=lambda name, *, allow_network: calls.append(name)
        or NameToStructureResult(name, None, "none", 0.0, message="not_found"),
    )

    assert calls == ["caffeine"]
    assert outcome.status is MolecularAssistantStatus.SUCCESS
    assert outcome.identity.status is MolecularIdentityStatus.UNVERIFIED
    assert outcome.default_candidate == "ai"
    assert outcome.origin is StructureOrigin.AI


def test_network_permission_is_passed_explicitly_and_only_extracted_name_is_resolved():
    for allow_network in (False, True):
        calls = []
        outcome = _run_job(
            prompt="Dibuja la cafeína",
            generator=lambda *_args: _success("CCO"),
            resolver=lambda name, *, allow_network: calls.append((name, allow_network))
            or _reference(name, "CCO"),
            allow_external=allow_network,
        )
        assert outcome.status is MolecularAssistantStatus.SUCCESS
        assert calls == [("cafeína", allow_network)]


def test_open_ended_prompt_and_ai_only_never_call_reference_resolver():
    for method, prompt in (
        (MolecularResolutionMethod.AI_REFERENCE, "Generate a molecule with three rings"),
        (MolecularResolutionMethod.AI, "Draw caffeine"),
    ):
        calls = []
        outcome = _run_job(
            prompt=prompt,
            generator=lambda *_args: _success("CCO"),
            resolver=lambda *_args, **_kwargs: calls.append(True),
            method=method,
        )
        assert outcome.status is MolecularAssistantStatus.SUCCESS
        assert calls == []
        if method is MolecularResolutionMethod.AI_REFERENCE:
            assert outcome.identity.status is MolecularIdentityStatus.NOT_APPLICABLE
        else:
            assert outcome.identity.status is MolecularIdentityStatus.UNVERIFIED
            assert outcome.origin is StructureOrigin.AI


def test_reference_only_skips_model_and_does_not_require_provider_configuration():
    generator_calls = []
    resolver_calls = []
    outcome = _run_job(
        prompt="Dibuja la tetrandrina",
        generator=lambda *_args: generator_calls.append(True),
        resolver=lambda name, *, allow_network: resolver_calls.append((name, allow_network))
        or _reference(name, "CCO"),
        method=MolecularResolutionMethod.REFERENCE,
        provider_id="invalid-provider-id-is-ignored",
    )

    assert generator_calls == []
    assert resolver_calls == [("tetrandrina", False)]
    assert outcome.status is MolecularAssistantStatus.SUCCESS
    assert outcome.default_candidate == "reference"
    assert outcome.origin is StructureOrigin.REFERENCE
    assert outcome.provider_id == "reference"


def test_reference_only_requires_explicit_molecule_name_without_lookup_or_model():
    calls = []
    outcome = _run_job(
        prompt="Genera una molécula con tres anillos",
        generator=lambda *_args: calls.append("model"),
        resolver=lambda *_args, **_kwargs: calls.append("reference"),
        method=MolecularResolutionMethod.REFERENCE,
    )

    assert calls == []
    assert outcome.status is MolecularAssistantStatus.PROVIDER_ERROR
    assert outcome.reason_code == "reference_name_required"
    assert outcome.graph is None


def test_invalid_reference_smiles_is_rejected_before_it_can_be_inserted():
    outcome = _run_job(
        prompt="Draw caffeine",
        generator=lambda *_args: pytest.fail("reference-only route must not call AI"),
        resolver=lambda name, *, allow_network: NameToStructureResult(
            name,
            _graph("C"),
            "pubchem",
            0.95,
            smiles="not a valid SMILES",
            resolved_name=name,
        ),
        method=MolecularResolutionMethod.REFERENCE,
    )

    assert outcome.status is MolecularAssistantStatus.PROVIDER_ERROR
    assert outcome.reason_code == "reference_invalid"
    assert outcome.reference_result.graph is None
    assert outcome.graph is None


def test_transform_request_remains_ai_only_even_when_its_text_names_a_molecule():
    reference_calls = []
    observed = []
    source_graph = _graph("C")
    outcome = _run_job(
        prompt="Replace the selected molecule with ethanol",
        generator=lambda request, _config: observed.append(request)
        or _success("CCO"),
        resolver=lambda *_args, **_kwargs: reference_calls.append(True),
        method=MolecularResolutionMethod.AI,
        source_graph=source_graph,
    )

    assert outcome.status is MolecularAssistantStatus.SUCCESS
    assert isinstance(observed[0], MolecularTransformationRequest)
    assert reference_calls == []


def test_abandoned_ai_job_does_not_start_a_later_reference_lookup():
    app = QApplication.instance()
    assert app is not None
    started = threading.Event()
    release = threading.Event()
    reference_calls = []

    def generator(*_args):
        started.set()
        release.wait(timeout=3.0)
        return _success("CCO")

    controller = MolecularAssistantController(
        generator=generator,
        reference_resolver=lambda name, *, allow_network: reference_calls.append(name)
        or _reference(name, "CCO"),
    )
    finished = QSignalSpy(controller.job_finished)
    job_id = controller.start_job(
        "Draw ethanol",
        base_url="https://provider.example/v1",
        model="offline-test-model",
        resolution_method=MolecularResolutionMethod.AI_REFERENCE,
    )
    assert job_id is not None and started.wait(timeout=2.0)
    controller.abandon_job(job_id)
    release.set()
    deadline = time.monotonic() + 4.0
    while controller.active_jobs() and time.monotonic() < deadline:
        app.processEvents()
        time.sleep(0.005)

    assert controller.active_jobs() == ()
    assert reference_calls == []
    assert len(finished) == 0
    controller.deleteLater()
    app.processEvents()


def test_reference_result_graph_is_revalidated_in_chemio_and_origin_is_closed_enum():
    outcome = _run_job(
        prompt="Draw ethanol",
        generator=lambda *_args: _failure(
            "timeout", MolecularAssistantStatus.PROVIDER_ERROR
        ),
        resolver=lambda name, *, allow_network: _reference(name, "CCO"),
    )

    assert isinstance(outcome.origin, StructureOrigin)
    assert outcome.origin is StructureOrigin.REFERENCE
    assert outcome.reference_result.graph is not None
    assert len(outcome.reference_result.graph.atoms) == len(_graph("CCO").atoms)
