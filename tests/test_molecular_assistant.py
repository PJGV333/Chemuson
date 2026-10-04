from __future__ import annotations

import json
import socket
import urllib.error
from dataclasses import dataclass

import pytest

from chemuson.chemio.rdkit_io import smiles_to_molgraph
from chemuson.chemio.rdkit_safe import smiles_to_molgraph_isolated
from chemuson.core.model import MolGraph
from chemuson.molecular_assistant import (
    MolecularAssistant,
    MolecularAssistantRequest,
    MolecularAssistantStatus,
    OpenAICompatibleConfig,
    OpenAICompatibleProvider,
    ProviderCancelled,
    ProviderError,
    ProviderErrorCode,
    ProviderResponse,
)
from chemuson.molecular_assistant import provider as provider_module
from chemuson.molecular_assistant import service as service_module
from chemuson.molecular_assistant.limits import (
    CHEMIO_VALIDATION_TIMEOUT_S,
    MAX_HTTP_RESPONSE_BYTES,
    MAX_MODEL_CONTENT_BYTES,
    MAX_PROMPT_BYTES,
)


@dataclass
class FakeProvider:
    content: str = '{"smiles":"CCO"}'
    provider_id: str = "test-provider"
    model_id: str | None = "test-model"
    failure: Exception | None = None

    def __post_init__(self):
        self.calls: list[MolecularAssistantRequest] = []

    def generate(self, request: MolecularAssistantRequest) -> ProviderResponse:
        self.calls.append(request)
        if self.failure is not None:
            raise self.failure
        return ProviderResponse(content=self.content, model_id=self.model_id)


@dataclass
class FakeTransport:
    response: provider_module.HttpResponse | None = None
    failure: Exception | None = None

    def __post_init__(self):
        self.calls: list[dict[str, object]] = []

    def post(self, url, *, headers, body, timeout_s, max_response_bytes):
        self.calls.append(
            {
                "url": url,
                "headers": dict(headers),
                "body": body,
                "timeout_s": timeout_s,
                "max_response_bytes": max_response_bytes,
            }
        )
        if self.failure is not None:
            raise self.failure
        assert self.response is not None
        return self.response


def _openai_body(content: str, *, model: str = "served-model") -> bytes:
    return json.dumps(
        {"choices": [{"message": {"content": content}}], "model": model}
    ).encode("utf-8")


def _single_atom_graph() -> MolGraph:
    graph = MolGraph()
    graph.add_atom("C", 0.0, 0.0)
    return graph


def _graph_chemistry_signature(graph: MolGraph):
    atoms = tuple(
        sorted(
            (
                atom.id,
                atom.element,
                atom.formal_charge,
                atom.isotope,
                atom.radical_electrons,
                atom.stereo_cip,
            )
            for atom in graph.atoms.values()
        )
    )
    bonds = tuple(
        sorted(
            (
                bond.a1_id,
                bond.a2_id,
                bond.order,
                bond.is_aromatic,
                bond.stereo.value,
                bond.stereo_ez,
            )
            for bond in graph.bonds.values()
        )
    )
    return atoms, bonds


def test_valid_fake_provider_uses_isolated_chemio_and_matches_standard_importer():
    parsed_graph, error = smiles_to_molgraph_isolated("CCO")
    if parsed_graph is None:
        pytest.skip(f"ChemIO RDKit worker unavailable: {error}")
    provider = FakeProvider()

    result = MolecularAssistant(provider).generate(
        MolecularAssistantRequest("Draw ethanol")
    )

    assert result.status is MolecularAssistantStatus.SUCCESS
    assert result.graph is not None and result.graph.atoms
    assert result.validation_passed is True
    assert result.reason_code is None
    assert result.proposed_smiles == "CCO"
    assert result.provider_id == "test-provider"
    assert result.model_id == "test-model"
    assert len(provider.calls) == 1
    assert _graph_chemistry_signature(result.graph) == _graph_chemistry_signature(
        smiles_to_molgraph("CCO")
    )


def test_real_invalid_smiles_is_rejected_by_isolated_chemio_when_available():
    _graph, error = smiles_to_molgraph_isolated("this is not a SMILES")
    if error == "rdkit_unavailable":
        pytest.skip("ChemIO RDKit worker unavailable")

    result = MolecularAssistant(FakeProvider(content='{"smiles":"this is not a SMILES"}')).generate(
        MolecularAssistantRequest("invalid structure")
    )

    assert result.status is MolecularAssistantStatus.INVALID_STRUCTURE
    assert result.reason_code == "invalid_smiles"
    assert result.graph is None
    assert result.validation_passed is False


def test_parser_receives_exact_smiles_and_configured_isolated_timeout(monkeypatch):
    provider = FakeProvider(content='{"smiles":" CCO "}')
    observed = {}
    graph = MolGraph()
    graph.add_atom("C", 0.0, 0.0)

    def fake_parser(smiles, *, timeout_s):
        observed["smiles"] = smiles
        observed["timeout_s"] = timeout_s
        return graph, None

    monkeypatch.setattr(service_module, "smiles_to_molgraph_isolated", fake_parser)
    result = MolecularAssistant(provider).generate(MolecularAssistantRequest("ethanol"))

    assert result.status is MolecularAssistantStatus.SUCCESS
    assert result.proposed_smiles == " CCO "
    assert observed == {"smiles": " CCO ", "timeout_s": CHEMIO_VALIDATION_TIMEOUT_S}


@pytest.mark.parametrize(
    ("parser_error", "status", "reason", "validation_passed"),
    [
        ("invalid_input", MolecularAssistantStatus.INVALID_STRUCTURE, "invalid_smiles", False),
        ("timeout", MolecularAssistantStatus.VALIDATION_ERROR, "parser_timeout", None),
        ("rdkit_unavailable", MolecularAssistantStatus.VALIDATION_ERROR, "parser_unavailable", None),
        ("invalid_worker_json", MolecularAssistantStatus.VALIDATION_ERROR, "parser_error", None),
        ("worker_exit_code:2", MolecularAssistantStatus.VALIDATION_ERROR, "parser_error", None),
        ("parser exception says invalid_input", MolecularAssistantStatus.VALIDATION_ERROR, "parser_error", None),
        ("unexpected worker detail", MolecularAssistantStatus.VALIDATION_ERROR, "parser_error", None),
    ],
)
def test_chemio_error_identifiers_map_exactly_and_fail_closed(
    monkeypatch, parser_error, status, reason, validation_passed
):
    monkeypatch.setattr(
        service_module,
        "smiles_to_molgraph_isolated",
        lambda _smiles, *, timeout_s: (None, parser_error),
    )
    result = MolecularAssistant(FakeProvider()).generate(MolecularAssistantRequest("ethanol"))

    assert result.status is status
    assert result.reason_code == reason
    assert result.validation_passed is validation_passed
    assert result.proposed_smiles == "CCO"
    assert result.graph is None


def test_fake_provider_response_parsing_is_deterministic(monkeypatch):
    monkeypatch.setattr(
        service_module,
        "smiles_to_molgraph_isolated",
        lambda _smiles, *, timeout_s: (_single_atom_graph(), None),
    )
    provider = FakeProvider()
    assistant = MolecularAssistant(provider)

    first = assistant.generate(MolecularAssistantRequest("ethanol"))
    second = assistant.generate(MolecularAssistantRequest("ethanol"))

    assert first.status is second.status is MolecularAssistantStatus.SUCCESS
    assert first.reason_code == second.reason_code is None
    assert first.proposed_smiles == second.proposed_smiles == "CCO"
    assert _graph_chemistry_signature(first.graph) == _graph_chemistry_signature(second.graph)
    assert len(provider.calls) == 2


def test_missing_or_unknown_parser_error_is_generic(monkeypatch):
    monkeypatch.setattr(
        service_module,
        "smiles_to_molgraph_isolated",
        lambda _smiles, *, timeout_s: (None, None),
    )
    result = MolecularAssistant(FakeProvider()).generate(MolecularAssistantRequest("ethanol"))

    assert result.status is MolecularAssistantStatus.VALIDATION_ERROR
    assert result.reason_code == "parser_error"
    assert result.graph is None
    assert result.validation_passed is None


def test_parser_exception_text_is_not_exposed(monkeypatch):
    def broken_parser(_smiles, *, timeout_s):
        raise RuntimeError("invalid_input: secret parser diagnostic")

    monkeypatch.setattr(service_module, "smiles_to_molgraph_isolated", broken_parser)
    result = MolecularAssistant(FakeProvider()).generate(MolecularAssistantRequest("ethanol"))

    assert result.status is MolecularAssistantStatus.VALIDATION_ERROR
    assert result.reason_code == "parser_error"
    assert result.proposed_smiles == "CCO"
    assert result.graph is None


@pytest.mark.parametrize(
    ("content", "reason"),
    [
        ("", "invalid_json"),
        ("not json", "invalid_json"),
        ('{"smiles":"CCO"', "invalid_json"),
        ('{"smiles":"CCO","name":"ethanol"}', "unexpected_fields"),
        ('{"name":"ethanol"}', "missing_smiles"),
        ('{"smiles":"   "}', "empty_smiles"),
        ('{"smiles":4}', "invalid_json"),
        ('["CCO"]', "invalid_json"),
        ('```json\n{"smiles":"CCO"}\n```', "invalid_json"),
        ('prefix {"smiles":"CCO"}', "invalid_json"),
        ('{"smiles":"CCO"} suffix', "invalid_json"),
        ('{"smiles":"CCO","smiles":"CC"}', "invalid_json"),
        ('{"smiles":NaN}', "invalid_json"),
    ],
)
def test_response_decoder_rejects_non_contract_content_without_chemio(
    monkeypatch, content, reason
):
    def parser_must_not_run(*_args, **_kwargs):
        pytest.fail("ChemIO must not receive malformed model output")

    monkeypatch.setattr(service_module, "smiles_to_molgraph_isolated", parser_must_not_run)
    result = MolecularAssistant(FakeProvider(content=content)).generate(
        MolecularAssistantRequest("ethanol")
    )

    assert result.status is MolecularAssistantStatus.MALFORMED_RESPONSE
    assert result.reason_code == reason
    assert result.graph is None
    assert result.validation_passed is None


def test_model_content_limit_is_checked_before_json_parsing(monkeypatch):
    content = '{"smiles":"' + ("C" * MAX_MODEL_CONTENT_BYTES) + '"}'
    monkeypatch.setattr(
        service_module,
        "smiles_to_molgraph_isolated",
        lambda *_args, **_kwargs: pytest.fail("oversized response must not reach ChemIO"),
    )
    result = MolecularAssistant(FakeProvider(content=content)).generate(
        MolecularAssistantRequest("ethanol")
    )
    assert result.status is MolecularAssistantStatus.MALFORMED_RESPONSE
    assert result.reason_code == "response_too_large"
    assert result.graph is None


@pytest.mark.parametrize(
    ("description", "reason"),
    [(" \n\t", "empty_prompt"), ("x" * (MAX_PROMPT_BYTES + 1), "request_too_large")],
)
def test_invalid_prompt_is_rejected_before_provider_call(description, reason):
    provider = FakeProvider()

    result = MolecularAssistant(provider).generate(MolecularAssistantRequest(description))

    assert result.status is MolecularAssistantStatus.INVALID_REQUEST
    assert result.reason_code == reason
    assert result.graph is None
    assert provider.calls == []


def test_wrong_request_type_and_invalid_unicode_are_controlled():
    provider = FakeProvider()
    assistant = MolecularAssistant(provider)

    wrong_type = assistant.generate("draw ethanol")  # type: ignore[arg-type]
    invalid_unicode = assistant.generate(MolecularAssistantRequest("\ud800"))

    assert wrong_type.status is MolecularAssistantStatus.INVALID_REQUEST
    assert wrong_type.reason_code == "invalid_prompt"
    assert invalid_unicode.status is MolecularAssistantStatus.INVALID_REQUEST
    assert invalid_unicode.reason_code == "invalid_prompt"
    assert provider.calls == []


@pytest.mark.parametrize(
    ("failure", "reason"),
    [
        (ProviderError(ProviderErrorCode.TIMEOUT), "timeout"),
        (ProviderError(ProviderErrorCode.NETWORK_ERROR), "network_error"),
        (ProviderError(ProviderErrorCode.HTTP_ERROR), "http_error"),
        (ProviderError(ProviderErrorCode.PROVIDER_ERROR), "provider_error"),
        (RuntimeError("secret transport diagnostic"), "provider_error"),
    ],
)
def test_provider_failures_are_sanitized(failure, reason):
    result = MolecularAssistant(FakeProvider(failure=failure)).generate(
        MolecularAssistantRequest("ethanol")
    )

    assert result.status is MolecularAssistantStatus.PROVIDER_ERROR
    assert result.reason_code == reason
    assert result.graph is None
    assert result.validation_passed is None
    assert "secret" not in repr(result)


def test_explicit_cancellation_returns_no_partial_graph():
    result = MolecularAssistant(FakeProvider(failure=ProviderCancelled())).generate(
        MolecularAssistantRequest("ethanol")
    )

    assert result.status is MolecularAssistantStatus.CANCELLED
    assert result.reason_code == "cancelled"
    assert result.graph is None
    assert result.validation_passed is None


def test_fake_provider_path_does_not_open_network(monkeypatch):
    def network_must_not_open(*_args, **_kwargs):
        pytest.fail("fake-provider tests must not open a network connection")

    monkeypatch.setattr(provider_module.urllib.request, "build_opener", network_must_not_open)
    monkeypatch.setattr(
        service_module,
        "smiles_to_molgraph_isolated",
        lambda _smiles, *, timeout_s: (_single_atom_graph(), None),
    )
    result = MolecularAssistant(FakeProvider()).generate(MolecularAssistantRequest("ethanol"))

    assert result.status is MolecularAssistantStatus.SUCCESS


def test_openai_adapter_posts_explicit_request_and_extracts_transport_metadata():
    transport = FakeTransport(
        response=provider_module.HttpResponse(200, _openai_body('{"smiles":"CCO"}'))
    )
    config = OpenAICompatibleConfig(
        base_url="https://provider.example/v1/",
        model="configured-model",
        api_key="test-secret",
        timeout_s=12.5,
        supports_json_output=True,
        provider_id="test-gateway",
    )
    provider = OpenAICompatibleProvider(config, transport=transport)

    response = provider.generate(MolecularAssistantRequest("Draw ethanol"))

    call = transport.calls[0]
    payload = json.loads(call["body"])
    assert call["url"] == "https://provider.example/v1/chat/completions"
    assert call["headers"]["Authorization"] == "Bearer test-secret"
    assert call["headers"]["Content-Type"] == "application/json"
    assert call["timeout_s"] == 12.5
    assert call["max_response_bytes"] == MAX_HTTP_RESPONSE_BYTES
    assert payload["model"] == "configured-model"
    assert payload["stream"] is False
    assert payload["response_format"] == {"type": "json_object"}
    assert payload["messages"][1] == {"role": "user", "content": "Draw ethanol"}
    assert response == ProviderResponse('{"smiles":"CCO"}', "served-model")
    assert provider.provider_id == "test-gateway"
    assert provider.model_id == "configured-model"
    assert "test-secret" not in repr(config)


def test_openai_adapter_omits_optional_json_mode_by_default():
    transport = FakeTransport(
        response=provider_module.HttpResponse(200, _openai_body('{"smiles":"CCO"}'))
    )
    provider = OpenAICompatibleProvider(
        OpenAICompatibleConfig("http://localhost:8000", "local-model"),
        transport=transport,
    )

    provider.generate(MolecularAssistantRequest("Draw ethanol"))

    payload = json.loads(transport.calls[0]["body"])
    assert "response_format" not in payload
    assert "Authorization" not in transport.calls[0]["headers"]


@pytest.mark.parametrize(
    ("failure", "reason"),
    [
        (TimeoutError("secret"), "timeout"),
        (socket.timeout("secret"), "timeout"),
        (OSError("secret"), "network_error"),
        (urllib.error.URLError("secret"), "network_error"),
        (urllib.error.URLError(TimeoutError("secret")), "timeout"),
    ],
)
def test_http_transport_errors_map_to_stable_codes(failure, reason):
    provider = OpenAICompatibleProvider(
        OpenAICompatibleConfig("https://provider.example", "model"),
        transport=FakeTransport(failure=failure),
    )
    result = MolecularAssistant(provider).generate(MolecularAssistantRequest("ethanol"))

    assert result.status is MolecularAssistantStatus.PROVIDER_ERROR
    assert result.reason_code == reason
    assert result.graph is None
    assert "secret" not in repr(result)


def test_http_status_and_invalid_provider_envelope_are_controlled():
    http_error_provider = OpenAICompatibleProvider(
        OpenAICompatibleConfig("https://provider.example", "model"),
        transport=FakeTransport(response=provider_module.HttpResponse(429, b"secret")),
    )
    protocol_error_provider = OpenAICompatibleProvider(
        OpenAICompatibleConfig("https://provider.example", "model"),
        transport=FakeTransport(response=provider_module.HttpResponse(200, b"{}")),
    )

    http_result = MolecularAssistant(http_error_provider).generate(
        MolecularAssistantRequest("ethanol")
    )
    protocol_result = MolecularAssistant(protocol_error_provider).generate(
        MolecularAssistantRequest("ethanol")
    )

    assert http_result.status is MolecularAssistantStatus.PROVIDER_ERROR
    assert http_result.reason_code == "http_error"
    assert protocol_result.status is MolecularAssistantStatus.PROVIDER_ERROR
    assert protocol_result.reason_code == "provider_error"
    assert "secret" not in repr(http_result)


def test_oversized_http_body_maps_to_malformed_response():
    provider = OpenAICompatibleProvider(
        OpenAICompatibleConfig("https://provider.example", "model"),
        transport=FakeTransport(
            response=provider_module.HttpResponse(
                200,
                b" " * (MAX_HTTP_RESPONSE_BYTES + 1),
            )
        ),
    )
    result = MolecularAssistant(provider).generate(MolecularAssistantRequest("ethanol"))

    assert result.status is MolecularAssistantStatus.MALFORMED_RESPONSE
    assert result.reason_code == "response_too_large"
    assert result.graph is None


def test_real_urllib_transport_reads_at_most_limit_plus_one(monkeypatch):
    class FakeResponse:
        status = 200

        def __init__(self):
            self.requested_read_size = None

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return None

        def read(self, size):
            self.requested_read_size = size
            return b"x" * size

    class FakeOpener:
        def __init__(self, response):
            self.response = response

        def open(self, _request, *, timeout):
            assert timeout == 1.0
            return self.response

    response = FakeResponse()
    monkeypatch.setattr(
        provider_module.urllib.request,
        "build_opener",
        lambda *_handlers: FakeOpener(response),
    )
    with pytest.raises(ProviderError) as error:
        provider_module._UrllibTransport().post(
            "https://provider.example/v1/chat/completions",
            headers={},
            body=b"{}",
            timeout_s=1.0,
            max_response_bytes=10,
        )

    assert error.value.code is ProviderErrorCode.RESPONSE_TOO_LARGE
    assert response.requested_read_size == 11


def test_redirect_handler_does_not_follow_server_selected_destination():
    handler = provider_module._NoRedirectHandler()
    assert handler.redirect_request(None, None, 302, "Found", {}, "http://elsewhere/") is None


@pytest.mark.parametrize(
    "config_kwargs",
    [
        {"base_url": "", "model": "model"},
        {"base_url": "file:///tmp/model", "model": "model"},
        {"base_url": "https://user:pass@provider.example", "model": "model"},
        {"base_url": "https://provider.example?token=secret", "model": "model"},
        {"base_url": " https://provider.example", "model": "model"},
        {"base_url": "https://provider.example", "model": " "},
        {"base_url": "https://provider.example", "model": " model"},
        {"base_url": "https://provider.example", "model": "model", "timeout_s": 0},
        {"base_url": "https://provider.example", "model": "model", "timeout_s": True},
        {"base_url": "https://provider.example", "model": "model", "timeout_s": float("inf")},
        {"base_url": "https://provider.example", "model": "model", "timeout_s": 10**1000},
    ],
)
def test_provider_configuration_requires_explicit_safe_endpoint_and_finite_timeout(
    config_kwargs,
):
    with pytest.raises(ValueError):
        OpenAICompatibleConfig(**config_kwargs)


@pytest.mark.parametrize("api_key", ["first\nsecond", "first\rsecond", "first\r\nsecond"])
def test_api_key_rejects_actual_line_break_characters(api_key):
    with pytest.raises(ValueError, match="api_key"):
        OpenAICompatibleConfig("https://provider.example", "model", api_key=api_key)


@pytest.mark.parametrize("api_key", ["first\\nsecond", "first\\rsecond"])
def test_api_key_allows_literal_backslash_sequences(api_key):
    config = OpenAICompatibleConfig("https://provider.example", "model", api_key=api_key)
    assert config.api_key == api_key
