from __future__ import annotations

import io
import json
import urllib.error

import pytest

from chemuson.core.model import MolGraph
from chemuson.molecular_assistant import (
    MolecularAssistant,
    MolecularAssistantRequest,
    MolecularAssistantStatus,
    MolecularTransformationRequest,
    OpenAICompatibleConfig,
    OpenAICompatibleProvider,
    ProviderResponse,
    StructuredOutputCapability,
)
from chemuson.molecular_assistant import provider as provider_module
from chemuson.molecular_assistant import service as service_module
from chemuson.molecular_assistant.limits import MAX_FORMAT_REPAIR_CONTENT_BYTES


def _body(content: str, *, reasoning: str | None = None) -> bytes:
    message = {"content": content}
    if reasoning is not None:
        message["reasoning_content"] = reasoning
    return json.dumps(
        {"choices": [{"message": message}], "model": "served-model"}
    ).encode("utf-8")


def _unsupported_body() -> bytes:
    return json.dumps(
        {"error": {"code": "response_format_not_supported", "message": "ignored"}}
    ).encode("utf-8")


def _graph() -> MolGraph:
    graph = MolGraph()
    graph.add_atom("C", 0.0, 0.0)
    return graph


class QueueTransport:
    def __init__(self, responses: list[provider_module.HttpResponse | Exception]) -> None:
        self.responses = list(responses)
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
        assert self.responses, "unexpected extra HTTP request"
        response = self.responses.pop(0)
        if isinstance(response, Exception):
            raise response
        return response


class QueueProvider:
    provider_id = "queue-provider"
    model_id = "queue-model"

    def __init__(self, contents: list[str]) -> None:
        self.contents = list(contents)
        self.calls: list[MolecularAssistantRequest] = []

    def generate(self, request: MolecularAssistantRequest) -> ProviderResponse:
        self.calls.append(request)
        assert self.contents, "unexpected extra provider request"
        return ProviderResponse(self.contents.pop(0), self.model_id)


def _native_config() -> OpenAICompatibleConfig:
    return OpenAICompatibleConfig(
        "https://provider.example/v1",
        "configured-model",
        api_key="transient-test-key",
        supports_json_output=True,
    )


def test_explicitly_unsupported_native_json_falls_back_once_to_text_mode():
    transport = QueueTransport(
        [
            provider_module.HttpResponse(400, _unsupported_body()),
            provider_module.HttpResponse(200, _body('{"smiles":"CCO"}')),
        ]
    )
    provider = OpenAICompatibleProvider(_native_config(), transport=transport)

    response = provider.generate(MolecularAssistantRequest("Draw ethanol"))

    payloads = [json.loads(call["body"]) for call in transport.calls]
    assert len(payloads) == 2
    assert payloads[0]["response_format"] == {"type": "json_object"}
    assert "response_format" not in payloads[1]
    assert payloads[0]["messages"][0]["content"] == payloads[1]["messages"][0]["content"]
    assert response.structured_output_requested is True
    assert response.structured_output_native is False
    assert response.structured_output_fallback_used is True


def test_exact_unsupported_error_code_in_urllib_http_error_also_falls_back():
    error = urllib.error.HTTPError(
        "https://provider.example/v1/chat/completions",
        400,
        "Bad Request",
        None,
        io.BytesIO(_unsupported_body()),
    )
    transport = QueueTransport(
        [
            error,
            provider_module.HttpResponse(200, _body('{"smiles":"CCO"}')),
        ]
    )
    response = OpenAICompatibleProvider(_native_config(), transport=transport).generate(
        MolecularAssistantRequest("Draw ethanol")
    )
    assert len(transport.calls) == 2
    assert response.structured_output_fallback_used is True


@pytest.mark.parametrize(
    ("status", "error_body"),
    [
        (400, b'{"error":{"code":"invalid_request_error"}}'),
        (401, _unsupported_body()),
        (403, _unsupported_body()),
        (429, _unsupported_body()),
        (500, _unsupported_body()),
    ],
)
def test_only_http_400_with_exact_unsupported_code_triggers_fallback(status, error_body):
    transport = QueueTransport([provider_module.HttpResponse(status, error_body)])
    result = MolecularAssistant(OpenAICompatibleProvider(_native_config(), transport=transport)).generate(
        MolecularAssistantRequest("Draw ethanol")
    )
    assert result.status is MolecularAssistantStatus.PROVIDER_ERROR
    assert result.reason_code == "http_error"
    assert len(transport.calls) == 1
    assert result.structured_output_fallback_used is False


def test_native_output_prompt_backslash_escape_and_reasoning_contract():
    transport = QueueTransport(
        [
            provider_module.HttpResponse(
                200,
                _body('{"smiles":"CCO"}', reasoning="private reasoning that must not escape"),
            )
        ]
    )
    response = OpenAICompatibleProvider(_native_config(), transport=transport).generate(
        MolecularAssistantRequest("Draw ethanol")
    )
    system_prompt = json.loads(transport.calls[0]["body"])["messages"][0]["content"]

    assert r'{"smiles":"F/C=C\\F"}' in system_prompt
    assert "Do not reveal reasoning" in system_prompt
    assert "private reasoning" not in response.content
    assert response.structured_output_native is True


def test_strict_decoder_accepts_only_json_escaped_stereochemical_backslash():
    smiles = r"F/C=C\F"
    encoded = json.dumps({"smiles": smiles})

    assert service_module.MolecularAssistant._decode_smiles(encoded) == (smiles, None)
    assert service_module.MolecularAssistant._decode_smiles(
        r'{"smiles":"F/C=C\F"}'
    ) == (None, "invalid_json")


def test_invalid_json_gets_one_data_bounded_format_repair_then_chemio(monkeypatch):
    original = r'{"smiles":"F/C=C\F"}'
    smiles = r"F/C=C\F"
    provider = QueueProvider([original, json.dumps({"smiles": smiles})])
    observed: list[str] = []
    monkeypatch.setattr(
        service_module,
        "smiles_to_molgraph_isolated",
        lambda value, *, timeout_s: observed.append(value) or (_graph(), None),
    )

    result = MolecularAssistant(provider).generate(MolecularAssistantRequest("draw E-1,2-difluoroethene"))

    assert result.status is MolecularAssistantStatus.SUCCESS
    assert result.proposed_smiles == smiles
    assert observed == [smiles]
    assert len(provider.calls) == 2
    assert json.dumps(original, ensure_ascii=True, separators=(",", ":")) in provider.calls[1].description
    assert "untrusted data" in provider.calls[1].description
    assert result.format_repair_used is True
    assert result.format_repair_succeeded is True


def test_format_repair_is_shared_by_transform_path(monkeypatch):
    provider = QueueProvider(["not JSON", json.dumps({"smiles": "CCCl"})])
    monkeypatch.setattr(
        service_module,
        "smiles_to_molgraph_isolated",
        lambda _value, *, timeout_s: (_graph(), None),
    )

    result = MolecularAssistant(provider).transform(
        MolecularTransformationRequest("CCO", "replace oxygen with chlorine")
    )

    assert result.status is MolecularAssistantStatus.SUCCESS
    assert result.proposed_smiles == "CCCl"
    assert len(provider.calls) == 2
    assert "Source SMILES: CCO" in provider.calls[0].description
    assert result.format_repair_used is True
    assert result.format_repair_succeeded is True


def test_second_invalid_json_fails_without_a_third_attempt_or_chemio(monkeypatch):
    provider = QueueProvider(["prose", "still prose", json.dumps({"smiles": "CCO"})])
    monkeypatch.setattr(
        service_module,
        "smiles_to_molgraph_isolated",
        lambda *_args, **_kwargs: pytest.fail("invalid JSON must not reach ChemIO"),
    )

    result = MolecularAssistant(provider).generate(MolecularAssistantRequest("ethanol"))

    assert result.status is MolecularAssistantStatus.MALFORMED_RESPONSE
    assert result.reason_code == "invalid_json"
    assert len(provider.calls) == 2
    assert provider.contents == [json.dumps({"smiles": "CCO"})]
    assert result.format_repair_used is True
    assert result.format_repair_succeeded is False


def test_oversized_original_content_is_not_sent_to_repair(monkeypatch):
    provider = QueueProvider(["x" * (MAX_FORMAT_REPAIR_CONTENT_BYTES + 1), "unused"])
    monkeypatch.setattr(
        service_module,
        "smiles_to_molgraph_isolated",
        lambda *_args, **_kwargs: pytest.fail("malformed output must not reach ChemIO"),
    )

    result = MolecularAssistant(provider).generate(MolecularAssistantRequest("ethanol"))

    assert result.status is MolecularAssistantStatus.MALFORMED_RESPONSE
    assert result.reason_code == "invalid_json"
    assert len(provider.calls) == 1
    assert result.format_repair_used is False
    assert result.format_repair_succeeded is None


def test_repaired_json_still_requires_chemio_validation(monkeypatch):
    smiles = "not a valid SMILES"
    provider = QueueProvider(["bad JSON", json.dumps({"smiles": smiles})])
    chemio_calls: list[str] = []

    def reject_chemistry(value, *, timeout_s):
        chemio_calls.append(value)
        return None, "invalid_input"

    monkeypatch.setattr(service_module, "smiles_to_molgraph_isolated", reject_chemistry)
    result = MolecularAssistant(provider).generate(MolecularAssistantRequest("test"))

    assert result.status is MolecularAssistantStatus.INVALID_STRUCTURE
    assert result.reason_code == "invalid_smiles"
    assert chemio_calls == [smiles]
    assert result.format_repair_used is True
    assert result.format_repair_succeeded is True


def test_capability_fallback_and_format_repair_share_a_three_request_ceiling(monkeypatch):
    transport = QueueTransport(
        [
            provider_module.HttpResponse(400, _unsupported_body()),
            provider_module.HttpResponse(200, _body("not JSON", reasoning="private")),
            provider_module.HttpResponse(200, _body('{"smiles":"CCO"}')),
        ]
    )
    provider = OpenAICompatibleProvider(_native_config(), transport=transport)
    monkeypatch.setattr(
        service_module,
        "smiles_to_molgraph_isolated",
        lambda _value, *, timeout_s: (_graph(), None),
    )

    result = MolecularAssistant(provider).generate(MolecularAssistantRequest("ethanol"))

    payloads = [json.loads(call["body"]) for call in transport.calls]
    assert len(payloads) == 3
    assert "response_format" in payloads[0]
    assert all("response_format" not in payload for payload in payloads[1:])
    assert result.status is MolecularAssistantStatus.SUCCESS
    assert result.structured_output_requested is True
    assert result.structured_output_native is False
    assert result.structured_output_fallback_used is True
    assert result.format_repair_used is True
    assert result.format_repair_succeeded is True
    assert "transient-test-key" not in repr(result)
    assert "private" not in repr(result)


def test_unknown_capability_is_distinct_from_prompt_only():
    unknown = OpenAICompatibleConfig("https://provider.example", "model")
    prompt_only = OpenAICompatibleConfig(
        "https://provider.example",
        "model",
        structured_output_capability=StructuredOutputCapability.PROMPT_ONLY,
    )
    assert unknown.structured_output_capability is StructuredOutputCapability.UNKNOWN
    assert prompt_only.structured_output_capability is StructuredOutputCapability.PROMPT_ONLY
