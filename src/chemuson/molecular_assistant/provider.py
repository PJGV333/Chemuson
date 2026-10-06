"""OpenAI-compatible transport adapter with injectable, bounded HTTP I/O."""

from __future__ import annotations

import ipaddress
import json
import math
import socket
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from enum import Enum
from typing import Mapping, Protocol
from urllib.parse import urlsplit

from chemuson.molecular_assistant.limits import (
    DEFAULT_MAX_OUTPUT_TOKENS,
    DEFAULT_PROVIDER_TIMEOUT_S,
    MAX_HTTP_RESPONSE_BYTES,
    MAX_MAX_OUTPUT_TOKENS,
    MAX_PROVIDER_TIMEOUT_S,
    MIN_MAX_OUTPUT_TOKENS,
    MIN_PROVIDER_TIMEOUT_S,
)
from chemuson.molecular_assistant.models import (
    MolecularAssistantRequest,
    ProviderResponse,
)


class ProviderErrorCode(str, Enum):
    """Finite provider failure vocabulary; never contains raw diagnostics."""

    TIMEOUT = "timeout"
    NETWORK_ERROR = "network_error"
    HTTP_ERROR = "http_error"
    PROVIDER_ERROR = "provider_error"
    RESPONSE_TOO_LARGE = "response_too_large"


class ProviderError(Exception):
    """A provider failure carrying only a stable, allowlisted code."""

    def __init__(self, code: ProviderErrorCode) -> None:
        self.code = ProviderErrorCode(code)
        super().__init__(self.code.value)


class ProviderCancelled(Exception):
    """Raised by a provider when an operation is explicitly cancelled."""


class MolecularStructureProvider(Protocol):
    """Replaceable provider contract consumed by the application service."""

    provider_id: str
    model_id: str | None

    def generate(self, request: MolecularAssistantRequest) -> ProviderResponse:
        """Return raw structured-response content without parsing chemistry."""


@dataclass(frozen=True)
class OpenAICompatibleConfig:
    """Explicit endpoint/model configuration for Chat Completions."""

    base_url: str
    model: str
    api_key: str | None = field(default=None, repr=False)
    timeout_s: float = DEFAULT_PROVIDER_TIMEOUT_S
    supports_json_output: bool = False
    provider_id: str = "openai-compatible"
    max_tokens: int = DEFAULT_MAX_OUTPUT_TOKENS

    def __post_init__(self) -> None:
        if (
            not isinstance(self.base_url, str)
            or not self.base_url.strip()
            or self.base_url != self.base_url.strip()
        ):
            raise ValueError("base_url must be an explicit HTTP(S) URL")
        try:
            parsed = urlsplit(self.base_url)
            _ = parsed.port
        except ValueError as exc:
            raise ValueError("base_url is invalid") from exc
        if (
            parsed.scheme not in {"http", "https"}
            or not parsed.hostname
            or parsed.username is not None
            or parsed.password is not None
            or parsed.query
            or parsed.fragment
        ):
            raise ValueError("base_url must be an explicit HTTP(S) URL without credentials")
        if (
            not isinstance(self.model, str)
            or not self.model.strip()
            or self.model != self.model.strip()
        ):
            raise ValueError("model must be explicitly configured")
        try:
            valid_timeout = (
                not isinstance(self.timeout_s, bool)
                and isinstance(self.timeout_s, (int, float))
                and math.isfinite(float(self.timeout_s))
                and self.timeout_s > 0
            )
        except (OverflowError, TypeError, ValueError):
            valid_timeout = False
        if not valid_timeout or not MIN_PROVIDER_TIMEOUT_S <= float(self.timeout_s) <= MAX_PROVIDER_TIMEOUT_S:
            raise ValueError("timeout_s must be finite and between 10 and 600 seconds")
        if self.api_key is not None and (
            not isinstance(self.api_key, str)
            or "\r" in self.api_key
            or "\n" in self.api_key
        ):
            raise ValueError("api_key must be an opaque single-line string")
        if (
            self.api_key
            and parsed.scheme == "http"
            and not _is_loopback_host(parsed.hostname or "")
        ):
            raise ValueError(
                "API keys require HTTPS or a loopback HTTP endpoint"
            )
        if not isinstance(self.supports_json_output, bool):
            raise ValueError("supports_json_output must be boolean")
        if (
            isinstance(self.max_tokens, bool)
            or not isinstance(self.max_tokens, int)
            or not MIN_MAX_OUTPUT_TOKENS <= self.max_tokens <= MAX_MAX_OUTPUT_TOKENS
        ):
            raise ValueError("max_tokens must be an integer between 64 and 8192")
        if (
            not isinstance(self.provider_id, str)
            or not self.provider_id.strip()
            or self.provider_id != self.provider_id.strip()
        ):
            raise ValueError("provider_id must be non-empty")

    @property
    def endpoint_url(self) -> str:
        """Build the fixed Chat Completions path from the configured base URL."""
        base = self.base_url.rstrip("/")
        if urlsplit(base).path.rstrip("/").endswith("/v1"):
            return f"{base}/chat/completions"
        return f"{base}/v1/chat/completions"


def _is_loopback_host(host: str) -> bool:
    """Recognize only literal loopback IPs and the exact localhost name."""
    normalized = host.casefold()
    if normalized == "localhost":
        return True
    if "%" in normalized:
        return False
    try:
        return ipaddress.ip_address(normalized).is_loopback
    except ValueError:
        return False


@dataclass(frozen=True)
class HttpResponse:
    """Minimal response shape returned by an HTTP transport."""

    status_code: int
    body: bytes


class HttpTransport(Protocol):
    """Small injectable POST interface for deterministic offline tests."""

    def post(
        self,
        url: str,
        *,
        headers: Mapping[str, str],
        body: bytes,
        timeout_s: float,
        max_response_bytes: int,
    ) -> HttpResponse:
        """POST bytes and return at most the configured response limit."""


class _NoRedirectHandler(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


class _UrllibTransport:
    """Standard-library POST implementation that reads at most limit + 1."""

    def post(
        self,
        url: str,
        *,
        headers: Mapping[str, str],
        body: bytes,
        timeout_s: float,
        max_response_bytes: int,
    ) -> HttpResponse:
        request = urllib.request.Request(
            url,
            data=body,
            headers=dict(headers),
            method="POST",
        )
        opener = urllib.request.build_opener(_NoRedirectHandler())
        with opener.open(request, timeout=timeout_s) as response:
            response_body = response.read(max_response_bytes + 1)
        if len(response_body) > max_response_bytes:
            raise ProviderError(ProviderErrorCode.RESPONSE_TOO_LARGE) from None
        return HttpResponse(status_code=response.status, body=response_body)


class OpenAICompatibleProvider:
    """Non-streaming `/v1/chat/completions` adapter; no SDK or implicit endpoint."""

    def __init__(
        self,
        config: OpenAICompatibleConfig,
        *,
        transport: HttpTransport | None = None,
    ) -> None:
        self.config = config
        self.provider_id = config.provider_id
        self.model_id = config.model
        self._transport = transport if transport is not None else _UrllibTransport()

    def generate(self, request: MolecularAssistantRequest) -> ProviderResponse:
        """Request one JSON SMILES proposal and return only message content."""
        payload: dict[str, object] = {
            "model": self.config.model,
            "messages": [
                {
                    "role": "system",
                    "content": (
                        'Return exactly one JSON object with the single string field '
                        '"smiles". Do not include Markdown, explanations, or extra fields.'
                    ),
                },
                {"role": "user", "content": request.description},
            ],
            "stream": False,
            "max_tokens": self.config.max_tokens,
        }
        if self.config.supports_json_output:
            payload["response_format"] = {"type": "json_object"}

        headers = {
            "Accept": "application/json",
            "Content-Type": "application/json",
        }
        if self.config.api_key:
            headers["Authorization"] = f"Bearer {self.config.api_key}"
        request_body = json.dumps(payload, ensure_ascii=False, separators=(",", ":")).encode("utf-8")

        try:
            response = self._transport.post(
                self.config.endpoint_url,
                headers=headers,
                body=request_body,
                timeout_s=float(self.config.timeout_s),
                max_response_bytes=MAX_HTTP_RESPONSE_BYTES,
            )
        except ProviderError:
            raise
        except urllib.error.HTTPError:
            raise ProviderError(ProviderErrorCode.HTTP_ERROR) from None
        except (TimeoutError, socket.timeout):
            raise ProviderError(ProviderErrorCode.TIMEOUT) from None
        except urllib.error.URLError as exc:
            if isinstance(exc.reason, (TimeoutError, socket.timeout)):
                raise ProviderError(ProviderErrorCode.TIMEOUT) from None
            raise ProviderError(ProviderErrorCode.NETWORK_ERROR) from None
        except OSError:
            raise ProviderError(ProviderErrorCode.NETWORK_ERROR) from None
        except Exception:
            raise ProviderError(ProviderErrorCode.PROVIDER_ERROR) from None

        if (
            not isinstance(response, HttpResponse)
            or isinstance(response.status_code, bool)
            or not isinstance(response.status_code, int)
            or not isinstance(response.body, bytes)
        ):
            raise ProviderError(ProviderErrorCode.PROVIDER_ERROR) from None
        if not 200 <= response.status_code < 300:
            raise ProviderError(ProviderErrorCode.HTTP_ERROR) from None
        if len(response.body) > MAX_HTTP_RESPONSE_BYTES:
            raise ProviderError(ProviderErrorCode.RESPONSE_TOO_LARGE) from None

        try:
            envelope = json.loads(response.body.decode("utf-8"))
            content = envelope["choices"][0]["message"]["content"]
        except (KeyError, IndexError, TypeError, ValueError, UnicodeError, RecursionError):
            raise ProviderError(ProviderErrorCode.PROVIDER_ERROR) from None
        if not isinstance(envelope, dict) or not isinstance(content, str):
            raise ProviderError(ProviderErrorCode.PROVIDER_ERROR) from None

        response_model = envelope.get("model")
        model_id = response_model.strip() if isinstance(response_model, str) else ""
        return ProviderResponse(content=content, model_id=model_id or self.config.model)
