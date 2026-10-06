# Design

## Decisions

1. Add a provider-neutral capability enum: `PROMPT_ONLY`, `OPENAI_JSON_OBJECT`, `UNKNOWN`. Preserve the existing `supports_json_output` argument as a compatibility shim (`false`→prompt-only, `true`→native JSON object); callers can use the enum when they need the unknown state.
2. Native JSON mode uses `response_format={"type":"json_object"}`. Recognize unsupported capability only from HTTP 400 with the exact JSON error code `response_format_not_supported`. Then issue one text-mode fallback using the same strict system contract. Generic HTTP/auth/rate-limit errors never trigger fallback. Cache the observed capability only on the provider instance.
3. Keep response decoding in M23 strict. Only `invalid_json` is eligible for one format-repair call. Pass a size-bounded, JSON-quoted copy of the original `message.content` as untrusted data; do not extract substrings, strip Markdown, repair JSON locally, or accept a different schema. Run the same strict decoder and isolated ChemIO validation on the retry. Generation and transform share this `generate` path.
4. Return only bounded boolean diagnostics (structured requested/native/fallback and repair used/succeeded) with operation results. Do not include or persist original model content, reasoning fields, credentials, headers, or exception text. The GUI may show a small diagnostic summary only.
5. Identity verification remains the existing downstream GUI stage after any M23 success; a repaired proposal receives no special trust or insertion path.

## Bounds and failure semantics

Format repair is limited to one retry, with a 4096-byte original-content cap and the existing 8 KiB prompt bound. A second invalid response remains `malformed_response/invalid_json`; no third request occurs. Oversized or non-`invalid_json` failures are not repair candidates. ChemIO remains authoritative and can still return `invalid_structure`.

## Files and verification

Expected changes: `molecular_assistant/{provider,models,service,limits,__init__}.py`, the existing dialog/controller presentation path if needed, focused tests, and this OpenSpec delta. Use fake transports/providers only for contract tests; run a few local Qwen requests only if the configured loopback server is already running. Never start a server or use external Internet.
