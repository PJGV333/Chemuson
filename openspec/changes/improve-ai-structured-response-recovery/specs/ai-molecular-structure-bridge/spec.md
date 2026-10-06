## MODIFIED Requirements

### Requirement: Model output SHALL conform to a strict molecular response schema

The assistant SHALL accept exactly one required top-level JSON field, `smiles`, containing a non-empty string within the configured size limit. If the first provider response is not valid JSON, the assistant MAY issue exactly one bounded format-only repair request only when the strict decoder reports `invalid_json`; the repaired response SHALL pass the same exact decoder before ChemIO is called.

#### Scenario: Minimal valid structured proposal
- **GIVEN** provider content containing `{"smiles":"CCO"}`
- **WHEN** the response decoder parses it
- **THEN** the decoder SHALL produce the SMILES value without requiring a generated name or metadata
- **AND** no repair request SHALL be made.

#### Scenario: Extra response fields are rejected
- **GIVEN** valid JSON containing `smiles` and another top-level field
- **WHEN** the response decoder validates it
- **THEN** it SHALL return `malformed_response/unexpected_fields`
- **AND** it SHALL NOT extract, delete, or rewrite fields or invoke format repair.

#### Scenario: Non-JSON or text-wrapped output is rejected
- **GIVEN** a response with prose, Markdown fences, truncated JSON, text before/after JSON, or otherwise malformed JSON
- **WHEN** strict decoding runs
- **THEN** it SHALL report `invalid_json` without extracting an embedded object or stripping fences
- **AND** the service MAY send one bounded format-repair request containing the original response as data
- **AND** any retry response SHALL pass the same strict JSON/schema decoder.

#### Scenario: A second invalid response fails closed
- **GIVEN** the initial and one repair response both fail strict JSON decoding
- **WHEN** the operation completes
- **THEN** the result SHALL be `malformed_response/invalid_json`
- **AND** no third provider request or ChemIO call SHALL occur.

#### Scenario: Escaped stereochemical backslash is strict JSON
- **GIVEN** a SMILES containing a literal stereochemical backslash
- **WHEN** it is represented in JSON with that backslash escaped as `\\`
- **THEN** strict decoding SHALL recover the exact original SMILES string.

## ADDED Requirements

### Requirement: Structured-output capability fallback SHALL be explicit and bounded

The provider contract SHALL represent prompt-only, OpenAI JSON-object, and unknown structured-output capability. When native JSON output was requested, the provider SHALL retry once without `response_format` only after HTTP 400 carries the exact error code `response_format_not_supported`. Generic HTTP, authentication, authorization, rate-limit, and server errors SHALL NOT trigger that fallback.

#### Scenario: Native JSON is supported
- **GIVEN** native JSON-object mode and a successful HTTP response with valid JSON content
- **WHEN** the provider generates a response
- **THEN** it SHALL send one request containing `response_format={"type":"json_object"}`
- **AND** diagnostics SHALL identify native structured output as used.

#### Scenario: Native JSON is explicitly unsupported
- **GIVEN** native JSON-object mode and HTTP 400 with code `response_format_not_supported`
- **WHEN** the provider handles the response
- **THEN** it SHALL send exactly one second request without `response_format` and with the same strict textual JSON contract
- **AND** a valid response SHALL continue through strict decoding and ChemIO
- **AND** diagnostics SHALL report fallback used and native output unsupported.

#### Scenario: Other HTTP errors do not trigger capability fallback
- **GIVEN** HTTP 400 with another code, HTTP 401/403, rate limiting, or a generic HTTP/server error
- **WHEN** the provider handles the response
- **THEN** it SHALL preserve the existing generic stable error mapping
- **AND** it SHALL NOT issue a fallback request.

#### Scenario: Repair and capability fallback remain bounded together
- **GIVEN** native JSON is unsupported, the single text fallback is malformed, and format repair is eligible
- **WHEN** the operation finishes
- **THEN** there SHALL be at most one native attempt, one text fallback, and one repair invocation (which may itself use at most one capability fallback only if the provider has not already learned prompt-only capability)
- **AND** the operation SHALL never loop or issue a third repair.

### Requirement: Diagnostics SHALL not expose model reasoning or secrets

Structured-output and repair diagnostics SHALL be non-persistent, bounded status metadata only. They MUST NOT include API keys, authorization headers, raw transport exceptions, or `reasoning_content`; the assistant SHALL consume only `message.content`.

#### Scenario: Diagnostic summary is safe
- **GIVEN** a provider response with separate reasoning content, or a failure involving a configured API key
- **WHEN** the UI presents response diagnostics
- **THEN** it SHALL show only capability/repair status booleans or fixed messages
- **AND** it SHALL NOT show or persist reasoning content, credentials, headers, or raw exception text.
