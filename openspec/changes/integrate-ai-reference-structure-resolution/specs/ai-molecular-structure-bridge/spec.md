## MODIFIED Requirements

### Requirement: Provider and validation provenance SHALL be observable without unnecessary persistence
Each operation result SHALL expose provider identity, model identity when available, proposed SMILES when extractable, validation outcome, state, and a stable failure reason without requiring persistence of prompt or provider payload. It MAY additionally expose an allowlisted finish reason and bounded numeric completion/reasoning token counts; it MUST NOT retain or expose reasoning content.

#### Scenario: Successful validation is traceable
- **GIVEN** a provider returns a valid SMILES and reports its provider/model identity
- **WHEN** ChemIO accepts the structure
- **THEN** the result SHALL expose those identities, the proposed SMILES, and `validation_passed=true`.

#### Scenario: Failed validation remains diagnosable
- **GIVEN** a provider response contains an extractable but chemically invalid SMILES
- **WHEN** ChemIO rejects it
- **THEN** the result SHALL expose the provider identity, proposed SMILES, `validation_passed=false`, and stable reason `invalid_smiles`
- **AND** no persistent log or document SHALL be written by default.

#### Scenario: Sensitive request data is not logged
- **GIVEN** a request contains user-supplied text and provider credentials
- **WHEN** diagnostics are emitted
- **THEN** prompt text, credentials, authorization headers, and raw HTTP bodies SHALL NOT be logged or persisted by default.

#### Scenario: Empty content exhausts the generation budget
- **GIVEN** an OpenAI-compatible response with empty `message.content`, `finish_reason="length"`, and bounded usage counts
- **WHEN** the Molecular Assistant processes it
- **THEN** the result reports the stable reason `generation_exhausted` and available bounded counts
- **AND** no format-repair request is sent because there is no content to repair.

#### Scenario: Reasoning is present separately from content
- **GIVEN** a response envelope includes `reasoning_content`
- **WHEN** provider metadata is extracted
- **THEN** `reasoning_content` is neither stored in `ProviderResponse`/result nor logged, displayed, or persisted.

#### Scenario: Provider supplies untrusted diagnostic values
- **GIVEN** finish-reason or usage values are malformed, negative, non-integral, or outside configured bounds
- **WHEN** the response is normalized
- **THEN** unsafe values are discarded or mapped to a fixed allowlisted value without failing open.
