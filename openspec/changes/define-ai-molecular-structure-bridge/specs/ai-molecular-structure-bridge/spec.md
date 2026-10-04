# AI Molecular Structure Bridge Specification

## Purpose

Define a provider-neutral, chemically validated boundary for proposing molecular structures from a user's natural-language description. This capability is independent of Clean2D and returns ChemUSON's existing `MolGraph` only after the existing ChemIO SMILES importer accepts the model's structured output.

## ADDED Requirements

### Requirement: Provider-neutral application boundary

The molecular assistant SHALL isolate natural-language structure interpretation and provider transport behind an application service and replaceable provider contract.

#### Scenario: Provider can be replaced without chemistry changes
- **GIVEN** a molecular assistant with a provider implementing the declared contract
- **WHEN** the provider implementation is replaced by another compatible provider
- **THEN** response decoding and ChemIO validation SHALL retain the same application contract
- **AND** the chemistry parser SHALL NOT depend on provider-specific SDK types

#### Scenario: Initial compatible endpoint contract
- **GIVEN** the future first HTTP provider adapter
- **WHEN** it sends a structure request
- **THEN** it SHALL use a configured OpenAI-compatible Chat Completions endpoint with non-streaming response content
- **AND** it SHALL expose transport metadata separately from model-generated JSON

### Requirement: Model output SHALL conform to a strict molecular response schema

The assistant SHALL accept exactly one required top-level JSON field, `smiles`, containing a non-empty string within the configured size limit.

#### Scenario: Minimal valid structured proposal
- **GIVEN** provider content containing `{"smiles":"CCO"}`
- **WHEN** the response decoder parses it
- **THEN** the decoder SHALL produce the SMILES value without requiring a generated name or metadata

#### Scenario: Extra response fields are rejected
- **GIVEN** valid JSON containing `smiles` and any other top-level field
- **WHEN** the response decoder validates the object
- **THEN** it SHALL return `malformed_response` rather than silently accepting unknown content

#### Scenario: Non-JSON or text-wrapped output is rejected
- **GIVEN** provider content with prose, Markdown fences, truncated JSON, a non-object JSON value, or an incorrect field type
- **WHEN** the response decoder validates it
- **THEN** it SHALL return a controlled `malformed_response` without extracting or repairing a substring

### Requirement: ChemIO SHALL decide whether a proposed SMILES becomes a molecular graph

The assistant SHALL validate extracted, untrusted SMILES through the existing isolated `chemuson.chemio.rdkit_safe.smiles_to_molgraph_isolated` path and SHALL publish a `MolGraph` only when that path returns a non-empty graph without error.

#### Scenario: Valid SMILES produces the normal ChemUSON graph
- **GIVEN** a structurally valid SMILES proposal
- **WHEN** the ChemIO importer accepts it
- **THEN** the result SHALL contain the graph produced by that importer
- **AND** its chemical atom/bond representation SHALL agree with the ordinary `chemuson.chemio.rdkit_io.smiles_to_molgraph` importer for the same accepted input

#### Scenario: Invalid SMILES is rejected before presentation
- **GIVEN** a schema-valid response whose SMILES is rejected by ChemIO
- **WHEN** validation completes
- **THEN** the result SHALL be `invalid_structure` with no graph
- **AND** the proposed text SHALL NOT be sent directly to a canvas or converted by a second ad-hoc parser

#### Scenario: Parser acceptance is not a claim of semantic correctness
- **GIVEN** ChemIO accepts a graph for a user request such as “draw vancomycin”
- **WHEN** the assistant reports success
- **THEN** it SHALL report parser acceptance only and SHALL NOT claim that the model's structure is scientifically or semantically proven to match the requested molecule

### Requirement: Operation outcomes and failure reasons SHALL be stable

The assistant SHALL expose stable result states `success`, `invalid_request`, `invalid_structure`, `validation_error`, `provider_error`, `malformed_response`, and `cancelled`, with a stable reason code for failures where applicable.

#### Scenario: Invalid request is rejected before transport
- **GIVEN** an empty request or one that exceeds the configured prompt limit
- **WHEN** the application service validates it
- **THEN** the state SHALL be `invalid_request` with a stable reason
- **AND** the provider SHALL NOT be called

#### Scenario: Provider timeout is controlled
- **GIVEN** a provider request exceeds its finite deadline
- **WHEN** the application service returns
- **THEN** its state SHALL be `provider_error` and its stable reason SHALL identify a timeout
- **AND** the raw transport exception SHALL NOT become the public reason code

#### Scenario: Chemistry worker failure is not misreported as invalid chemistry
- **GIVEN** an otherwise schema-valid SMILES and an unavailable or timed-out ChemIO parser worker
- **WHEN** isolated validation cannot finish
- **THEN** the state SHALL be `validation_error` with a stable parser-related reason
- **AND** the result SHALL contain no graph

#### Scenario: Empty and partial responses are controlled
- **GIVEN** the provider returns no content or incomplete structured content
- **WHEN** decoding completes
- **THEN** the state SHALL be `malformed_response` with an appropriate stable reason
- **AND** no partial graph SHALL be returned

#### Scenario: Explicit cancellation has no partial success
- **GIVEN** a request is explicitly cancelled before ChemIO validation succeeds
- **WHEN** the result is produced
- **THEN** its state SHALL be `cancelled` and it SHALL contain no graph

### Requirement: Provider and validation provenance SHALL be observable without unnecessary persistence

Each operation result SHALL expose provider identity, model identity when available, proposed SMILES when extractable, validation outcome, state, and a stable failure reason without requiring persistence of prompt or provider payload.

#### Scenario: Successful validation is traceable
- **GIVEN** a provider returns a valid SMILES and reports its provider/model identity
- **WHEN** ChemIO accepts the structure
- **THEN** the result SHALL expose those identities, the proposed SMILES, and `validation_passed=true`

#### Scenario: Failed validation remains diagnosable
- **GIVEN** a provider response contains an extractable but chemically invalid SMILES
- **WHEN** ChemIO rejects it
- **THEN** the result SHALL expose the provider identity, proposed SMILES, `validation_passed=false`, and stable reason `invalid_smiles`
- **AND** no persistent log or document SHALL be written by default

#### Scenario: Sensitive request data is not logged
- **GIVEN** a request contains user-supplied text and provider credentials
- **WHEN** diagnostics are emitted
- **THEN** prompt text, credentials, authorization headers, and raw HTTP bodies SHALL NOT be logged or persisted by default

### Requirement: Untrusted output SHALL be bounded and non-executable

The assistant SHALL apply strict schema checks, finite provider timeouts and size limits to untrusted request, provider response and SMILES content, and SHALL treat all model content solely as data. Isolated chemistry validation SHALL use the existing ChemIO worker timeout and SHALL fail closed if that worker cannot complete.

#### Scenario: Oversized response is rejected before parsing
- **GIVEN** a provider response exceeds a configured byte/character limit
- **WHEN** the service receives it
- **THEN** it SHALL return a controlled failure without truncating the response or invoking ChemIO

#### Scenario: Generated content is never executed
- **GIVEN** model content includes Python, shell, tool-call-like text, or instructions outside the response schema
- **WHEN** the service processes it
- **THEN** it SHALL reject non-contract content and SHALL NOT execute it or call arbitrary application methods

#### Scenario: Endpoint is explicitly configured
- **GIVEN** a user prompt or model response contains a URL or tool instruction
- **WHEN** a provider request is prepared
- **THEN** the service SHALL use only the explicitly configured provider endpoint and SHALL NOT derive destinations or tools from generated content

### Requirement: Failed generation SHALL not mutate a document or canvas

The molecular assistant SHALL have no document or GUI mutation capability and SHALL return no usable graph on any non-success state.

#### Scenario: Provider or parser failure preserves the active molecule
- **GIVEN** an existing molecule and a future caller that commits only successful results
- **WHEN** the provider fails, the response is malformed, the request is cancelled, or ChemIO rejects the SMILES
- **THEN** the caller SHALL not invoke graph insertion and the existing molecule, selection, undo state and dirty state SHALL remain unchanged

#### Scenario: Successful result is committed only after validation
- **GIVEN** a provider-produced response
- **WHEN** ChemIO has not yet accepted it
- **THEN** no graph insertion or document mutation SHALL occur

### Requirement: Clean2D SHALL remain independent of AI providers

M02 Clean2D SHALL remain usable without a model and SHALL have no dependency on the molecular assistant or provider infrastructure.

#### Scenario: Clean2D runs without AI infrastructure
- **GIVEN** a ChemUSON installation with no provider configured or reachable
- **WHEN** a user invokes Clean2D on an existing `MolGraph`
- **THEN** Clean2D SHALL operate through its existing deterministic pipeline without importing or contacting AI infrastructure

#### Scenario: Dependency direction remains downstream
- **GIVEN** the future module catalog and import-boundary tests
- **WHEN** they inspect AI, ChemIO and Clean2D dependencies
- **THEN** the assistant MAY depend on ChemIO/Core, GUI may later depend on the assistant, and M02 SHALL NOT depend on the assistant or provider adapters

### Requirement: Provider tests SHALL be deterministic and offline

The assistant's contract tests SHALL use fake providers/transports and SHALL NOT require a real API, local model server, network access, or credentials.

#### Scenario: Fake provider exercises a valid structure
- **GIVEN** a fake provider returning a valid structured SMILES proposal
- **WHEN** the assistant service is tested
- **THEN** the test SHALL exercise the response contract and normal ChemIO graph construction without Internet access

#### Scenario: Fake provider exercises controlled failures
- **GIVEN** fakes for invalid SMILES, malformed/empty/extra-text response, timeout, and network failure
- **WHEN** the assistant tests run
- **THEN** each case SHALL assert the stable state/reason and absence of graph on failure

#### Scenario: Import and dependency boundaries are checked statically
- **GIVEN** the future assistant module and architecture catalog
- **WHEN** architecture tests inspect imports
- **THEN** the assistant SHALL NOT import GUI, canvas, Clean2D, ChemName, `name2structure` or execute `tools` code, and Clean2D SHALL NOT import AI
