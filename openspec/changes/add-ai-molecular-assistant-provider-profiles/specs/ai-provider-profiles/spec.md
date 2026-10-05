# AI Provider Profiles Specification

## Purpose

Define editable endpoint presets for named OpenAI-compatible providers while preserving the existing generic Chat Completions adapter, strict output validation, and transient credential handling.

## ADDED Requirements

### Requirement: Explicit OpenAI-Compatible Profiles
The system SHALL expose profiles for a custom OpenAI-compatible endpoint, OpenAI, LM Studio, and llama.cpp server. Each named profile SHALL identify its display name, explicit default base URL, API-key requirement, and stable provider ID. Selecting a named profile SHALL NOT change the configured endpoint without user-visible UI state.

#### Scenario: Select a named profile
- **GIVEN** the user selects OpenAI, LM Studio, or llama.cpp
- **WHEN** the request dialog updates
- **THEN** it displays that profile's default base URL in an editable field and uses the profile's stable ID as provider provenance.

#### Scenario: Use a custom endpoint
- **GIVEN** the user selects the custom OpenAI-compatible profile
- **WHEN** the endpoint is entered
- **THEN** the existing generic Chat Completions adapter uses only that explicit endpoint.

### Requirement: Endpoint-Specific Model IDs
The system SHALL require a non-empty, user-supplied model ID and SHALL NOT assume or hard-code a provider model catalog.

#### Scenario: User configures a model
- **GIVEN** a selected endpoint exposes a model under its own identifier
- **WHEN** the user submits that exact identifier
- **THEN** it is sent in the existing Chat Completions `model` field without provider-specific rewriting.

#### Scenario: Model ID is missing
- **GIVEN** the model field is empty
- **WHEN** the user submits the request
- **THEN** no provider worker or network request starts.

### Requirement: Profile-Aware Transient Credentials
The OpenAI profile SHALL require a non-empty API key before worker startup. Local and custom profiles MAY accept an optional key. The key SHALL remain masked and transient, and changing profiles SHALL clear any entered key.

#### Scenario: OpenAI selected without a key
- **GIVEN** the OpenAI profile is selected and its API key is blank
- **WHEN** the user submits a request
- **THEN** the controller refuses to start and the dialog reports a controlled configuration error.

#### Scenario: Profile changes after a key is entered
- **GIVEN** an API key is entered for one profile
- **WHEN** the user selects a different profile
- **THEN** the previous key is cleared before another request can be submitted.

### Requirement: Existing Provider Contract Is Preserved
Every profile SHALL use the existing bounded, non-streaming OpenAI-compatible Chat Completions adapter and strict M23 response decoder. The generic configuration constructor and provider contract SHALL remain backward compatible.

#### Scenario: Profile request is sent
- **GIVEN** any valid profile and model ID
- **WHEN** an offline fake transport observes the request
- **THEN** it sees the profile's explicit base URL plus `/v1/chat/completions`, the exact model ID, and a Bearer header only when a key was provided.

#### Scenario: Response violates the M23 schema
- **GIVEN** a profile endpoint returns non-contract content
- **WHEN** M23 decodes the response
- **THEN** the existing strict decoder returns its stable failure and no graph, without profile-specific extraction or repair.

### Requirement: Offline Compatibility Verification
Tests SHALL verify profile metadata and shared protocol construction with fake transports. Tests SHALL NOT require credentials, internet, or a manually started model server; the profile list SHALL NOT claim live interoperability verification.

#### Scenario: Profile compatibility tests run
- **GIVEN** profile tests run
- **WHEN** fake transport responses are exercised
- **THEN** endpoint, model, authentication-header, provenance, and controlled failure semantics are deterministic and no network request is made.
