# AI Provider Runtime Controls Specification

## ADDED Requirements

### Requirement: Bounded provider runtime controls
The request UI SHALL allow per-profile timeout seconds in the inclusive range 10–600 and an output-token cap. The timeout SHALL be passed to the provider configuration. The Chat Completions request SHALL use the provider-neutral `max_tokens` field, defaulting to 4096 and constrained to 64–8192. No vendor-specific fields SHALL be introduced.

#### Scenario: Configure timeout and output bound
- **WHEN** a user selects a profile, sets a timeout and token cap, and submits a request
- **THEN** those values reach `OpenAICompatibleConfig` and the common Chat Completions payload
- **AND** out-of-range values are rejected before provider I/O.

### Requirement: Safe per-profile preferences
The application SHALL persist base URL, model ID, timeout, JSON-output support, and max tokens per profile using existing local settings. It SHALL NOT persist API keys, prompts, or provider responses. Profile endpoint defaults remain editable and llama.cpp's default stays on port 8080.

#### Scenario: Reopen a configured local profile
- **WHEN** a user edits and submits a llama.cpp profile, closes the dialog, and opens it again
- **THEN** the endpoint, model, timeout, JSON-output setting, and token cap are restored
- **AND** the API key field is empty.

### Requirement: Responsive progress and controlled errors
During a request, the UI SHALL remain responsive and show elapsed whole seconds. The UI SHALL translate known stable reason codes to helpful human-readable messages without displaying raw exceptions.

#### Scenario: Long request times out
- **WHEN** a provider returns the stable `timeout` reason
- **THEN** the UI shows the configured timeout duration in a human-readable message
- **AND** it does not display an exception string or freeze the event loop.
