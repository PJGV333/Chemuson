## Context

`OpenAICompatibleProvider` already sends bounded, non-streaming `POST /v1/chat/completions` requests with strict response decoding downstream. The configuration accepts an explicit base URL, model ID, transient optional key, and provider ID. The Phase 2 dialog currently asks users to enter all values manually. The initial product scope named OpenAI, LM Studio, and llama.cpp-compatible endpoints; these share the existing wire protocol.

## Goals / Non-Goals

**Goals:**
- Offer editable endpoint presets for the three named services plus a custom OpenAI-compatible endpoint.
- Carry a stable profile ID into M23 provenance while leaving the generic adapter and strict decoder in control.
- Keep model IDs endpoint-selected, editable, and not hard-coded to a changing provider catalog.
- Make API-key requirements visible and reject a missing OpenAI key before starting a worker.
- Keep model requests explicit, asynchronous, transient, and testable entirely offline.

**Non-Goals:**
- Implement native Anthropic, Ollama, or other non-Chat-Completions protocols.
- Call `/models`, enumerate models, claim that every model is compatible, or add model-ID suggestions that may become stale.
- Store endpoint preferences, API keys, prompts, or provider history.
- Change response schema, ChemIO validation, Clean2D, review, insertion, or undo contracts.
- Connect to a live service or make manual interoperability assertions.

## Decisions

1. **Use explicit profiles over multiple adapters.** Profiles are metadata plus defaults consumed by the existing OpenAI-compatible adapter; no per-vendor transport behavior is introduced. `custom` remains the default, preserving the current blank endpoint behavior. Profiles: `openai` → `https://api.openai.com/v1` (key required), `lm-studio` → `http://127.0.0.1:1234/v1` (local key optional), `llama-cpp` → `http://127.0.0.1:8080/v1` (local key optional), and `openai-compatible` → caller-supplied URL. Every endpoint remains editable.

2. **Do not hard-code model names or add model discovery.** Users enter the exact model ID exposed by their selected endpoint. Local servers can load arbitrary models, and hosted catalogs evolve; a universal list would be stale or misleading. A future discovery feature requires its own API/compatibility contract and asynchronous UI design.

3. **Make credentials profile-aware and transient.** OpenAI profile selection marks the key as required and the controller refuses to start if it is blank. Local/custom profiles retain optional keys. Changing profiles clears any entered key to prevent accidentally forwarding credentials to a different host. Keys remain masked and are erased after submission as before.

4. **Carry profile identity through existing config.** The selected profile ID populates `OpenAICompatibleConfig.provider_id`; `MolecularAssistantResult.provider_id` and the review's existing provenance label then identify the configured profile. `OpenAICompatibleConfig` remains backward compatible with `provider_id='openai-compatible'` by default.

5. **Measure protocol shape offline, not server availability.** Fake transports verify each configured profile produces the expected `/v1/chat/completions` URL, model field, optional bearer header, and strict response handling. This establishes conformance to the shared adapter contract without network access; actual server-version interoperability is explicitly unverified.

## Risks / Trade-offs

- [A provider changes its port or model ID format] → URL remains editable and model ID is always user-supplied; do not promise live discovery.
- [An OpenAI key is accidentally sent to another endpoint] → profile changes clear the key; endpoint is explicit/editable, and no key is persisted.
- [A server calls itself OpenAI-compatible but differs in envelope details] → strict decoder fails closed with existing stable errors; profiles do not add heuristic parsing.
- [Profile presets are mistaken for live-verified integrations] → UI and spec state that these are protocol presets; live server testing is not performed in this phase.

## Migration Plan

No migration. The existing UI's endpoint/model workflow is preserved through the `custom` profile, and old direct controller calls keep the default profile ID. Removing the selector or profile catalog leaves the generic OpenAI-compatible provider usable.

## Open Questions

None for OpenAI-compatible profiles. Native provider protocols and model discovery remain separate scope.
