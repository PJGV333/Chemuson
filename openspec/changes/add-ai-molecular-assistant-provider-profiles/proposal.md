## Why

M23 and the Phase 2 UI already support an explicitly configured OpenAI-compatible Chat Completions endpoint, but users must manually enter a URL without any provider context. The first UI can safely cover the intended local and hosted services through their shared protocol, provided profiles are clearly identified and model IDs remain endpoint-specific.

## What Changes

- Add a small, immutable catalog of OpenAI-compatible profiles for custom endpoints, OpenAI, LM Studio, and llama.cpp server.
- Let users choose a profile in the existing request dialog; populate its editable default URL and show whether the profile requires an API key.
- Preserve free-text model identifiers because loaded models and their IDs are controlled by each endpoint and change independently.
- Keep all profiles on the existing strict Chat Completions adapter; verify protocol construction offline with fake transports and controller/UI tests.

## Capabilities

### New Capabilities
- `ai-provider-profiles`: Explicit profiles and user-configured model IDs for supported OpenAI-compatible endpoints.

### Modified Capabilities
- None. Provider selection is a new additive capability integrated into the existing transient request dialog; its review and insertion contract remains unchanged.

## Impact

M23 gains profile metadata but no new dependency or provider wire protocol. M10 continues to depend one-way on M23. No API discovery request, model server, persistent setting, API key storage, external dependency, or live compatibility test is introduced. Profiles are compatibility presets, not claims that each server version/model has been live-verified.
