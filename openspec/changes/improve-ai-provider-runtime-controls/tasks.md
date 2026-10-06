## Tasks

- [x] Add provider-neutral `max_tokens` to config and Chat Completions payload with validation/tests.
- [x] Add timeout and token controls to the dialog and pass their values to the provider config.
- [x] Persist endpoint/model/timeout/JSON-mode/max-tokens per provider profile; prove API keys never persist.
- [x] Show elapsed request time without blocking the GUI.
- [x] Map stable provider reason codes to human-facing messages without raw exceptions.
- [x] Run offline provider/UI/settings tests, architecture, OpenSpec, Ruff, and diff checks under hard timeouts.
