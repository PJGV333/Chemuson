# Design

- Timeout is an integer seconds control, range 10–600, default 60. It is passed to `OpenAICompatibleConfig.timeout_s` and stored independently for each profile.
- Use the shared Chat Completions `max_tokens` field with default 4096 and editable range 64–8192. 4096 leaves headroom for large complete SMILES proposals while providing a finite generation bound; the existing 16 KiB response/SMILES validation limits remain authoritative. No vendor-specific token field is used.
- Persist only profile ID, base URL, model ID, timeout, JSON-output preference, and max tokens through the existing settings store. Never persist request text or API key.
- Show elapsed whole seconds using a GUI timer while the worker runs; do not invent progress percentages.
- Keep internal `reason_code` stable; map known failure codes to user-friendly text and avoid raw exception output.
- llama.cpp remains `127.0.0.1:8080/v1` by default; user changes persist locally per profile.
