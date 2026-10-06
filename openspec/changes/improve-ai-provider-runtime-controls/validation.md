# Validation — provider runtime controls

## Implementation

- Shared Chat Completions `max_tokens` defaults to 4096, range 64–8192; request timeout is 10–600 seconds.
- Dialog controls are non-blocking and pass values into the request configuration. Per-profile base URL, model, timeout, JSON mode, and token cap are stored in QSettings; API keys remain transient and are removed from legacy persistent settings.
- Elapsed request time and user-facing timeout/network/structured-response messages are covered without exposing raw provider exceptions.

## Checks

- Provider, preferences, evaluator, integration contracts: **20 passed**.
- Molecular Assistant service + identity tests: **81 passed**; assistant UI: **19 passed**; typed transform: **5 passed**.
- Architecture suite: **278 passed**.
- `python -m compileall -q src tests tools packaging`: **PASS**.
- Strict OpenSpec: **52 passed, 0 failed**.
- Changed-file Ruff using repository-selected checks: **PASS**. Repository-wide selected Ruff retains only baseline `F401 math` at `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.
- Offline tests contain no live provider/network calls. Local-only Qwen smoke results are recorded in `../stabilize-ai-molecular-assistant-integration/validation.md`.

No runtime dependency, Clean2D, persistence-format, or API-key storage changes were introduced.
