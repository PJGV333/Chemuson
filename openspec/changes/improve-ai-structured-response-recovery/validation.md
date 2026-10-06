# Validation

- Strict OpenSpec validation: `openspec validate improve-ai-structured-response-recovery --strict` — valid.
- Offline service/provider/identity tests: 113 passed (`tests/test_molecular_assistant.py`, `tests/test_molecular_assistant_recovery.py`, `tests/test_molecular_identity_verification.py`).
- GUI tests: 20 passed (`tests/test_molecular_assistant_ui.py`); transform tests: 5 passed; architecture suite: 279 passed.
- `python -m compileall -q src tests tools packaging`: passed.
- `pytest --collect-only -q`: 1994 tests collected in 0.64s.
- Scoped Ruff on all changed Python files: passed. Full scoped repository Ruff still reports only the pre-existing unrelated F401 `math` in `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.
- `git diff --check`: passed.
- Loopback availability: `127.0.0.1:1234/v1/models` was unavailable; the already-running `127.0.0.1:8080/v1/models` exposed `qwen3.8-27b`. No server was started and no external endpoint was contacted.
- Bounded live probes against the existing local `qwen3.8-27b`, each configured with a 45-second request timeout: ethanol succeeded (`validation_passed=true`, native mode explicitly rejected and one text fallback used); caffeine succeeded (`validation_passed=true`, prompt-only mode reused); tetrandrine timed out (`reason_code=timeout`). No retry was made for the live timeout. These are runtime observations, not chemistry/identity accuracy claims.
- No monolithic `pytest -q` was run, per the user's explicit bounded-test policy. The separate lifecycle correction discovered in this session is documented in `fix-molecular-assistant-dialog-lifecycle`.
