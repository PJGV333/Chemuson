# Phase 4 baseline

Captured before Phase 4 file changes on `ai/molecular-assistant-foundation`. `HEAD` and `origin/ai/molecular-assistant-foundation` were both `5f8f5386ab6dea4b084c7e00597ee34d533728c6`; `git status --short` was empty. Recovery checkpoint `checkpoint/ai-molecular-assistant-before-phase4` points to this commit. Command results and known findings are summarized in `validation.md`; verbose command output was not retained in the repository.

| Command | Baseline result |
|---|---|
| `git status --short` | no output (clean) |
| `python -m compileall src tests tools packaging` | exit 0 |
| `pytest --collect-only -q` | `1910 tests collected in 0.82s`; exit 0 |
| `timeout 300s pytest -q tests/test_molecular_assistant.py tests/test_molecular_assistant_ui.py tests/architecture/test_molecular_assistant_ui_boundary.py tests/architecture/test_molecular_assistant_boundary.py` | `85 passed, 2 skipped in 4.85s`; exit 0 |
| `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` | exit 1 on the known F401 `math` import in `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; unchanged/out of scope |
| Full `pytest -q` | not run; preceding captured baseline is 19:26 with an unrelated pre-existing CompChem failure and exceeds the requested 5–10 minute maximum. |

No tracked or untracked Phase 4 files existed during this baseline.
