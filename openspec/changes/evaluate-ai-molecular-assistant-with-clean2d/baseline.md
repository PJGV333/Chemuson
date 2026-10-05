# Phase 3 baseline

Captured before creating or changing Phase 3 files on branch `ai/molecular-assistant-foundation`. `HEAD` and `origin/ai/molecular-assistant-foundation` were both `aca8ea3557c2bce148a74e39842b60881cefda38`; `git status --short` was empty. A recovery checkpoint branch `checkpoint/ai-molecular-assistant-before-phase3` points to this commit. See [`baseline-output.log`](baseline-output.log) for command outcomes.

| Command | Baseline result |
|---|---|
| `git status --short` | no output (clean) |
| `python -m compileall src tests tools packaging` | exit 0 |
| `pytest --collect-only -q` | `1902 tests collected in 0.82s`; exit 0 |
| `timeout 300s pytest -q tests/test_molecular_assistant.py tests/test_clean2d_safety.py tests/test_clean2d_engine_candidates.py tests/test_clean2d_quality_reporting.py` | `129 passed, 2 skipped in 1.86s`; exit 0 |
| `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` | exit 1 due to the same known F401 `math` import at `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; out of scope |
| Full `pytest -q` | not run. The immediately preceding recorded full-suite baseline takes 19:26 and includes one unrelated CompChem failure; it exceeds the user's 5–10 minute ceiling. No Phase 3 result is inferred from that historical baseline. |

No production files had changed when these commands were captured. All newly introduced Phase 3 tests will be run with a hard timeout no greater than five minutes.
