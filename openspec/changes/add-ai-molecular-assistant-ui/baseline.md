# Phase 2 baseline

Captured on `ai/molecular-assistant-foundation` before production-code changes. `HEAD` and `origin/ai/molecular-assistant-foundation` were both `b9d4d9e67dc4f1e176d418a43fb9b3af63641eb6`; `git status --short` was empty. A recovery checkpoint branch `checkpoint/ai-molecular-assistant-before-phase2` points to this baseline.

Detailed baseline commands and exact outcomes are summarized in this file and `validation.md`; the original verbose transcript has been retired.

| Command | Baseline result |
|---|---|
| `git status --short` | no output (clean) |
| `python -m compileall src tests tools packaging` | exit 0 |
| `pytest --collect-only -q` | `1886 tests collected in 0.75s`; exit 0 |
| `pytest -q` | `1 failed, 1828 passed, 57 skipped in 1166.77s (0:19:26)`; exit 1. Existing failure: `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend`, assertion at `/home/ccachyavgp/Documentos/ChemUSON-UI/tests/test_compchem3d_dock.py:59` (`assert []`). This is a known pre-existing failure and is outside Phase 2. |
| `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` | exit 1 for the existing F401 `math` import in `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; unrelated to this phase. |

The baseline run took approximately 19 minutes. Do not fix the unrelated CompChem or Ruff findings as part of this change.