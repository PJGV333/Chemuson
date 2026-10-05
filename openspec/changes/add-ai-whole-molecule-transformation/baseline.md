# Phase 5 baseline

Captured on local branch `ai/molecular-assistant-foundation` at `e0a802cb8a6a439294f8fb209e88ca643166272e` (`fix(gui): wait for owned workers before close`) before creating this Phase 5 change. The recovery checkpoint `checkpoint/ai-molecular-assistant-before-phase5` points to the same commit. `git status --short` was empty. `origin/ai/molecular-assistant-foundation` still points to `2de88ce03d786bc386bd452ea828906f5f81e325`; the prior lifecycle push could not authenticate. No source changes were made between the recorded validation and this baseline.

The command output/logs remain under `/tmp` and are not versioned.

| Required baseline command | Result |
|---|---|
| `git status --short` | no output; clean at `e0a802cb8a6a439294f8fb209e88ca643166272e` |
| `python -m compileall src tests tools packaging` | exit 0; `/tmp/chemuson-worker-shutdown-final-compileall.log` |
| `pytest --collect-only -q` | `1918 tests collected in 0.61s`; `/tmp/chemuson-worker-shutdown-collect.log` |
| `pytest -q` | exit 139 (SIGSEGV) near 69%, at the same Molecular Assistant `QSignalSpy.wait()` test and Qt `QUndoStack`/widget destruction stack as the pre-lifecycle baseline; two prior known failures are the Clean2D candidate-source assertion and CompChem fake-backend timeout. Logs: `/tmp/chemuson-worker-shutdown-full-suite.log` and `/tmp/chemuson-worker-shutdown-full-verbose.log`; pre-change comparator: `/tmp/chemuson-phase45-baseline-pytest.log` |
| `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` | exit 1 for the existing F401 `math` import in `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; `/tmp/chemuson-worker-shutdown-final-ruff.log` |

Relevant baseline regression evidence: the original CompChem-export then Molecular-Assistant worker pair now passes five times in both orders after the lifecycle commit; the 111-test offline worker/provider/geometry/UI block passes. These are pre-Phase-5 checks, not claims of full-suite success. All provider and network paths remain fake/offline.
