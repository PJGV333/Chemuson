# Phase 2 validation record

## Baseline

- Branch `ai/molecular-assistant-foundation`; `HEAD` and `origin/ai/molecular-assistant-foundation`: `b9d4d9e67dc4f1e176d418a43fb9b3af63641eb6`.
- Tree clean before changes. Full baseline transcript: `baseline-output.log`.
- Baseline full suite: 1 unrelated CompChem failure, 1,828 passed, 57 skipped; runtime 19:26. Baseline repository Ruff: one unrelated F401 in the Clean2D test.

## Focused verification

- `pytest -q tests/architecture`: **277 passed** in 10.33s.
- `pytest -q tests/test_command_palette.py`: **27 passed** in 18.58s.
- `pytest -q tests/test_main_window_tabs.py::test_import_smiles_inserts_without_clearing_and_undo`: **1 passed** in 0.87s.
- Offline Phase 2 UI/controller/canvas test cases were run individually to respect the requested test-time cap: **14 passed** total, covering menu/palette discovery, modeless/masked configuration, worker-thread execution, invalid config, abandoned late result suppression, six failure states preserving editor snapshots, preview-before-insert dispatch, decline/late result preservation, and one-step canvas undo/redo.
- `python -m compileall src tests tools packaging`: **passed**.
- `pytest --collect-only -q`: **1902 tests collected**.
- Scoped Ruff on changed production/tests with `F401,F811,F821,E722,E741`: **passed**.
- `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`: retains the same baseline F401 at `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3` (`math`), unchanged and out of scope.
- `pytest -q tests/architecture/test_molecular_assistant_ui_boundary.py tests/architecture/test_import_boundaries.py tests/architecture/test_module_catalog.py tests/architecture/test_molecular_assistant_boundary.py`: **92 passed** in 6.90s.
- `openspec validate add-ai-molecular-assistant-ui --strict`: **valid**.
- `openspec validate --all --strict`: **45 passed, 0 failed**.
- `git diff --cached --check`: **passed** on the complete staged change.

## Time-bounded test decision

The full suite was not rerun: the recorded baseline itself takes 19:26 and contains the known unrelated failure, exceeding the operator's explicit 5–10 minute maximum. An initial combined window-level insertion/undo test was stopped; verification was split into a short fake-canvas route test plus a direct normal-canvas undo/redo test, both completing in under one second. No manual provider test, network request, model server, or UI service was run.