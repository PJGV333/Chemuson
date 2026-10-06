# Phase 3 validation record

## Baseline

- Branch `ai/molecular-assistant-foundation`, clean tree at `aca8ea3557c2bce148a74e39842b60881cefda38`; same as `origin`.
- Baseline: compileall passed; 1,902 tests collected; focused M23/Clean2D tests `129 passed, 2 skipped` in 1.86s.
- Repository-wide Ruff baseline has the existing unrelated F401 `math` import in `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.
- Full pytest was not run: the prior recorded baseline takes 19:26 and contains a pre-existing CompChem failure, exceeding the requested 5–10 minute maximum. See `baseline.md`; its summary contains the recorded duration and baseline failure.

## Focused verification

- `timeout 300s pytest -q tests/test_ai_clean2d_evaluation.py tests/test_molecular_assistant.py tests/test_clean2d_safety.py tests/test_clean2d_engine_candidates.py tests/test_clean2d_quality_reporting.py`: **136 passed, 2 skipped** in 2.10s.
- `timeout 300s pytest -q tests/architecture/test_molecular_assistant_boundary.py tests/architecture/test_import_boundaries.py tests/architecture/test_module_catalog.py`: **91 passed** in 6.97s.
- The isolated evaluator+boundary tests also passed: **14 passed** in 0.82s, including an actual small-graph Clean2D engine invocation with a 0.25s RDKit backend timeout.
- `python -m compileall src tests tools packaging`: **passed**.
- `pytest --collect-only -q`: **1,910 tests collected** in 0.76s.
- Scoped Ruff on the evaluator and changed test files: **passed**.
- `openspec validate evaluate-ai-molecular-assistant-with-clean2d --strict`: **valid**.
- `openspec validate --all --strict`: **46 passed, 0 failed**.
- `git diff --check`: **passed**.

## Architectural outcome and test limits

The evaluator is only in `tools/`; it composes existing M23 and M02 APIs, while M02 remains free of M23 imports and neither production module changes. No module catalog, production behavior, external dependency, network service, or GUI/document state was changed. All evaluator provider tests are offline fakes. The full suite was deliberately not rerun due to its recorded 19:26 duration and known baseline failure.
