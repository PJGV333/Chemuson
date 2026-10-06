# Validation — AI Molecular Assistant integration closure

## Implementation scope

- Removed the three redundant tracked `baseline-output.log` artifacts; retained concise baselines and validation records.
- Removed the obsolete Phase 5 design-block from `AGENT_REPORT.md`, retaining product-history context in `docs/history/CAMPAIGNS.md`.
- Aligned runtime dependencies with `pyproject.toml`, separated development requirements, documented editable installation, and added a consistency test.
- Added selected/candidate/rejected layout summaries to the Clean2D evaluator, without modifying `src/chemuson/clean2d/`; `--api-key-env` is opt-in and never reads ambient keys by default.
- Replaced the generic transformation hook with typed `MolecularTransformationRequest`/`MolecularAssistant.transform()`, added Insert Variant and preserved atomic Replace Original.
- Added configurable provider timeout and output-token controls, non-secret QSettings profile preferences, elapsed status and human-facing failures.
- Added conservative, offline-testable identity verification via the existing Name→Structure resolver and isolated ChemIO InChI worker. Identity is separate from M23 and Clean2D; mismatch needs a second explicit confirmation.
- Cataloged the M10→M16 dependency. No new packages or runtime dependencies; `src/chemuson/clean2d/` is unchanged.

## Offline verification

All pytest commands were bounded to ≤8 minutes; no monolithic pytest run was started (baseline full run was 19:26; operator cap is 10 minutes).

- Current `pytest --collect-only -q`: **1946 collected** in 0.84s (`/tmp/chemuson-session-final-collect.log`).
- `python -m compileall -q src tests tools packaging`: **PASS**.
- `pytest -q tests/architecture`: **278 passed** in 10.34s.
- Focused provider/settings/evaluator/integration block: **20 passed**.
- M23 + identity block: **81 passed**; UI block: **19 passed**; transform block: **5 passed**.
- Ordered Assistant → transform → UI shard around the historic ~69% crash point: **99 passed** in 12.35s.
- Shutdown lifecycle: **2 passed**. Original CompChem→Assistant crash pair: **5/5 forward and 5/5 reverse**, with every invocation separately bounded.
- Broad shards covered all 1946 collected test items in separate bounded invocations: **1922 passed, 20 skipped, 4 failed**. The four failures are the one already-recorded Clean2D candidate baseline and three stereo-import assertions described below.
- Strict OpenSpec: **52 passed, 0 failed**.
- Changed-file Ruff with repository-selected rules (`F401,F811,F821,E722,E741`): **PASS**. Repository-wide same rules reports only the recorded baseline unused `math` import at `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.
- `git diff --check`: PASS.

## Residual test observations

- The untouched `tests/test_clean2d_engine_candidates.py::test_generate_candidates_attempts_rdkit_for_cyclic_graphs` remains failing as already recorded in the GUI-shutdown baseline. No Clean2D source or exception catalogue was altered.
- The installed RDKit 2026.03.6 environment additionally exposes three failures in untouched `tests/test_smiles_stereo_import.py` (two atom-CIP assertions and tetrandrine visual wedge preservation); the vancomycin case passes. These chemistry expectations are outside scope and were not changed. They are not claimed as fixed or conclusively baseline-proven in this continuation.
- One broad `test_[s-z]*.py` invocation reached the 8-minute guard around 87%; it was not repeated. Smaller file-family shards completed. A narrowed ordered S/T/UI/update shard then completed in 7:24 (353 passed, 4 skipped, the same three stereo failures). No SIGSEGV occurred in any bounded current shard. The historical monolithic abort is therefore not reproduced, but is not claimed ruled out.
- A prior combined UI/transform/lifecycle invocation reached its guard. Inspection suggested the new mismatch test's teardown closed an intentionally dirty document and may have opened an unhandled save prompt; the test now marks its undo stack clean before teardown. UI, transform, and lifecycle shards pass independently, and Assistant→transform→UI passes in collection order. The exact timed-out combined command was not repeated, so this is a plausible test-teardown cause, not a confirmed second production SIGSEGV trigger.
- Test-environment quirk: `pytest` is `/usr/bin/pytest` (system Python without RDKit); current RDKit tests used `PYTHONPATH=$PWD/chem/lib/python3.14/site-packages:$PWD/src`. `python` is the project venv and has RDKit but no pytest.

## Local Qwen-only smoke tests

`/v1/models` and `/props` were checked before generation. The only model was `local-model`, path `/mnt/data_sata/models_hf/qwen3.8-27b/Qwen3.8-27B-Ollama-Q4_K_M.gguf`. Requests used only `http://127.0.0.1:8081/v1`, no API key, no remote resolver, a 4096-token cap and ≤5-minute per-request timeout. Summaries under `/tmp/chemuson-qwen-local-smoke/` contain no raw model output or secrets.

| Request | Result | ChemIO | Offline identity |
|---|---|---|---|
| Draw ethanol | success, 3.045 s | accepted | verified against offline fixture |
| Draw caffeine (first attempt) | provider timeout at 120 s | — | — |
| Draw caffeine (retry with 300 s bound) | success, 13.451 s | accepted | verified against offline fixture |
| Draw cholesterol | `malformed_response / invalid_json`, 125.086 s | not reached | not checked |
| Draw vancomycin | `malformed_response / invalid_json`, 87.327 s | not reached | not checked |
| Draw erythromycin | `malformed_response / invalid_json`, 84.380 s | not reached | not checked |

No raw exception, model reasoning, or proposal text was logged. No non-local service, server start, or Clean2D behavior was used.

## Architecture and release state

- No Clean2D, GUI-controller hierarchy, `.cmsn` persistence, chemistry heuristics, or external dependency changes.
- M10→M16 is recorded in `architecture/modules.yml`; architecture tests pass.
- Commits, normal push result, final HEAD and tree status are recorded in the accompanying task report after release operations.
- Status is **implementation complete with explicit verification residuals**; this does not claim the historical monolithic SIGSEGV is eliminated because the complete suite cannot be run within the mandated cap.
