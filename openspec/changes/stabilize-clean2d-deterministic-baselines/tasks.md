# Tasks: Stabilize Clean 2D Deterministic Baselines

## 1. OpenSpec Artifacts
- [x] Create proposal, design, tasks, and spec delta.
- [x] Validate `stabilize-clean2d-deterministic-baselines` with `openspec validate --strict` before implementation.

## 2. Baseline Harness
- [x] Add test-only baseline record creation helpers.
- [x] Canonicalize dictionaries, sets, sequences, floats, metrics, and snapshots.
- [x] Preserve observable result-state fields without recomputing them from metrics.
- [x] Keep metrics diagnostic-only.

## 3. Determinism Tests
- [x] Assert each corpus case emits a JSON-stable baseline record.
- [x] Assert baseline records contain required metadata.
- [x] Assert repeated runs of the same case compare equivalent after canonicalization.
- [x] Assert float comparison uses explicit tolerances.
- [x] Assert corpus case order is stable.
- [x] Assert candidate source labels are stable or documented as known observational instability.

## 4. Non-Behavioural Guards
- [x] Verify the change does not modify production layout, ranking, selection, UI, backends, or public document format.
- [x] Document any observed instability as follow-up rather than fixing geometry or ranking.

## 5. Validation
- [x] Run `openspec validate stabilize-clean2d-deterministic-baselines --strict`.
- [x] Run `pytest tests/test_clean2d_deterministic_baselines.py tests/test_clean2d_regression_corpus.py tests/test_clean2d_geometry_metric_contracts.py -q`.
- [x] Run requested existing Clean 2D tests.
- [x] Run `openspec validate --all --strict` if feasible.

## 6. Commit
- [x] Commit as `test(clean2d): add deterministic baseline harness`.
