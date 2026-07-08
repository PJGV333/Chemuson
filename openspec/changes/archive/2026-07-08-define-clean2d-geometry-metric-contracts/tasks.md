# Tasks: Define Clean 2D Geometry Metric Contracts

## 1. OpenSpec Artifacts
- [x] Create `proposal.md` with metric-contract scope and non-goals.
- [x] Create `design.md` covering registry, metric classes, before/after records, optional values, tolerances, and diagnostic-only status.
- [x] Create `tasks.md` with small implementation steps.
- [x] Create `specs/clean-2d-geometry-metrics/spec.md` with verifiable requirements.
- [x] Run `openspec validate define-clean2d-geometry-metric-contracts --strict` before implementation.

## 2. Test-Only Metric Contract Helpers
- [x] Add a stable metric definition dataclass or equivalent structure.
- [x] Add a registry for all current corpus metrics.
- [x] Add JSON-stable metric record validation.
- [x] Add explicit numeric comparison tolerances.
- [x] Keep helpers diagnostic-only and test-only.

## 3. Corpus Integration Tests
- [x] Verify all corpus-captured metrics have definitions.
- [x] Verify every definition has name, type, polarity, scope, required, and diagnostic-only semantics.
- [x] Verify corpus before/after metric records conform to the contract.
- [x] Verify optional metrics can be `None` only when declared nullable or unavailable.
- [x] Verify JSON serialization uses primitive values only.

## 4. Non-Behavioural Guard Tests
- [x] Verify metric contract tests do not require aesthetic thresholds.
- [x] Verify result-state contract remains independent of geometry metric values.
- [x] Verify this change does not touch ranking, layout, selection, UI, or backends.

## 5. Validation
- [x] Run `openspec validate define-clean2d-geometry-metric-contracts --strict` after implementation.
- [x] Run `pytest tests/test_clean2d_geometry_metric_contracts.py tests/test_clean2d_regression_corpus.py -q`.
- [x] Run requested existing Clean 2D tests when feasible.

## 6. Follow-Up Notes
- [x] Document any geometry-quality issues as baseline, known-delicate, or future threshold work rather than fixing them.
