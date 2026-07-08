# Tasks: Persist Clean 2D Baseline Reports

## 1. OpenSpec Artifacts
- [x] Create proposal, design, tasks, and spec delta.
- [x] Validate `persist-clean2d-baseline-reports` before implementation.

## 2. Report Helpers
- [x] Add test-only baseline report dataclass or equivalent structure.
- [x] Build reports for the full regression corpus.
- [x] Sort records by stable case name.
- [x] Add schema/version metadata and caller metadata.
- [x] Add aggregate summary counts.
- [x] Add JSON serialization helper.

## 3. Report Diffing
- [x] Compare reports using deterministic baseline equivalence rules.
- [x] Report added and removed cases.
- [x] Report changed result state, stable reason, selected source, candidate sources, metrics, snapshot, and policy evidence.
- [x] Treat differences as observational in this change.

## 4. Tests
- [x] Verify report JSON stability.
- [x] Verify one record per corpus case.
- [x] Verify stable record ordering.
- [x] Verify summary counts.
- [x] Verify consecutive reports compare equivalent.
- [x] Verify artificial differences are reported.
- [x] Verify metric differences within tolerance are ignored.
- [x] Verify implementation is test-only.

## 5. Validation And Commit
- [x] Run requested OpenSpec validations and pytest suites.
- [x] Commit as `test(clean2d): add persistent baseline reports`.
