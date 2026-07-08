# Tasks: Define Clean 2D Baseline Diff Review Policy

## 1. OpenSpec Artifacts
- [x] Create proposal, design, tasks, and spec delta.
- [x] Validate `define-clean2d-baseline-diff-review-policy` before implementation.

## 2. Review Policy Helpers
- [x] Add test-only review item structure.
- [x] Classify baseline diff fields into stable categories and severities.
- [x] Label known-delicate changes without hiding them.
- [x] Produce JSON-stable review summaries.

## 3. Tests
- [x] Verify empty diff produces no review items.
- [x] Verify result-state changes classify as contract-risk.
- [x] Verify reason changes classify as needs-review.
- [x] Verify selected source and candidate source changes classify as routing review.
- [x] Verify metric changes classify as geometry-diagnostic and observational-only.
- [x] Verify snapshot changes classify as diagnostic review.
- [x] Verify policy evidence changes classify as complex-policy review.
- [x] Verify known-delicate changes remain visible and labeled.
- [x] Verify policy helpers do not modify reports or production code.

## 4. Validation And Commit
- [x] Run requested OpenSpec validations and pytest suites.
- [x] Commit as `test(clean2d): add baseline diff review policy`.
