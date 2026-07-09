## 1. Policy Helper

- [x] 1.1 Add test-only Clean 2D diff blocking policy helper under `tests/clean2d_regression/blocking_policy.py`.
- [x] 1.2 Return a structured, JSON-stable decision without mutating the input review summary.
- [x] 1.3 Keep the helper advisory-only and free of production imports or gate enforcement.

## 2. Test Coverage

- [x] 2.1 Add tests for empty summaries, contract-risk, needs-review, geometry-risk, known-delicate-only, and mixed known-delicate plus contract-risk cases.
- [x] 2.2 Add tests for deterministic JSON serialization and non-mutation of input summaries.
- [x] 2.3 Add tests that confirm no production behavior or automatic gate is introduced.

## 3. Validation

- [x] 3.1 Run `openspec validate define-clean2d-diff-blocking-policy --strict` before implementation.
- [x] 3.2 Run targeted policy/report pytest coverage after implementation.
- [x] 3.3 Run full requested Clean 2D regression pytest coverage.
- [x] 3.4 Run `openspec validate --all --strict`.
