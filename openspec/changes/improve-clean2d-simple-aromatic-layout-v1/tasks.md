## 1. Baseline And Investigation

- [x] 1.1 Validate `improve-clean2d-simple-aromatic-layout-v1` with OpenSpec strict before implementation.
- [x] 1.2 Generate `/tmp/clean2d-before.json` before production changes.
- [x] 1.3 Inspect current Clean 2D candidate selection and simple aromatic metrics to choose the smallest safe implementation point.

## 2. Production Improvement

- [x] 2.1 Implement a narrow production change that improves simple aromatic layout selection or acceptance.
- [x] 2.2 Preserve chemical identity, selection metadata, stereo metadata, aromaticity, and existing safety checks.
- [x] 2.3 Avoid behavior changes for high-risk complex-policy and tetrandrine-like protected cases.

## 3. Tests And Baseline Evidence

- [x] 3.1 Add focused tests for measurable simple aromatic improvement and invariant preservation.
- [x] 3.2 Add or reuse tests proving high-risk complex-policy and tetrandrine-like behavior is unchanged.
- [x] 3.3 Generate `/tmp/clean2d-after.json`, run baseline compare, and run baseline review.

## 4. Validation

- [x] 4.1 Run the requested OpenSpec validation and targeted Clean 2D pytest commands.
- [x] 4.2 Run the requested baseline/report and quality/debug/protected pytest commands.
- [x] 4.3 Run a reasonable full test suite if feasible.
