## 1. Baseline And Scope

- [x] 1.1 Validate the OpenSpec change with `openspec validate improve-clean2d-monosubstituted-aromatic-layout-v1 --strict` before implementation.
- [x] 1.2 Generate `/tmp/clean2d-before-monosub-aromatic.json` before production changes.
- [x] 1.3 Inspect current metrics for toluene-like, phenol-like, and aniline-like corpus cases.

## 2. Production Candidate

- [x] 2.1 Add a narrowly scoped monosubstituted aromatic template candidate for quick/publication modes only.
- [x] 2.2 Detect only one neutral 5/6-member aromatic ring with one neutral terminal substituent and no stereo, fused, interaction, or selection-boundary risk.
- [x] 2.3 Preserve chemical identity and rely on existing safety/ranking checks.

## 3. Tests And Baseline Evidence

- [x] 3.1 Add focused tests for measurable monosubstituted aromatic improvement and non-current selected source.
- [x] 3.2 Add tests for outward substituent orientation and invariant preservation.
- [x] 3.3 Add tests that the new candidate does not appear in propose, simple aromatics regressions, fused/high-risk, or tetrandrine-like coverage.
- [x] 3.4 Generate after baseline report, compare, and review.

## 4. Validation

- [x] 4.1 Run requested OpenSpec and directed pytest validation.
- [x] 4.2 Run requested baseline/report and quality/debug/protected pytest validation.
- [x] 4.3 Run full `pytest -q` and confirm only the documented preexisting failures remain.
