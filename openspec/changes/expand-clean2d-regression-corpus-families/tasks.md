# Tasks: Expand Clean 2D Regression Corpus Families

## 1. OpenSpec Artifacts
- [x] Create `proposal.md` with test-only scope and non-goals.
- [x] Create `design.md` documenting metadata, family coverage, assertions, metrics, gaps, and follow-up roadmap.
- [x] Create `tasks.md` with small implementation phases.
- [x] Create `specs/clean-2d-regression-corpus/spec.md` with verifiable requirements.
- [x] Run `openspec validate expand-clean2d-regression-corpus-families --strict` before touching tests.

## 2. Case Metadata
- [x] Add stable case tags while preserving existing case metadata.
- [x] Mark existing cases with appropriate baseline, known-delicate, selection-boundary, stereo-sensitive, or complex-policy tags.
- [x] Add tests that require each case to declare tags.

## 3. New Family Builders
- [x] Add substituted aromatic builders.
- [x] Add additional fused aromatic builder(s).
- [x] Add heteroaromatic builders.
- [x] Add charged molecule builders.
- [x] Add stereo-sensitive builders.
- [x] Add selection-boundary builders.
- [x] Add multi-block builders.
- [x] Add simple macrocycle/cavity/cyclophane baseline if viable, otherwise document the gap.

## 4. Registry Expansion
- [x] Add 10-20 new cases to the registry.
- [x] Prioritize diversity over quantity.
- [x] Keep expected states observational and aligned with current controlled behaviour.
- [x] Mark delicate or failure cases without adding geometry fixes.

## 5. Assertion Expansion
- [x] Add minimum coverage tests by family/tag.
- [x] Assert charged cases preserve formal charges.
- [x] Assert stereo-sensitive cases preserve stereo metadata.
- [x] Assert selection-boundary cases preserve boundary bonds.
- [x] Assert complex-policy-sensitive cases do not require global redraw when current policy blocks it.
- [x] Keep metric capture diagnostic-only.

## 6. Validation
- [x] Run `openspec validate expand-clean2d-regression-corpus-families --strict` after implementation.
- [x] Run `pytest tests/test_clean2d_regression_corpus.py -q`.
- [x] Run relevant existing Clean 2D tests when feasible.

## 7. Follow-Up Notes
- [x] Document any geometry-quality failures as baseline, known-delicate, known-failure, or follow-up work rather than fixing them.
