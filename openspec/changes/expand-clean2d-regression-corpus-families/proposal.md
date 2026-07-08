# Proposal: Expand Clean 2D Regression Corpus Families

## Context
`create-clean2d-regression-corpus` established a test-only Clean 2D regression corpus with stable case names, contract assertions, identity preservation checks, diagnostic metric capture, and JSON-serializable snapshots. The initial corpus is intentionally small and needs broader chemical coverage before future changes define strict geometry thresholds or improve layout algorithms.

## Objective
Expand the Clean 2D regression corpus with representative molecule families while keeping the change observational and test-only. The goal is to improve measurement coverage for future deterministic and geometric quality work, not to change Clean 2D behaviour.

## Scope
- Add 10-20 new test-only regression cases to the existing corpus registry.
- Cover substituted aromatics, additional fused aromatics, heteroaromatics, charged molecules, stereo-sensitive molecules, selection-boundary cases, multi-block systems, and simple complex-policy baselines where viable.
- Add or extend case metadata to classify cases as `baseline`, `known_delicate`, `known_failure`, `complex_policy_guard`, `selection_boundary`, and/or `stereo_sensitive`.
- Reuse existing diagnostic metric capture and debug snapshot generation.
- Add tests for minimum family coverage and tag-specific invariants.

## Non-Goals
- No layout algorithm changes.
- No candidate ranking changes.
- No RDKit, CoordGen, or backend changes.
- No UI changes.
- No public document format changes.
- No new strict aesthetic thresholds.
- No opportunistic fixes for benzene, naphthalene, tetrandrine, macrocycles, or partial-selection geometry.

## Success Criteria
- OpenSpec validates in strict mode before and after implementation.
- The expanded corpus runs without uncontrolled Clean 2D exceptions.
- Every case declares stable metadata and tags.
- Charged cases preserve formal charges.
- Stereo-sensitive cases preserve stereo metadata.
- Selection-boundary cases preserve boundary connectivity under the current contract.
- Complex-policy-sensitive cases record current policy behaviour without requiring global redraw.
- All cases can emit JSON-serializable debug snapshots.
