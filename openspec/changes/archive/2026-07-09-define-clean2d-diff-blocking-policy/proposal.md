## Why

Clean 2D baseline diff review can already surface observable changes, but reviewers need a stable way to interpret their severity without turning those observations into automatic gates. This change defines a test-only, documented policy for mapping review summary items to human-review decisions.

## What Changes

- Add a new Clean 2D diff blocking policy capability that classifies existing diff review summary items into structured decisions.
- Define severity categories for `blocking-candidate`, `review-required`, `informational`, `allowed-known-delicate`, and `future-gate-candidate` use.
- Add test-only helper coverage that consumes existing review summaries/items and produces JSON-stable decisions for reviewer guidance.
- Keep `blocking-candidate` advisory only; no CI gate, production behavior, layout, ranking, selection, UI, backend, public format, or geometry changes are introduced.

## Capabilities

### New Capabilities
- `clean-2d-diff-blocking-policy`: Defines the stable, test-only policy for classifying Clean 2D baseline diff review summaries into structured human-review decisions.

### Modified Capabilities

## Impact

- Affected test-only code: `tests/clean2d_regression/blocking_policy.py`, `tests/test_clean2d_diff_blocking_policy.py`.
- Affected planning docs: OpenSpec change artifacts under `openspec/changes/define-clean2d-diff-blocking-policy/`.
- No production modules, public document formats, Clean 2D layout/ranking/selection behavior, UI code, RDKit/CoordGen integrations, or CI gates are changed.
