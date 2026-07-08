# Proposal: Define Clean 2D Baseline Diff Review Policy

## Why
Clean 2D baseline reports can now be generated and compared, but report differences need a stable test-owned interpretation policy. Before introducing any gates, reviewers need structured severity and category labels for observable differences.

## What Changes
- Add a test-only policy for classifying diffs produced by `compare_baseline_reports`.
- Classify contract, routing, geometry-diagnostic, snapshot, and complex-policy changes.
- Keep known-delicate case changes visible while labeling them separately.
- Produce JSON-stable review items and summaries.
- Document which categories may be candidates for future gates without activating gates now.

## Non-Goals
- No production changes.
- No layout, ranking, selection, UI, backend, RDKit, or CoordGen changes.
- No public document format changes.
- No geometry gates or CI failures for visual quality.
- No geometry repairs.
- No large `engine.py` refactor.

## Success Criteria
- OpenSpec validates before implementation.
- Empty diffs produce equivalent review summaries with no items.
- Artificial field changes are classified into stable categories and severities.
- Known-delicate changes remain visible and are tagged as known-delicate related.
- Review summaries serialize with JSON-stable primitives.
