# Proposal: Persist Clean 2D Baseline Reports

## Why
Clean 2D now has deterministic baseline records for individual regression corpus cases. To compare behaviour across repeated runs, commits, or saved snapshots, the project needs a test-owned report format that aggregates all case records and can produce structured diffs without becoming a geometry gate.

## What Changes
- Add a test-only Clean 2D baseline report builder for the full regression corpus.
- Serialize reports as JSON-stable data with schema and version metadata.
- Sort records by stable case name.
- Include aggregate summaries for case counts, result states, and classification tags.
- Compare two reports using existing deterministic baseline equivalence and metric tolerances.
- Produce structured observational diffs for result state, stable reason, selected source, candidate sources, metrics, snapshot, and policy evidence.

## Non-Goals
- No production changes.
- No layout, ranking, or candidate selection changes.
- No UI changes.
- No RDKit, CoordGen, or backend changes.
- No public document format changes.
- No strict geometry gates.
- No geometry repairs.
- No large `engine.py` refactor.

## Success Criteria
- OpenSpec validates before implementation.
- A full-corpus report serializes with `json.dumps(..., allow_nan=False, sort_keys=True)`.
- Consecutive reports from the same corpus compare equivalent.
- Artificial differences are reported structurally and observationally.
- Existing Clean 2D tests continue to pass.
