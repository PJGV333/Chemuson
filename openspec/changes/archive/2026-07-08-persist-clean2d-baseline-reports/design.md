# Design: Clean 2D Baseline Reports

## Overview
This change adds a test-only report layer on top of deterministic baseline records. Reports are intended for reproducible comparison and developer diagnostics, not for production use or geometric acceptance gates.

## Report Shape
A report contains:
- `schema`: stable schema identifier.
- `version`: integer schema version.
- `metadata`: run metadata supplied by tests or callers.
- `summary`: aggregate counts.
- `records`: one canonical baseline record per corpus case.

The report is JSON-stable and sorted by case name.

## Summary
The summary includes:
- total case count;
- result-state counts;
- tag counts for known classifications such as `known_delicate`, `complex_policy_guard`, `selection_boundary`, `stereo_sensitive`, and `charged`;
- family counts.

## Diffing
Report comparison produces a structured diff with:
- `equivalent`: boolean;
- `observational_only`: true for this change;
- `added_cases` and `removed_cases`;
- `changed_cases`, each with changed field names and before/after values.

The diff checks result state, stable reason, selected source, candidate sources, metrics, snapshot, and policy evidence. Metrics use existing tolerance-aware comparison from deterministic baselines.

## JSON Stability
Reports use canonicalized records from the baseline harness. No timestamps, temporary paths, or absolute paths are required. Caller-provided metadata is canonicalized and must be JSON-stable.

## Test-Only Boundary
The implementation belongs under `tests/clean2d_regression/` and does not modify production code, layout, ranking, selection, UI, backends, or public document format.

## Follow-Up Roadmap
- Add an optional developer script to write reports to disk if needed.
- Add CI artifact generation only after report shape stabilizes.
- Promote selected comparisons into gates only in a future OpenSpec change.
