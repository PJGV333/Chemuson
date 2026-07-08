# Design: Clean 2D Baseline Diff Review Policy

## Overview
This change adds a test-only interpretation layer for baseline report diffs. It does not change report generation, baseline records, Clean 2D production behaviour, or any geometry decisions.

## Review Item
Each changed field in a diff becomes a review item containing:
- `case_name`
- `field`
- `category`
- `severity`
- `known_delicate_related`
- `message`

The item is JSON-stable and observational-only.

## Categories And Severities
- `result_state`: category `contract`, severity `contract-risk`.
- `stable_reason`: category `contract`, severity `needs-review`.
- `selected_source`: category `routing`, severity `needs-review`.
- `candidate_sources`: category `routing`, severity `needs-review`.
- `metrics`: category `geometry-diagnostic`, severity `geometry-risk`, observational-only.
- `snapshot`: category `diagnostic`, severity `needs-review`.
- `policy_evidence`: category `complex-policy`, severity `needs-review`.

If the affected case has tag `known_delicate`, the item remains visible and gets `known_delicate_related=True`; its severity may be reported as `known-delicate-change` in summaries while retaining original category context.

## Summary
The summary aggregates:
- total item count;
- counts by severity;
- counts by category;
- known-delicate-related count;
- `observational_only=True`.

## Test Boundary
The implementation belongs under `tests/clean2d_regression/`. It consumes existing report diffs and does not modify reports, baseline records, production layout, ranking, selection, UI, backends, or public document format.

## Future Gates
Future OpenSpec changes may promote selected contract-risk classifications into gates. Geometry-diagnostic differences remain observational in this change.
