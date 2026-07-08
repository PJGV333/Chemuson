# Design: Clean 2D Deterministic Baselines

## Overview
This change introduces a test-only harness that executes regression corpus cases and turns each run into a canonical baseline record. The baseline answers whether two equivalent executions produce the same observable Clean 2D report without changing production behaviour.

## Baseline Record
The baseline record contains:
- `case_name`
- `family`
- `tags`
- `mode`
- `target`
- `result_state`
- `stable_reason`
- `selected_source`
- `candidate_sources`
- `metrics`
- `snapshot`
- `policy_evidence`, when present in candidate metadata or snapshot diagnostics

The record is test-owned and JSON-serializable. It is not a public document format.

## Canonicalization
The harness canonicalizes:
- dictionaries by sorted string keys;
- sets and frozensets as sorted lists;
- tuples and lists recursively;
- floats by finite-value validation and metric-aware tolerance comparison;
- `None` as a stable value;
- candidate sources as the observed sequence converted to immutable tuples;
- snapshot metadata with test-owned case metadata preserved and path/timestamp-like fields excluded if present.

`NaN` and infinity are invalid in baseline records.

## Comparison
Two runs of the same case are compared field-by-field after canonicalization. Numeric metric fields use explicit tolerances from `tests.clean2d_regression.metrics`. Non-metric floats in snapshots are rounded to a stable precision for canonical comparison.

Changes in `result_state`, `stable_reason`, `selected_source`, or candidate source ordering are not hidden. If such instability is observed, the case should be marked as known observational instability or documented as follow-up.

## Snapshot Handling
Snapshots are included after canonicalization. If future snapshot fields contain timestamps, temporary paths, absolute paths, or ephemeral identifiers, the baseline helper may exclude those non-semantic fields and must document the excluded key names in tests.

## Test Scope
The tests verify JSON stability, minimum metadata, repeated execution equivalence, float tolerance use, stable corpus ordering, candidate source stability, and diagnostic-only metric behaviour. The tests must not fail solely because current geometry is visually imperfect.

## Follow-Up Roadmap
- Add persisted baseline report files only after the record shape proves stable.
- Investigate any known observational instability in separate OpenSpec changes.
- Define stricter geometry gates only after deterministic baselines are established.
