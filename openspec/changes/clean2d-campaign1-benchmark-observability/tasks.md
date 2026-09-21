# Tasks: Clean2D Campaign 1 Benchmark and Observability

## OpenSpec and baseline

- [x] Read the master Clean2D campaign policy and current Clean2D contracts.
- [x] Capture the pre-Campaign-1 baseline in this change.
- [x] Validate this OpenSpec strictly before implementation.

## Benchmark foundation

- [x] Add explicit size class to every stable regression case without renaming existing IDs.
- [x] Derive stable topology metadata from each built graph.
- [x] Extend baseline records with metric-vector and runtime evidence while preserving current fields.
- [x] Preserve deterministic canonicalization and exclude runtime-only fields from equivalence.
- [x] Keep candidate sources and result states auditable.

## Tests and evidence

- [x] Add failing contract tests for taxonomy and metadata.
- [x] Add failing contract tests for metric-vector completeness and runtime evidence.
- [x] Add failing contract tests for deterministic report serialization and CLI behavior.
- [x] Run focused tests, full architecture tests, and the full suite.
- [x] Produce a machine-readable before/after-compatible evidence report.
- [x] Review scope and run compileall, Ruff, strict OpenSpec validation, and diff check.

## Promotion

- [ ] Commit exactly the Campaign 1 scope as `Build Clean2D benchmark and observability foundation`.
- [ ] Continue to Campaign 2 only after Campaign 1 gates pass and its commit is verified.
