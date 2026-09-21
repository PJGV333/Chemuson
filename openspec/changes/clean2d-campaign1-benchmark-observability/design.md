# Design: Clean2D Campaign 1 Benchmark and Observability

## Existing surfaces

Campaign 1 reuses `tests.clean2d_regression.cases`, `metrics`, `assertions`, `baselines`, `reports`, and `review_policy`, plus `scripts/clean2d_baseline_report.py`. Production debug snapshots remain the source of candidate and decision evidence; the benchmark layer only canonicalizes and aggregates data already exposed by those surfaces.

## Case identity and taxonomy

`Clean2DRegressionCase` remains the stable case identity. A case declares a size class explicitly, while topology metadata is derived from its built graph at execution time. Family tags remain orthogonal to size class. Existing tags and names are preserved.

Derived metadata is observational and includes atom count, heavy-atom count, bond count, ring count, connected components, and conservative connector/block counts when derivable. Missing or unsupported values are represented as `None`, never guessed.

## Metric vector

The report retains the existing geometry metric registry and adds campaign-level evidence in the same baseline record rather than creating a second telemetry system. The vector contains safety/quality metrics when computable, candidate count/sources, result state, stable reason, and runtime evidence. Values are JSON primitives; unavailable values are explicit `None`; no new acceptance threshold is introduced.

## Determinism

Case ordering, candidate source ordering, metadata ordering, float canonicalization, and report serialization are stable. Runtime is evidence only and is excluded from equivalence decisions because wall-clock timing is inherently ephemeral. Repeated executions must otherwise compare equivalent.

## CLI

The existing developer-only CLI remains the entrypoint. `write` emits one canonical report, `compare` emits a machine-readable diff and returns `0` for equivalent or `1` for changed reports, and `review` classifies observable changes without silently accepting them. Errors return `2`.

## Scope guard

Allowed paths are this OpenSpec, test-owned regression infrastructure, the developer baseline-report script, and focused tests. No production Clean2D file or architecture catalog is changed by Campaign 1.
