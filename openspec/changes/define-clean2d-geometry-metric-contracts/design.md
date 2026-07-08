# Design: Clean 2D Geometry Metric Contracts

## Overview
This change defines metric semantics for test and diagnostic reporting. It does not promote metrics into pass/fail geometry gates. The implementation should live in test-only corpus helpers unless a small reporting helper already exists and can be extended without changing production behaviour.

## Metric Definition Registry
A small registry will define one `Clean2DMetricDefinition` per stable metric. Each definition includes:
- `name`: stable serialized metric name.
- `value_type`: `string`, `integer`, `number`, or `nullable-number`.
- `unit`: physical or normalized unit, or `None` for categorical/informational values.
- `polarity`: `lower-is-better`, `higher-is-better`, `categorical`, or `informational`.
- `scope`: whole molecule, selected subgraph, candidate, or final result.
- `required`: whether the metric key is expected in corpus records.
- `diagnostic_only`: true for all metrics in this change.
- `optional_when`: short explanation of when `None` is allowed.

## Metric Classes
### A. Contract State Metrics
- `quality_class`: categorical state reported by existing quality helper.
- `reason`: informational reason string; may be empty or `None` when no reason applies.

### B. Topological Safety Metrics
- `crossings`: integer count of bond crossings in the evaluated scope.

### C. Collision / Spacing Metrics
- `min_nonbonded_distance`: nullable number in canvas pixels. `None` means no finite non-bonded distance was applicable.

### D. Ring-Shape Metrics
- `min_ring_degeneracy`: nullable dimensionless number. `None` means no ring was applicable.

### E. Bond-Length Metrics
- `length_rms_error`: number, normalized to target bond length.
- `length_max_error`: number, normalized to target bond length.

### F. Angle Metrics
- `angle_rms_deviation`: number in degrees.
- `angle_max_deviation`: number in degrees.

### G. Aggregate Visual Diagnostic
- `visual_score`: number where higher values are better. It remains diagnostic-only.

## Before/After Records
Corpus metric records should use this shape:

```python
{
    "before": {"quality_class": "...", ...},
    "after": {"quality_class": "...", ...} | None,
}
```

`after` may be `None` when Clean 2D returns no selected candidate. This is a controlled execution outcome, not a metric failure.

## Optional, Absent, NaN, And Non-Applicable Values
Metric keys should be present when the metric is required by the registry. A value may be `None` only when the definition allows a non-applicable or unavailable state. `NaN`, positive infinity, and negative infinity must not be serialized; test helpers normalize them to `None` only for metrics whose definitions allow nullable values.

## Floating-Point Tolerances
Metric comparison tests will use explicit tolerances for numeric values:
- Absolute tolerance: `1e-9` for normalized scores/errors.
- Pixel tolerance: `1e-6` for canvas-coordinate distances.
- Degree tolerance: `1e-6` for angle values.

These tolerances are for comparison and serialization stability only. They are not aesthetic thresholds.

## Diagnostic-Only Status
All metrics in this change are diagnostic-only. They may be used to inspect before/after behaviour, but they must not fail a regression case solely because geometry is visually imperfect. Future changes may promote selected metrics into gates only after a separate OpenSpec proposal defines thresholds and migration policy.

## Future Metric Table
Future metrics may include displacement, mean bond length, non-bonded collision counts, ring area ratio, candidate novelty, and backend provenance. This change documents that roadmap without implementing those metrics.
