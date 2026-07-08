## ADDED Requirements

### Requirement: Stable Metric Names
The system SHALL define stable Clean 2D geometry metric names.

#### Scenario: Captured metric has definition
- **GIVEN** the regression corpus emits a metric key
- **WHEN** metric contract validation runs
- **THEN** the key is present in the stable Clean 2D geometry metric registry.

### Requirement: Metric Semantics
The system SHALL define type, unit, polarity and availability semantics for each metric.

#### Scenario: Metric definition is complete
- **GIVEN** a registered geometry metric definition
- **WHEN** the definition is inspected
- **THEN** it declares value type, unit, polarity, scope, required status, optional availability, and diagnostic-only status.

### Requirement: Diagnostic Metrics Are Not Gates
The system SHALL distinguish diagnostic metrics from pass/fail geometric gates.

#### Scenario: Diagnostic metric is visually imperfect
- **GIVEN** a regression corpus case produces visually imperfect diagnostic metric values
- **WHEN** the metric contract test evaluates the case
- **THEN** the test does not fail solely because the diagnostic geometry value is aesthetically poor.

### Requirement: Declared Optionality
The system SHALL allow metrics to be absent only when declared non-applicable or unavailable.

#### Scenario: Optional metric is unavailable
- **GIVEN** a metric value is unavailable for a controlled result
- **WHEN** the metric record is serialized
- **THEN** the value is represented as `None` only if the metric definition allows nullable or unavailable values.

### Requirement: JSON-Stable Metric Records
The system SHALL serialize metric records using JSON-stable primitive values.

#### Scenario: Metric record serializes
- **GIVEN** a corpus case emits before/after metric records
- **WHEN** `json.dumps(..., allow_nan=False)` is called
- **THEN** serialization succeeds without custom encoders, `NaN`, or infinity.

### Requirement: Explicit Float Tolerances
The system SHALL compare floating-point metric values using explicit tolerances.

#### Scenario: Numeric metric comparison uses tolerance
- **GIVEN** two numeric metric records differ by less than the declared tolerance
- **WHEN** the metric comparison helper evaluates them
- **THEN** the values compare equal for regression stability purposes.

### Requirement: No Aesthetic Failures In This Change
The system SHALL NOT fail regression corpus cases solely because diagnostic geometry metrics are visually imperfect.

#### Scenario: Current baseline has poor geometry metric
- **GIVEN** a known-delicate corpus case has poor current geometry metrics
- **WHEN** the metric contract tests run
- **THEN** the case may pass if execution is controlled, identity invariants hold, and metrics serialize correctly.

### Requirement: Result State Independence
The system SHALL preserve the existing Clean 2D result-state contract independently from geometry metric values.

#### Scenario: Metric values do not change result state
- **GIVEN** metric contract helpers validate a corpus metric record
- **WHEN** the Clean 2D result state is inspected
- **THEN** the state remains the one returned by the existing Clean 2D engine contract and is not recomputed from diagnostic metric values.
