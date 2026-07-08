# clean-2d-deterministic-baselines Specification

## Purpose
TBD - created by archiving change stabilize-clean2d-deterministic-baselines. Update Purpose after archive.
## Requirements
### Requirement: JSON-Stable Baseline Records
The system SHALL produce JSON-stable Clean 2D baseline records for regression corpus cases.

#### Scenario: Corpus case emits baseline record
- **GIVEN** a Clean 2D regression corpus case
- **WHEN** the deterministic baseline harness executes the case
- **THEN** the resulting baseline record serializes with `json.dumps(..., allow_nan=False)` without custom encoders.

### Requirement: Canonical Baseline Comparison
The system SHALL canonicalize baseline records before comparison.

#### Scenario: Equivalent records canonicalize consistently
- **GIVEN** two baseline records from equivalent executions
- **WHEN** canonicalization is applied
- **THEN** non-semantic ordering differences do not affect comparison.

### Requirement: Explicit Metric Tolerances
The system SHALL compare repeated executions of the same case using explicit metric tolerances.

#### Scenario: Float metric differs within tolerance
- **GIVEN** two numeric diagnostic metric values differ within the declared tolerance
- **WHEN** baseline comparison runs
- **THEN** the values compare equivalent.

### Requirement: Result-State Contract Preservation
The system SHALL preserve the existing Clean 2D result-state contract and SHALL NOT recompute result state from metrics.

#### Scenario: Baseline record keeps engine result state
- **GIVEN** Clean 2D returns a result state
- **WHEN** a baseline record is created
- **THEN** the record stores the returned state without deriving a replacement state from geometry metrics.

### Requirement: Diagnostic Metrics Remain Diagnostic
The system SHALL keep diagnostic metrics diagnostic-only.

#### Scenario: Visual metric is imperfect
- **GIVEN** a known-delicate case has visually imperfect metrics
- **WHEN** deterministic baseline tests run
- **THEN** the test does not fail solely because the current geometry is visually imperfect.

### Requirement: Observable Source Labels
The system SHALL record selected source and candidate source labels when available.

#### Scenario: Candidate sources are available
- **GIVEN** Clean 2D returns candidate source labels
- **WHEN** a baseline record is created
- **THEN** the selected source and observed candidate source sequence are recorded.

### Requirement: Snapshot Canonicalization
The system SHALL make snapshot comparison independent of non-semantic ordering and ephemeral fields.

#### Scenario: Snapshot contains ordered dictionaries
- **GIVEN** a Clean 2D debug snapshot
- **WHEN** the baseline harness canonicalizes it
- **THEN** dictionary key ordering does not affect comparison and documented ephemeral fields are excluded.

### Requirement: Known Instability Documentation
The system SHALL document known observational instability instead of hiding it.

#### Scenario: Stable result fields differ between executions
- **GIVEN** repeated executions differ in result state, stable reason, selected source, or candidate source ordering
- **WHEN** baseline comparison detects the difference
- **THEN** the case is marked as known observational instability or the difference is reported as follow-up.

### Requirement: No Behavioural Changes
The system SHALL NOT change layout, candidate ranking, candidate selection, UI, backends, or public document format.

#### Scenario: Baseline harness runs
- **GIVEN** the deterministic baseline harness executes a corpus case
- **WHEN** Clean 2D production code is invoked
- **THEN** existing production behaviour is observed but not changed by the harness.

