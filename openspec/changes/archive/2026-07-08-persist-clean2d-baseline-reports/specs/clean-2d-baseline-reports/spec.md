## ADDED Requirements

### Requirement: JSON-Stable Baseline Reports
The system SHALL produce JSON-stable Clean 2D baseline reports.

#### Scenario: Full corpus report serializes
- **GIVEN** the Clean 2D regression corpus
- **WHEN** a baseline report is generated
- **THEN** `json.dumps(..., allow_nan=False, sort_keys=True)` serializes the report without custom encoders.

### Requirement: One Record Per Case
The system SHALL include one canonical baseline record per regression corpus case.

#### Scenario: Corpus coverage is complete
- **GIVEN** the regression corpus contains registered cases
- **WHEN** a baseline report is generated
- **THEN** the report contains exactly one record for each registered case.

### Requirement: Stable Record Ordering
The system SHALL order records by stable case name.

#### Scenario: Records are sorted
- **GIVEN** a baseline report is generated
- **WHEN** record names are inspected
- **THEN** they are sorted lexicographically by case name.

### Requirement: Schema And Version Metadata
The system SHALL include schema and version metadata.

#### Scenario: Report metadata is present
- **GIVEN** a baseline report is generated
- **WHEN** its top-level fields are inspected
- **THEN** it includes a stable schema identifier and integer version.

### Requirement: Aggregate Summary
The system SHALL include an aggregate summary of case counts, result states and classification tags.

#### Scenario: Summary reflects records
- **GIVEN** a baseline report is generated
- **WHEN** its summary is inspected
- **THEN** case count, result-state counts, family counts, and tag counts match the records.

### Requirement: Report Comparison
The system SHALL compare two baseline reports using the deterministic baseline equivalence rules.

#### Scenario: Consecutive reports are equivalent
- **GIVEN** two reports generated consecutively from the same corpus
- **WHEN** they are compared
- **THEN** the diff reports them as equivalent.

### Requirement: Structured Difference Reporting
The system SHALL report differences in result state, reason, selected source, candidate sources, metrics, snapshot and policy evidence.

#### Scenario: Artificial report changes are detected
- **GIVEN** a baseline report is modified in an observable field
- **WHEN** report comparison runs
- **THEN** the diff identifies the changed case and field.

### Requirement: Observational Differences
The system SHALL keep baseline report differences observational in this change.

#### Scenario: Difference is detected
- **GIVEN** two baseline reports differ
- **WHEN** the diff is produced
- **THEN** the diff records the difference without promoting it into a geometry gate.

### Requirement: No Metric Gates
The system SHALL NOT promote diagnostic metric differences into geometry gates.

#### Scenario: Metric differs outside tolerance
- **GIVEN** a diagnostic metric differs outside tolerance
- **WHEN** reports are compared
- **THEN** the difference is reported observationally rather than failing as a geometry quality gate.

### Requirement: No Production Behaviour Changes
The system SHALL NOT change layout, ranking, selection, UI, backends or public document format.

#### Scenario: Report generation runs
- **GIVEN** the baseline report helper runs in tests
- **WHEN** Clean 2D production code is invoked
- **THEN** production behaviour is observed but not modified by report generation.
