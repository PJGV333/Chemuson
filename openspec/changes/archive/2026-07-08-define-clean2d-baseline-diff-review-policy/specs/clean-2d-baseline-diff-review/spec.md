## ADDED Requirements

### Requirement: Stable Diff Review Policy
The system SHALL classify Clean 2D baseline report diffs using a stable review policy.

#### Scenario: Diff is classified
- **GIVEN** a Clean 2D baseline report diff
- **WHEN** the review policy is applied
- **THEN** each changed field receives a stable category and severity.

### Requirement: Difference Categories
The system SHALL distinguish contract changes, routing changes, geometry-diagnostic changes, snapshot changes and policy-evidence changes.

#### Scenario: Field categories are assigned
- **GIVEN** changed fields in a baseline diff
- **WHEN** review items are generated
- **THEN** result-state and reason changes are contract changes, source changes are routing changes, metric changes are geometry-diagnostic changes, snapshot changes are diagnostic changes, and policy-evidence changes are complex-policy changes.

### Requirement: Result-State Risk
The system SHALL classify result-state changes as contract-risk.

#### Scenario: Result state changes
- **GIVEN** a baseline diff contains a result-state change
- **WHEN** it is classified
- **THEN** the review item has severity `contract-risk`.

### Requirement: Metric Changes Are Observational
The system SHALL classify metric-only changes as geometry-diagnostic and observational-only.

#### Scenario: Metric differs outside tolerance
- **GIVEN** a baseline diff contains only a metric change
- **WHEN** it is classified
- **THEN** the review item has category `geometry-diagnostic` and remains observational-only.

### Requirement: Known-Delicate Visibility
The system SHALL keep known-delicate changes visible while labeling them separately.

#### Scenario: Known-delicate case changes
- **GIVEN** a changed case has the `known_delicate` tag
- **WHEN** review items are generated
- **THEN** the items remain present and are marked as known-delicate related.

### Requirement: JSON-Stable Review Summary
The system SHALL produce JSON-stable review summaries.

#### Scenario: Summary serializes
- **GIVEN** review items are summarized
- **WHEN** `json.dumps(..., allow_nan=False, sort_keys=True)` is called
- **THEN** serialization succeeds without custom encoders.

### Requirement: No Gates In This Change
The system SHALL NOT promote review classifications into CI gates in this change.

#### Scenario: Contract risk is reported
- **GIVEN** a review item has severity `contract-risk`
- **WHEN** the review summary is produced
- **THEN** it reports the risk observationally without enforcing a gate.

### Requirement: No Production Behaviour Changes
The system SHALL NOT change layout, ranking, selection, UI, backends or public document format.

#### Scenario: Review policy runs
- **GIVEN** the review policy classifies a diff
- **WHEN** the classification completes
- **THEN** baseline reports and Clean 2D production behaviour are not modified.
