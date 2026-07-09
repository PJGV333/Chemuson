# clean-2d-diff-blocking-policy Specification

## Purpose
TBD - created by archiving change define-clean2d-diff-blocking-policy. Update Purpose after archive.
## Requirements
### Requirement: Stable Clean 2D diff blocking policy evaluation
The system SHALL evaluate Clean 2D diff review summaries using a stable blocking policy.

#### Scenario: Empty review summary is clean
- **WHEN** the blocking policy evaluates a Clean 2D diff review summary with no review items
- **THEN** the policy decision SHALL be `clean`

### Requirement: Contract risk classification
The system SHALL classify `contract-risk` items as blocking candidates for human review.

#### Scenario: Contract risk item is blocking candidate
- **WHEN** the blocking policy evaluates a review summary containing a non-known-delicate `contract-risk` item
- **THEN** the policy decision SHALL be `blocking-candidate`

### Requirement: Needs review classification
The system SHALL classify routing, snapshot, and policy-evidence changes as `review-required`.

#### Scenario: Needs review item requires review
- **WHEN** the blocking policy evaluates a review summary containing a `needs-review` item for routing, snapshot, or policy evidence
- **THEN** the policy decision SHALL be `review-required`

### Requirement: Geometry risk remains observational
The system SHALL classify `geometry-risk` as `review-required` and observational-only.

#### Scenario: Geometry risk is review required and observational
- **WHEN** the blocking policy evaluates a review summary containing a `geometry-risk` item
- **THEN** the policy decision SHALL be `review-required`
- **AND** the policy decision SHALL indicate `observational_only`

### Requirement: Known delicate changes remain visible
The system SHALL keep known-delicate changes visible while allowing them to receive an `allowed-known-delicate` decision.

#### Scenario: Only known delicate changes are allowed known delicate
- **WHEN** the blocking policy evaluates a review summary containing only `known-delicate-change` items
- **THEN** the policy decision SHALL be `allowed-known-delicate`
- **AND** the policy output SHALL retain reasons for the known-delicate items

#### Scenario: Known delicate changes do not hide stronger risks
- **WHEN** the blocking policy evaluates a review summary containing both `known-delicate-change` and non-known-delicate `contract-risk` items
- **THEN** the policy decision SHALL be `blocking-candidate`
- **AND** the policy output SHALL retain reasons for both items

### Requirement: JSON-stable policy decisions
The system SHALL produce JSON-stable blocking policy decisions.

#### Scenario: Policy decision serializes deterministically
- **WHEN** the blocking policy decision is converted to a JSON-compatible mapping
- **THEN** it SHALL serialize with `json.dumps(..., allow_nan=False, sort_keys=True)`

### Requirement: Advisory-only policy
The system SHALL NOT enforce CI gates in this change.

#### Scenario: Blocking candidate does not enforce a gate
- **WHEN** the blocking policy returns `blocking-candidate`
- **THEN** it SHALL remain an advisory recommendation for human review
- **AND** it SHALL NOT fail CI or otherwise block execution automatically

### Requirement: No Clean 2D production behavior changes
The system SHALL NOT change Clean 2D production behaviour, layout, ranking, selection, UI, backends, or public document format.

#### Scenario: Policy runs only in test infrastructure
- **WHEN** the blocking policy helper is added
- **THEN** it SHALL live in test-only infrastructure
- **AND** it SHALL NOT modify production Clean 2D modules
