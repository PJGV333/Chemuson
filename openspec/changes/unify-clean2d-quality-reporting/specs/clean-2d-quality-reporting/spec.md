## ADDED Requirements

### Requirement: Diagnostic Contract Schema
The system SHALL expose a single stable contract shape for Clean 2D reporting diagnostics. The contract defines the following fields:
- **state** (string): one of `applied`, `rejected`, `preserve-only`, `no-op`, or `failed-controlled`.
- **reason** (string, optional): one of the stable reasons defined in document-clean2d-decision-contract (`invalid-coordinates`, `invariant-violation`, `stereo-risk`, `boundary-bond-risk`, `new-crossing-risk`, `collision-risk`, `collapsed-ring-risk`, `excessive-displacement`, `worse-quality`, `backend-failure`). The field is required for `rejected` and `failed-controlled`, optional for `no-op` and `preserve-only`, and absent for `applied` unless provided by an explicit non-contract diagnostic path.
- **score** (float, 0-1): a normalized reporting score where 1 is best and 0 indicates failure. This value is for reporting and diagnostics only and SHALL NOT change ranking, candidate selection, or layout decisions.
- **source** (string): identifier of the originating module (`engine`, `safety`, `depictor`, `controller`).
- **internal** (object, optional): module-specific metadata that may contain raw diagnostics or helper data. Tests should ignore this field unless explicitly testing internal logic.

#### Scenario: Successful diagnostic
- **GIVEN** a molecule is cleaned and passes all checks
- **WHEN** the diagnostic is generated
- **THEN** the JSON includes `state=applied`, no `reason`, a positive `score` between 0 and 1, and a non-empty `source`.

#### Scenario: Rejected with reason
- **GIVEN** a molecule fails safety checks due to stereochemistry risk
- **WHEN** the diagnostic is generated
- **THEN** the JSON includes `state=rejected`, `reason=stereo-risk`, a low `score` near 0, and `source=safety`.

#### Scenario: Internal metadata preserved
- **GIVEN** a module produces additional internal data during diagnosis
- **WHEN** the diagnostic is generated
- **THEN** the JSON includes an `internal` field containing that data without affecting the stable contract fields.

### Requirement: Score Normalization
The system SHALL normalize reporting quality scores to the 0-1 range regardless of source module's original diagnostic scale. This normalized score SHALL be used only for reporting and diagnostics in this change.

#### Scenario: Engine score normalization
- **GIVEN** engine produces a raw quality metric of 85 (on its internal scale)
- **WHEN** the diagnostic is generated
- **THEN** the `score` field contains a normalized value between 0 and 1, not the raw 85.

#### Scenario: Safety score normalization
- **GIVEN** safety produces a penalty-based metric of -2 (on its internal scale)
- **WHEN** the diagnostic is generated
- **THEN** the `score` field contains a normalized value between 0 and 1, not the raw -2.

#### Scenario: Reporting score does not affect selection
- **GIVEN** existing candidate ranking and selection logic uses its current internal score inputs
- **WHEN** a normalized reporting score is added to diagnostics
- **THEN** ranking and candidate selection continue to use the existing internal inputs unchanged.

### Requirement: Source Identification
The system SHALL include a `source` field identifying which module produced the diagnostic, enabling traceability in multi-module pipelines.

#### Scenario: Engine source identification
- **GIVEN** engine generates a diagnostic for a cleaned molecule
- **WHEN** the diagnostic is generated
- **THEN** the `source` field equals `"engine"`.

#### Scenario: Safety source identification
- **GIVEN** safety module detects a collision risk and rejects the candidate
- **WHEN** the diagnostic is generated
- **THEN** the `source` field equals `"safety"` and `reason` equals `"collision-risk"`.

### Requirement: Contract Test Validation
The system SHALL provide contract tests that validate diagnostics conform to the unified contract shape, regardless of source module. Tests SHALL verify state values, reason validity when required, reporting score range, and source identification.

#### Scenario: Contract test validates applied diagnostic
- **GIVEN** a diagnostic with `state=applied` is produced by any module
- **WHEN** the contract test runs
- **THEN** it passes without requiring `reason`, verifies `score` is between 0 and 1, and confirms `source` is non-empty.

#### Scenario: Contract test validates rejected diagnostic
- **GIVEN** a diagnostic with `state=rejected` is produced by any module
- **WHEN** the contract test runs
- **THEN** it fails if `reason` is missing or not one of the stable reasons, and verifies `score` is between 0 and 1.

#### Scenario: Contract test validates failed-controlled diagnostic
- **GIVEN** a diagnostic with `state=failed-controlled` is produced by any module
- **WHEN** the contract test runs
- **THEN** it fails if `reason` is missing or not one of the stable reasons, and verifies `score` is between 0 and 1.

#### Scenario: Contract test validates optional reasons
- **GIVEN** a diagnostic with `state=no-op` or `state=preserve-only`
- **WHEN** the contract test runs
- **THEN** it passes whether `reason` is absent or contains one of the stable reasons.

#### Scenario: Contract test ignores internal metadata
- **GIVEN** a diagnostic includes arbitrary data in the `internal` field
- **WHEN** the contract test runs
- **THEN** it passes without validating the contents of `internal`, treating it as opaque module-specific data.

### Requirement: Minimal Adapter Scope
The system SHALL add adapters or wrappers only where a clear diagnostic surface already exists. This change MUST NOT perform a full migration of engine.py, safety.py, depiction_candidates.py, or Clean2DController to a new internal representation.

#### Scenario: Existing diagnostic surface adapted
- **GIVEN** a module already exposes diagnostic information at a clear boundary
- **WHEN** the reporting contract is connected
- **THEN** a minimal adapter maps that information to the stable contract without changing layout, geometry, ranking, or selection behavior.

#### Scenario: No clear diagnostic surface
- **GIVEN** a module does not expose a clear diagnostic boundary
- **WHEN** the reporting contract is connected
- **THEN** the implementation records the gap as follow-up work instead of refactoring the module in this change.
