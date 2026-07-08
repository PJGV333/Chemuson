# clean-2d-debug-snapshots Specification

## Purpose
TBD - created by archiving change add-clean2d-debug-snapshots. Update Purpose after archive.
## Requirements
### Requirement: Opt-in Clean 2D debug snapshot capture
The system SHALL provide Clean 2D debug snapshot capture only when explicitly enabled by a supported opt-in mechanism.

#### Scenario: Normal Clean 2D use does not capture snapshots
- **GIVEN** Clean 2D is executed without the debug environment variable, without an explicit debug parameter, and without a test snapshot helper
- **WHEN** the Clean 2D run completes
- **THEN** the system SHALL NOT write or return a debug snapshot
- **AND** the Clean 2D final state, stable reason, candidate ranking, selected candidate, and geometry SHALL match the behavior of the same run before this capability existed

#### Scenario: Environment variable enables snapshots
- **GIVEN** the Clean 2D debug snapshot environment variable is enabled
- **WHEN** Clean 2D executes
- **THEN** the system SHALL capture a Clean 2D debug snapshot for that run

#### Scenario: Explicit debug parameter enables snapshots
- **GIVEN** a Clean 2D caller passes the explicit debug snapshot parameter
- **WHEN** Clean 2D executes
- **THEN** the system SHALL capture a Clean 2D debug snapshot for that run

#### Scenario: Test helper enables snapshots
- **GIVEN** a test uses the Clean 2D snapshot helper
- **WHEN** Clean 2D executes inside the helper scope
- **THEN** the system SHALL capture a Clean 2D debug snapshot for the test-controlled location or return path

### Requirement: Stable JSON snapshot schema
The system SHALL serialize each Clean 2D debug snapshot as stable JSON containing a schema identifier and schema version.

#### Scenario: Snapshot is JSON serializable
- **GIVEN** Clean 2D debug snapshot capture is enabled
- **WHEN** a snapshot is produced
- **THEN** the snapshot SHALL be serializable using standard JSON primitives
- **AND** the serialized JSON SHALL NOT require Python object reprs, Qt objects, RDKit objects, memory addresses, or machine-specific paths

#### Scenario: Snapshot declares schema and version
- **GIVEN** a Clean 2D debug snapshot exists
- **WHEN** the snapshot is read by test code
- **THEN** the snapshot SHALL contain a stable schema identifier
- **AND** the snapshot SHALL contain a schema version suitable for future compatibility checks

### Requirement: Snapshot records Clean 2D run context
The system SHALL record the Clean 2D run context needed to reproduce the target of the run.

#### Scenario: Snapshot records mode and target
- **GIVEN** Clean 2D debug snapshot capture is enabled
- **WHEN** Clean 2D is executed for a whole structure
- **THEN** the snapshot SHALL record the Clean 2D mode
- **AND** the snapshot SHALL record that the target was the whole structure

#### Scenario: Snapshot records targeted atoms
- **GIVEN** Clean 2D debug snapshot capture is enabled
- **WHEN** Clean 2D is executed for a subset of atoms
- **THEN** the snapshot SHALL record the Clean 2D mode
- **AND** the snapshot SHALL record the target atom IDs
- **AND** the target atom IDs SHALL be stable JSON values

#### Scenario: Snapshot records initial selection
- **GIVEN** Clean 2D debug snapshot capture is enabled
- **WHEN** Clean 2D starts with an existing atom or bond selection
- **THEN** the snapshot SHALL record the initial selection separately from the target atom IDs

### Requirement: Snapshot records molecule topology and coordinates
The system SHALL record the molecule topology and relevant coordinates used by the Clean 2D run.

#### Scenario: Snapshot records atoms and bonds
- **GIVEN** Clean 2D debug snapshot capture is enabled
- **WHEN** a snapshot is produced
- **THEN** the snapshot SHALL include atom IDs
- **AND** the snapshot SHALL include bond IDs
- **AND** the snapshot SHALL include bond connectivity by atom ID
- **AND** the snapshot SHALL include bond orders

#### Scenario: Snapshot records initial coordinates
- **GIVEN** Clean 2D debug snapshot capture is enabled
- **WHEN** Clean 2D starts
- **THEN** the snapshot SHALL record initial coordinates for atoms participating in the run context

#### Scenario: Snapshot records final coordinates when available
- **GIVEN** Clean 2D debug snapshot capture is enabled
- **WHEN** Clean 2D produces final coordinates
- **THEN** the snapshot SHALL record final coordinates

#### Scenario: Snapshot omits final coordinates when unavailable
- **GIVEN** Clean 2D debug snapshot capture is enabled
- **WHEN** Clean 2D ends before final coordinates exist
- **THEN** the snapshot SHALL remain valid JSON
- **AND** the snapshot SHALL represent final coordinates as absent or null according to the snapshot schema

### Requirement: Snapshot records candidate evaluation evidence
The system SHALL record candidate source and diagnostic evidence observed during Clean 2D candidate evaluation without changing the evaluation process.

#### Scenario: Snapshot records candidate sources
- **GIVEN** Clean 2D debug snapshot capture is enabled
- **WHEN** Clean 2D evaluates candidate sources
- **THEN** the snapshot SHALL record the candidate sources that were evaluated
- **AND** recording candidate sources SHALL NOT add, remove, reorder, or re-score candidates

#### Scenario: Snapshot records quality diagnostics when present
- **GIVEN** Clean 2D debug snapshot capture is enabled
- **WHEN** a candidate has `quality_diagnostic` data
- **THEN** the snapshot SHALL include that candidate's `quality_diagnostic` data in JSON-serializable form

#### Scenario: Snapshot tolerates missing quality diagnostics
- **GIVEN** Clean 2D debug snapshot capture is enabled
- **WHEN** a candidate has no `quality_diagnostic` data
- **THEN** the snapshot SHALL remain valid
- **AND** the snapshot SHALL still preserve the candidate source information when available

### Requirement: Snapshot records final Clean 2D decision
The system SHALL record the final Clean 2D decision state and stable reason when available.

#### Scenario: Snapshot records final state
- **GIVEN** Clean 2D debug snapshot capture is enabled
- **WHEN** Clean 2D completes with any final state
- **THEN** the snapshot SHALL record the final state using the stable Clean 2D decision contract

#### Scenario: Snapshot records stable final reason when present
- **GIVEN** Clean 2D debug snapshot capture is enabled
- **WHEN** Clean 2D completes with a stable final reason
- **THEN** the snapshot SHALL record the stable final reason

#### Scenario: Snapshot remains valid without final reason
- **GIVEN** Clean 2D debug snapshot capture is enabled
- **WHEN** Clean 2D completes without a stable final reason
- **THEN** the snapshot SHALL remain valid JSON
- **AND** the snapshot SHALL represent the final reason as absent or null according to the snapshot schema

### Requirement: Optional internal metadata
The system SHALL allow optional internal metadata in Clean 2D debug snapshots without making that metadata required for tests or normal behavior.

#### Scenario: Snapshot includes optional metadata
- **GIVEN** Clean 2D debug snapshot capture is enabled
- **WHEN** implementation-specific metadata is available
- **THEN** the snapshot MAY include internal metadata only in a JSON-serializable metadata section
- **AND** tests SHALL NOT depend on unstable metadata such as object reprs, memory addresses, timestamps, or machine-specific paths

### Requirement: Snapshot test helpers
The system SHALL provide helpers for writing and reading Clean 2D debug snapshots in tests.

#### Scenario: Test writes snapshot
- **GIVEN** a test enables Clean 2D debug snapshot capture through the test helper
- **WHEN** Clean 2D completes
- **THEN** the helper SHALL write a JSON snapshot to a test-controlled path or return an equivalent JSON-serializable object

#### Scenario: Test reads snapshot
- **GIVEN** a test has a Clean 2D debug snapshot JSON file
- **WHEN** the test reads the snapshot through the helper
- **THEN** the helper SHALL return the parsed snapshot data
- **AND** the helper SHALL preserve the stable schema fields required by this specification

#### Scenario: Snapshot tests do not assert algorithm fixes
- **GIVEN** snapshot tests cover known complex Clean 2D failure shapes
- **WHEN** those tests validate snapshots
- **THEN** the tests SHALL assert snapshot structure and captured evidence
- **AND** the tests SHALL NOT require geometry, ranking, selection, or algorithmic fixes for the known out-of-scope failures

