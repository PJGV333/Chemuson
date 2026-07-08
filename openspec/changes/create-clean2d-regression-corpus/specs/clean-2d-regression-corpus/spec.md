## ADDED Requirements

### Requirement: Stable Named Regression Cases
The corpus SHALL provide stable named regression cases.

#### Scenario: Registry exposes stable names
- **GIVEN** the regression corpus registry is loaded
- **WHEN** the cases are enumerated
- **THEN** every case has a non-empty unique name suitable for pytest ids and debug snapshots.

### Requirement: Case Metadata Contract
Each case SHALL declare family, mode, target and expected contract states.

#### Scenario: Case metadata is complete
- **GIVEN** a regression case is registered
- **WHEN** its metadata is inspected
- **THEN** it declares a non-empty family, a Clean 2D mode, a target descriptor, and at least one expected contract state.

### Requirement: Valid Molecule Graph Builders
Each case SHALL build a valid molecule graph.

#### Scenario: Builder returns a graph executable by Clean 2D
- **GIVEN** a regression case builder
- **WHEN** the builder is executed
- **THEN** it returns a molecule graph with atoms and bonds that can be passed to the Clean 2D engine.

### Requirement: Controlled Clean 2D Execution
Each case SHALL execute through the Clean 2D engine without uncontrolled exceptions.

#### Scenario: Corpus case runs through engine
- **GIVEN** a valid regression case
- **WHEN** Clean 2D is executed for the declared mode and target
- **THEN** execution either returns a controlled result state declared by the case or records a declared known controlled failure.

### Requirement: Chemical Identity Preservation
Each case SHALL preserve chemical identity invariants.

#### Scenario: Identity signature is preserved
- **GIVEN** a regression case molecule before Clean 2D execution
- **WHEN** Clean 2D returns a result
- **THEN** atom ids, bond ids, element symbols, bond endpoints, bond orders, and stable stereo metadata present in the input are preserved.

### Requirement: JSON-Serializable Debug Snapshots
Each case SHALL be able to produce a JSON-serializable debug snapshot.

#### Scenario: Snapshot can be serialized
- **GIVEN** a regression case has been executed
- **WHEN** its debug snapshot is generated
- **THEN** `json.dumps` can serialize the snapshot without custom encoders.

### Requirement: Known-Delicate Case Recording
Known-delicate cases SHALL be allowed to record current behavior without forcing geometry fixes.

#### Scenario: Known-delicate case captures current state
- **GIVEN** a case is marked known-delicate
- **WHEN** the regression test executes it
- **THEN** the test verifies controlled execution, identity preservation, and snapshot serializability without requiring a geometry improvement.

### Requirement: No Behavioural Geometry Changes
The corpus SHALL NOT change candidate ranking, geometry generation, UI behavior or document format.

#### Scenario: Corpus remains test-only infrastructure
- **GIVEN** the regression corpus is added
- **WHEN** Clean 2D runs outside the tests
- **THEN** candidate ranking, geometry generation, UI behavior, and public document format remain unchanged.
