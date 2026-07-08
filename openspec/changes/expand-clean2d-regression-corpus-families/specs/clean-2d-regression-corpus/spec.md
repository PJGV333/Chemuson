## ADDED Requirements

### Requirement: Expanded Representative Families
The regression corpus SHALL include representative simple, aromatic, fused, heteroaromatic, charged, stereo-sensitive, selection-boundary and multi-block cases.

#### Scenario: Minimum families are represented
- **GIVEN** the regression corpus registry is loaded
- **WHEN** family labels and tags are inspected
- **THEN** the registry includes cases covering simple, aromatic, fused, heteroaromatic, charged, stereo-sensitive, selection-boundary, and multi-block inputs.

### Requirement: Stable Family Labels
Each case SHALL declare a stable family label.

#### Scenario: Case has family label
- **GIVEN** a registered corpus case
- **WHEN** its metadata is inspected
- **THEN** it has a non-empty stable family label.

### Requirement: Stable Case Classification Tags
Each case SHALL declare whether it is baseline, known-delicate, known-failure, stereo-sensitive, selection-boundary or complex-policy-sensitive.

#### Scenario: Case classification is explicit
- **GIVEN** a registered corpus case
- **WHEN** its metadata is inspected
- **THEN** it has at least one stable classification tag.

### Requirement: Identity Invariants Preserved
Each case SHALL preserve chemical identity invariants after Clean 2D execution.

#### Scenario: Expanded case preserves identity
- **GIVEN** an expanded corpus case
- **WHEN** Clean 2D executes with the declared mode and target
- **THEN** atom ids, bond ids, element symbols, bond endpoints, bond orders, aromaticity, charges, and stable stereo metadata remain unchanged when a candidate is selected.

### Requirement: Stereo Metadata Preservation
Each stereo-sensitive case SHALL preserve stereo metadata, even when geometry is not improved.

#### Scenario: Stereo-sensitive case preserves metadata
- **GIVEN** a stereo-sensitive corpus case
- **WHEN** Clean 2D executes
- **THEN** atom and bond stereo metadata present before execution is still present after execution.

### Requirement: Formal Charge Preservation
Each charged case SHALL preserve formal charges.

#### Scenario: Charged case preserves charges
- **GIVEN** a charged corpus case
- **WHEN** Clean 2D executes
- **THEN** formal charges on atoms are unchanged.

### Requirement: Selection Boundary Preservation
Each selection-boundary case SHALL preserve boundary connectivity and not move non-target atoms unexpectedly beyond the current contract.

#### Scenario: Selection-boundary case preserves boundary bonds
- **GIVEN** a selection-boundary corpus case with target atoms
- **WHEN** Clean 2D executes
- **THEN** bonds connecting target and non-target atoms keep their ids, endpoints, orders, aromaticity, and stereo metadata, and non-target atoms are not asserted as layout targets.

### Requirement: Complex Policy Guard Recording
Each complex-policy-sensitive case SHALL not require global redraw when the existing policy classifies it as high-risk.

#### Scenario: Complex policy guard remains observational
- **GIVEN** a complex-policy-sensitive corpus case
- **WHEN** the current Clean 2D policy returns preserve-only or failed-controlled behaviour
- **THEN** the corpus records that controlled state without requiring an applied global redraw.

### Requirement: Observational Corpus Only
The expanded corpus SHALL remain observational and SHALL NOT impose new aesthetic thresholds.

#### Scenario: Metrics do not become aesthetic gates
- **GIVEN** an expanded corpus case emits diagnostic metrics
- **WHEN** regression tests evaluate the case
- **THEN** metrics are captured for diagnostics but do not fail the test solely because current geometry is visually imperfect.
