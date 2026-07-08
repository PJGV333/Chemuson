# clean-2d-complex-policy Specification

## Purpose
TBD - created by archiving change stabilize-clean2d-complex-policy. Update Purpose after archive.
## Requirements
### Requirement: Deterministic complex-policy classification
The system SHALL provide a deterministic Clean 2D complex-policy classification for structures evaluated by `classify_clean2d_complexity` or equivalent complex-policy entry points.

#### Scenario: Classification uses existing structure signals
- **GIVEN** a molecule with existing graph, multilayer, and layout-quality signals
- **WHEN** Clean 2D classifies complexity for routing
- **THEN** the policy SHALL classify the structure using those existing signals
- **AND** the policy SHALL NOT generate Clean 2D candidates as part of classification

#### Scenario: Classification exposes stable route decisions
- **GIVEN** Clean 2D classifies any structure for complexity routing
- **WHEN** the classification is returned to the engine
- **THEN** the classification SHALL expose stable information indicating whether preserve-only is required
- **AND** the classification SHALL expose stable information indicating whether local repair is allowed
- **AND** the classification SHALL expose stable information indicating whether global redraw candidates are allowed
- **AND** the classification SHALL expose a stable reason or policy evidence for these route decisions

### Requirement: High-risk complex structures block global redraw candidates
The system SHALL classify structures with high-risk complexity signals as not eligible for global redraw candidates unless the policy explicitly marks that route safe.

#### Scenario: Hierarchical blocks block global redraw
- **GIVEN** a Clean 2D structure with hierarchical block-layout signals
- **WHEN** the complex policy classifies the structure
- **THEN** the policy SHALL mark global redraw candidates as disallowed
- **AND** Clean 2D SHALL NOT call global redraw candidate sources for that run

#### Scenario: Macrocycles and cyclophanes block global redraw by default
- **GIVEN** a Clean 2D structure classified as a macrocycle or cyclophane
- **WHEN** the complex policy classifies the structure
- **THEN** the policy SHALL mark global redraw candidates as disallowed by default
- **AND** Clean 2D SHALL preserve existing safety gates for any allowed alternative route

#### Scenario: Internal cavities and complex bridges block global redraw by default
- **GIVEN** a Clean 2D structure with internal cavity or complex bridged-system signals
- **WHEN** the complex policy classifies the structure
- **THEN** the policy SHALL mark global redraw candidates as disallowed by default
- **AND** Clean 2D SHALL NOT use global redraw as a fallback for that high-risk structure

### Requirement: Preserve-only is mandatory for high-risk structures without approved repair routes
The system SHALL use preserve-only behavior for high-risk complex structures when neither local repair nor internal non-redraw candidates are policy-approved.

#### Scenario: High-risk imported structure preserves geometry
- **GIVEN** a high-risk complex structure with existing coordinates
- **AND** the complex policy disallows global redraw
- **AND** the complex policy does not approve local repair or another non-redraw repair route
- **WHEN** Clean 2D executes in quick or publication mode
- **THEN** Clean 2D SHALL preserve the existing geometry
- **AND** Clean 2D SHALL expose a stable final state and stable reason for the preserve-only decision

#### Scenario: Preserve-only does not alter simple structures
- **GIVEN** a simple structure without high-risk complex-policy signals
- **WHEN** Clean 2D executes in quick or publication mode
- **THEN** the complex policy SHALL NOT force preserve-only
- **AND** Clean 2D simple-structure geometry behavior SHALL remain unchanged

### Requirement: Local repair eligibility is policy-gated
The system SHALL use `local_graph_cleaner` only when the complex policy explicitly allows local repair for the classified structure.

#### Scenario: Hierarchical blocks do not select local graph route
- **GIVEN** a tetrandrine-like structure with hierarchical block-layout signals
- **WHEN** Clean 2D executes in quick or publication mode
- **THEN** the complex policy SHALL reject local-graph routing unless a specific safe local defect route is approved
- **AND** Clean 2D SHALL NOT select `local_graph` as the final candidate source for that structure

#### Scenario: Local repair remains available for approved local defects
- **GIVEN** a complex structure with a local defect classified as safe for local repair
- **WHEN** Clean 2D executes in quick or publication mode
- **THEN** the complex policy MAY allow `local_graph_cleaner`
- **AND** any accepted local repair SHALL still pass existing Clean 2D safety gates

### Requirement: Internal candidates are allowed only when policy-approved
The system SHALL allow internal non-redraw candidates for complex structures only when the complex policy marks them safe for the classified structure.

#### Scenario: Policy-approved internal route may run
- **GIVEN** a complex structure where the policy disallows global redraw but approves a specific internal non-redraw route
- **WHEN** Clean 2D executes
- **THEN** Clean 2D MAY evaluate that internal route
- **AND** Clean 2D SHALL continue to apply existing invariant and safety validation before accepting it

#### Scenario: Unapproved internal route is skipped
- **GIVEN** a high-risk complex structure where the policy disallows an internal route
- **WHEN** Clean 2D executes
- **THEN** Clean 2D SHALL NOT use that route as a fallback around global redraw blocking

### Requirement: Policy decisions are observable through diagnostics
The system SHALL expose complex-policy routing decisions through stable diagnostic evidence without requiring debug snapshots for normal execution.

#### Scenario: Quality diagnostic includes policy evidence when available
- **GIVEN** Clean 2D classifies a complex structure and produces a result diagnostic
- **WHEN** `quality_diagnostic` is available for the result or selected candidate
- **THEN** the diagnostic SHALL include JSON-serializable policy evidence when available
- **AND** that evidence SHALL NOT alter ranking or selection

#### Scenario: Debug snapshot records policy evidence when enabled
- **GIVEN** Clean 2D debug snapshots are explicitly enabled
- **WHEN** Clean 2D executes for a complex-policy decision
- **THEN** the snapshot SHALL include available policy evidence in JSON-serializable metadata or diagnostics
- **AND** disabling snapshots SHALL NOT change the policy decision

### Requirement: In-scope regression coverage
The system SHALL include regression tests for the in-scope complex-policy failures without asserting fixes for out-of-scope fused-ring or propose-mode failures.

#### Scenario: Tetrandrine-like hierarchical blocks avoid local graph selection
- **GIVEN** the existing tetrandrine-like hierarchical blocks regression case
- **WHEN** Clean 2D executes after this change
- **THEN** the test SHALL verify that hierarchical complexity is classified as high risk
- **AND** the test SHALL verify that `local_graph` is not selected as the final candidate source

#### Scenario: Complex quick clean does not call global redraw candidates
- **GIVEN** the existing complex quick-clean regression case
- **WHEN** Clean 2D executes after this change
- **THEN** the test SHALL verify that global redraw candidate sources are not called
- **AND** the test SHALL verify that the final outcome follows the complex-policy route instead of global redraw fallback

#### Scenario: Out-of-scope failures remain out of scope
- **GIVEN** naphthalene fused aromatic systems or propose-mode smart-clean failures
- **WHEN** tests are added for this change
- **THEN** those tests SHALL NOT require fused-ring geometry fixes
- **AND** those tests SHALL NOT require propose-mode behavior fixes

