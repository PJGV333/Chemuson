# Special Projection Geometry

## ADDED Requirements

### Requirement: Projection geometry is renderer-style independent
Haworth, Fischer and chair projection builders SHALL expose deterministic geometry and typed projection-specific parameters independently of Qt widgets, dialogs, canvas/items, and global `DrawingStyle`. They MAY mark projection-specific front edges using molecular bond semantics. Conversion to `MolGraph` SHALL preserve deterministic atom/bond identity and SHALL NOT invent stereo metadata.

#### Scenario: Build projection without GUI dependencies
- **GIVEN** one projection-specific parameter value object
- **WHEN** its geometry builder is called repeatedly
- **THEN** it returns equivalent geometry and graphs
- **AND** importing the projection package does not load dialogs, canvas, graphics items, or Qt widgets.

### Requirement: Haworth geometry follows explicit projection rules
The Haworth builder SHALL create a stable ring with an explicit ring oxygen, canonical atom ordering, explicit front-edge policy, and substituent placement derived from declared projection parameters. Alpha/beta and D/L orientation SHALL be explicit inputs or documented template data; the builder SHALL NOT encode inferred `BondStereo` in the graph.

#### Scenario: Build beta-D-glucopyranose template
- **GIVEN** default beta-D-glucopyranose parameters and a bond length
- **WHEN** the Haworth geometry is built
- **THEN** ring and substituent atoms/bonds have stable topology, the declared front arc is bold consistently, and labels/bonds avoid obvious overlap in the production renderer.

### Requirement: Legacy special-template APIs remain compatible
Existing `build_haworth_template`, `build_linear_chain_template`, and `build_chair_template` entry points SHALL remain callable with their current signatures and return `MolGraph`, delegating to the projection builders where applicable.

#### Scenario: Existing template caller
- **GIVEN** a caller using a legacy template function
- **WHEN** it constructs the template as before
- **THEN** it receives a compatible deterministic `MolGraph` without opting into new UI or style behavior.

### Requirement: Fischer and chair have extensible parameterized bases
The projection package SHALL provide usable parameter dataclasses and deterministic base constructors for Fischer and chair geometry, including their canonical backbone/terminal or conformation/axial-equatorial policy surfaces, without requiring a complete GUI.

#### Scenario: Extend a non-Haworth projection
- **GIVEN** valid Fischer or chair parameters
- **WHEN** the corresponding base constructor runs
- **THEN** it returns stable geometry suitable for future renderer/UI integration and adds no implicit stereo metadata.
