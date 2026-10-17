# Projection stereochemistry contract

## ADDED Requirements

### Requirement: Projection semantics are explicit and distinct from 2D coordinates
The projection API SHALL keep renderer geometry separate from immutable typed semantic descriptors. Haworth descriptors SHALL preserve alpha/beta, D/L, numbered ring sites and each declared substituent orientation; Fischer descriptors SHALL state horizontal bonds point toward the observer and vertical bonds recede, preserving indexed left/right substituents without inferring undeclared R/S; chair descriptors SHALL preserve conformation and each substituent's up/down plus axial/equatorial relation. No projection SHALL encode these conventions using indiscriminate `BondStereo` or alkene E/Z annotations.

#### Scenario: Invert anomer and carbohydrate configuration
- **GIVEN** typed beta-D and alpha-D or beta-L Haworth geometries
- **WHEN** their semantic descriptors are compared
- **THEN** alpha/beta inversion changes the declared anomeric orientation, while D/L inversion mirrors the declared substituent pattern
- **AND** conversion to a graph invents no `BondStereo`.

### Requirement: In-memory projection instances preserve semantics without changing CMSN
A projection-aware instance SHALL carry both its graphical `MolGraph` and immutable semantic descriptor. Legacy `MolGraph` template adapters SHALL remain callable and are rendering-only. The implementation MUST NOT claim that plain `MolGraph` or current `.cmsn` serialization persists projection descriptors; no serialization schema change is permitted in this scope.

#### Scenario: Adapt projection for drawing
- **GIVEN** a Fischer, Haworth or chair projection instance
- **WHEN** its graphical graph is used by the existing renderer
- **THEN** the typed projection semantics remain available on the instance
- **AND** serializing the graph alone is documented/tested as a rendering-only, semantics-losing adaptation.

### Requirement: Chair flips preserve up/down while changing conformation relation
Chair projection descriptors SHALL distinguish conformation, axial/equatorial and up/down for each numbered carbon. Tests SHALL verify the geometric/conformational relation under chair flip without inferring absolute configuration not declared by input.

#### Scenario: Invert chair conformation
- **GIVEN** the same numbered substituent declared with up/down and axial/equatorial placement in a chair conformation
- **WHEN** the chair is flipped
- **THEN** the new geometry and descriptor state the resulting conformation explicitly
- **AND** the declared chemical up/down relation is not silently inverted or discarded.
