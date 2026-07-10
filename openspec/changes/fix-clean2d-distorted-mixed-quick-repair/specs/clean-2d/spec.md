## ADDED Requirements

### Requirement: Quick repair rejects individual structural-bond extremes
When a Quick or Publication Clean2D run repairs a baseline classified as `needs_rebuild`, the engine SHALL reject a candidate if any structural bond is shorter than 70% or longer than 145% of its order-specific desired length.

#### Scenario: One short bond hidden by a correct average
- **WHEN** a repair candidate has an acceptable average bond length but a structural bond below 70% of its desired length
- **THEN** the engine rejects that candidate and records the short bond ID

#### Scenario: Safe fallback is available
- **WHEN** a candidate is rejected for an individual bond extreme and another candidate passes all safety gates
- **THEN** ranking selects the safe candidate without forcing its source

### Requirement: Block layouts follow the common repair gate
`block_layout` candidates SHALL pass chemical invariants, stereo preservation, crossing, collision, ring-degeneracy, final-quality, and individual-bond-range checks before acceptance.

#### Scenario: Block layout remains needs_rebuild
- **WHEN** a `block_layout` candidate on a `needs_rebuild` baseline still has final quality `needs_rebuild`
- **THEN** the engine rejects it and considers the remaining candidates

### Requirement: Propose policy remains unchanged
Direct `PROPOSE` SHALL reject any baseline that is not `good`.

#### Scenario: Needs-polish base
- **WHEN** direct `PROPOSE` receives a `needs_polish` baseline
- **THEN** it returns the existing optimization-required rejection without selecting a repair candidate
