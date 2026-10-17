# Haworth editorial geometry

## ADDED Requirements

### Requirement: Canonical Haworth ring and vertex-aware substituents
The Haworth builder SHALL define separate canonical, non-regular piranose and furanose geometry, stable ring numbering, validated aspect/perspective/tilt parameters, and an explicit canonical front-edge policy. Substituent vectors SHALL depend on ring position and declared up/down orientation rather than a generic radial rule. The beta-D-glucopyranose reference SHALL depict C1 OH up, C2 OH down, C3 OH up, C4 OH down and C5 CH2OH up.

#### Scenario: C3 hydroxyl points out of beta-D-glucose ring
- **GIVEN** default beta-D-glucopyranose parameters
- **WHEN** renderer geometry is created at supported scales
- **THEN** C3 OH is directed up/out, its label is outside the ring, and no label/bond/front-edge collision is introduced.

### Requirement: Haworth visual refinement does not alter chemistry or AB3
Haworth refinement SHALL change only projection geometry and local front-edge emphasis; it SHALL preserve declared topology, anomer/D-L orientation and the exact opt-in AB3 parameter vector. Visual comparison SHALL use equal scale and the real renderer, without declaring owner aesthetic acceptance.

#### Scenario: Compare before and after at equal scale
- **GIVEN** baseline and corrected C3 exports using the same style and viewport dimensions
- **WHEN** they are compared
- **THEN** the evidence records both source/provenance and matching scale
- **AND** the change does not claim final owner visual approval.
