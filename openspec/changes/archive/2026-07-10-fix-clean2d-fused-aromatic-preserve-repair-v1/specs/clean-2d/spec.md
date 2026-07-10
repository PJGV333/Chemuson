## ADDED Requirements

### Requirement: Simple fused aromatic preserve repair
Clean2D SHALL provide a safe deterministic repair candidate for an unsubstituted neutral naphthalene-like system of exactly two aromatic six-membered rings sharing one edge.

#### Scenario: Collapsed naphthalene
- **WHEN** a whole simple fused aromatic system is classified for preserve-only handling with collapsed geometry
- **THEN** Clean2D SHALL select a safe fused-aromatic template with non-degenerate rings and valid individual bond lengths.

#### Scenario: Unsupported fused system
- **WHEN** a fused system has substituents, stereo, non-aromatic bonds, more rings, or ambiguous topology
- **THEN** Clean2D SHALL NOT generate the fused-aromatic template.
