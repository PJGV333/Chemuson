## ADDED Requirements

### Requirement: Clean 2D improves monosubstituted aromatic layout selection
Clean 2D SHALL be able to select a safe, measurably better layout candidate for simple monosubstituted aromatic molecules without changing propose behavior or protected complex behavior.

#### Scenario: Monosubstituted aromatic candidate improves measurably
- **GIVEN** a simple monosubstituted aromatic molecule with one neutral 5/6-member aromatic ring and one small neutral terminal substituent
- **WHEN** Clean 2D runs in quick or publication mode
- **THEN** Clean 2D SHALL be able to select a non-current monosubstituted aromatic candidate that improves at least one tracked geometry or visual metric without increasing crossings or ring-degeneracy risk

#### Scenario: Substituent is placed outward
- **GIVEN** Clean 2D selects a monosubstituted aromatic template candidate
- **WHEN** the candidate coordinates are inspected
- **THEN** the external substituent SHALL be positioned outward from the aromatic ring anchor in a reasonable radial direction

#### Scenario: Chemical identity remains preserved
- **GIVEN** Clean 2D improves a monosubstituted aromatic layout
- **WHEN** the result is compared with the input molecule
- **THEN** atom IDs, bond IDs, atom count, bond count, elements, charges, bond endpoints, bond orders, aromaticity, stereo metadata, and selection metadata SHALL remain unchanged

#### Scenario: Propose and protected structures are excluded
- **GIVEN** propose mode, a fused aromatic, a high-risk complex-policy molecule, or a protected tetrandrine-like molecule
- **WHEN** Clean 2D candidates are generated or selected
- **THEN** the monosubstituted aromatic template candidate SHALL NOT be used
