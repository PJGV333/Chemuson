## ADDED Requirements

### Requirement: Clean 2D improves simple aromatic layout selection
Clean 2D SHALL be able to select a safe, measurably better layout candidate for simple isolated or substituted aromatic molecules without weakening chemical identity or complex-structure protections.

#### Scenario: Simple aromatic candidate improves measurably
- **GIVEN** a simple isolated or substituted aromatic molecule with an available safe candidate that improves at least one tracked geometry or visual metric without worsening crossing or ring-degeneracy safety
- **WHEN** Clean 2D selects and applies a candidate
- **THEN** at least one tracked simple-aromatic corpus case SHALL show a measurable improvement such as higher visual score, lower bond-length error, lower angle deviation, or safer nonbonded spacing

#### Scenario: Chemical identity remains preserved
- **GIVEN** Clean 2D improves a simple aromatic layout
- **WHEN** the result is compared with the input molecule
- **THEN** atom IDs, bond IDs, atom count, bond count, elements, charges, bond endpoints, bond orders, aromaticity, stereo metadata, and selection metadata SHALL remain unchanged

#### Scenario: Protected complex behavior remains unchanged
- **GIVEN** a high-risk complex-policy molecule or protected tetrandrine-like molecule
- **WHEN** Clean 2D runs after the simple aromatic improvement
- **THEN** the result state and protection behavior SHALL NOT change because of the simple aromatic improvement

#### Scenario: Baseline report changes remain reviewable
- **GIVEN** before and after Clean 2D baseline reports for the simple aromatic improvement
- **WHEN** the reports are compared and reviewed
- **THEN** changed cases SHALL be attributable to expected simple-aromatic behavior or observational geometry metrics
- **AND** the change SHALL NOT introduce unexpected contract-risk outside the intended family
