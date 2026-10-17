## ADDED Requirements

### Requirement: Local templates reach 3D generation without hidden graph mutation
ChemUSON SHALL calculate conformers from a defensive, diagnostically normalized graph copy, SHALL preserve the authored MolGraph, and SHALL never infer absent stereochemistry. Supported abbreviations SHALL be expanded using the existing declared mappings. A contradictory aromatic bridge SHALL only be represented as its stored single order when that order is explicitly single; ambiguous cases SHALL fail with a stage-specific diagnostic.

#### Scenario: Tetrandrina template reaches ETKDG
- **GIVEN** the exact user-local `Productos naturales / tetrandrina` template has 43 atoms, 48 bonds, supported abbreviations, and aromatic bridge bond 9 with explicit single order
- **WHEN** conformer generation begins
- **THEN** the saved template graph remains byte/field-equivalent and the calculation graph is expanded on a copy
- **AND** the invalid acyclic aromatic flag is reported and normalized only for RDKit calculation
- **AND** ETKDG coordinates are published for document atoms and scene-only expanded atoms
- **AND** no missing atom/bond stereochemistry is assigned

#### Scenario: Unambiguous aromatic rings remain aromatic
- **GIVEN** a graph whose aromatic edges participate in cycles and a graph with an aromatic bridge whose stored order is not explicit single
- **WHEN** calculation normalization runs
- **THEN** cycle aromaticity is preserved
- **AND** ambiguous bridge data returns a diagnostic rather than a guessed bond order

#### Scenario: Scene publication is distinct from computation
- **GIVEN** successful ETKDG coordinates but a nonconverged force-field minimization
- **WHEN** a scene is constructed
- **THEN** finite coordinates may be displayed with the nonconvergence diagnostic
- **AND** the UI does not label the geometry as minimized or visually accepted without real GUI verification
