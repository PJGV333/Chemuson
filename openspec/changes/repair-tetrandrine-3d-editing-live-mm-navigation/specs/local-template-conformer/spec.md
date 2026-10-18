## ADDED Requirements

### Requirement: Local templates reach 3D generation without hidden graph mutation
ChemUSON SHALL expand supported abbreviations only on a defensive calculation copy, preserve the authored MolGraph and specified stereochemistry, and never infer absent stereochemistry. ChemUSON SHALL reject any aromatic bond that is a graph bridge with a stage-specific diagnostic rather than silently removing its aromatic flag. Valid aromatic rings SHALL retain their aromaticity.

#### Scenario: Authorized tetrandrina entry matches reference connectivity
- **GIVEN** the locally saved `Productos naturales / tetrandrina` entry has the authorized repair: remove terminal duplicate atom 11 and replace bond 5–11 with the aromatic closure 5–6
- **WHEN** the MolBlock is imported and supported `CH3`/`OCH3` groups are expanded on a calculation copy
- **THEN** the expanded graph has formula `C38H42N2O6`, 46 heavy atoms, 52 heavy-atom bonds, 88 atoms with hydrogens, one component, and cycle rank 7
- **AND** it has 24 aromatic atoms and 24 aromatic bonds and is connectivity-isomorphic to PubChem CID 73078 ignoring unspecified stereochemistry
- **AND** all remaining authored 2D coordinates and specified attributes are preserved
- **AND** no stereochemistry is inferred

#### Scenario: Historical malformed fixture is rejected
- **GIVEN** the preserved pre-repair fixture contains aromatic bridge bond 9 from atom 5 to terminal atom 11 and lacks aromatic closure 5–6
- **WHEN** conformer generation begins
- **THEN** it fails before RDKit sanitization with a stage-specific invalid-aromatic-bridge diagnostic naming bond 9
- **AND** no aromatic flag or source connectivity is changed to hide the defect

#### Scenario: Corrected local graph reaches ETKDG and scene construction
- **GIVEN** the repaired local graph passes RDKit sanitization and ETKDG
- **WHEN** the worker returns finite coordinates
- **THEN** the scene contains 88 atoms and 94 bonds, including expanded abbreviation atoms and hydrogens
- **AND** MMFF94/UFF convergence, missing parameters, and unoptimized geometry remain distinctly reported

#### Scenario: Scene publication is distinct from computation
- **GIVEN** successful ETKDG coordinates but a nonconverged force-field minimization
- **WHEN** a scene is constructed
- **THEN** finite coordinates may be displayed with the nonconvergence diagnostic
- **AND** the UI does not label the geometry as minimized or visually accepted without real GUI verification
