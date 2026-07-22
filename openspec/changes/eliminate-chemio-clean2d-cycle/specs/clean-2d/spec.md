## ADDED Requirements

### Requirement: Clean2D Owns Imported Depiction Selection

Clean2D SHALL own imported-depiction quality, candidate construction, ranking,
scaffold and block-unwrap integration, best-depiction selection, and reports.

#### Scenario: SMILES depiction is requested
- **WHEN** a consumer imports a SMILES depiction through the deliberate Clean2D API
- **THEN** Clean2D ranks candidates and returns the same selected graph or controlled fallback

### Requirement: Imported Depiction Behavior Is Preserved

Imported depiction SHALL preserve candidate source names, deterministic
rejected-score-source ordering, score values, metadata, rejection reasons,
target bond length, worker timeouts, stereo metadata, and fallback errors.

#### Scenario: Candidate set is ranked
- **WHEN** worker candidates, scaffold candidates, or block-unwrap candidates are available
- **THEN** their accepted and rejected ordering and report metadata match the pre-migration behavior

### Requirement: Imported Depiction API Is Deliberate

Clean2D SHALL reexport `DepictionCandidate`,
`smiles_to_depiction_candidates`, `smiles_to_molgraph_best_depiction`, and
`smiles_to_molgraph_best_depiction_with_report` from `chemuson.clean2d`.

#### Scenario: Consumer uses the new path
- **WHEN** an external-to-M02 repository consumer needs imported depiction
- **THEN** it imports a declared Clean2D symbol rather than an M01 compatibility shim

### Requirement: Imported Candidates Do Not Mutate Their Base Graph

Scaffold and block-unwrap candidate construction SHALL operate on copies and
shall not mutate the base imported candidate graph.

#### Scenario: Derived candidate is evaluated
- **WHEN** a scaffold or block-unwrap coordinate set is scored
- **THEN** the original candidate graph retains its coordinates and metadata
