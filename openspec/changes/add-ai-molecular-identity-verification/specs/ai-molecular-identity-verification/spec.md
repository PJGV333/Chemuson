# AI Molecular Identity Verification Specification

## ADDED Requirements

### Requirement: Structural validity and identity are separate
The system SHALL never treat ChemIO/RDKit parsing as proof that a proposal matches its requested name. Identity status SHALL be one of `not_applicable`, `unverified`, `verified`, `mismatch`, or `reference_error`.

#### Scenario: Valid parse but different named molecule
- **WHEN** ChemIO accepts a proposal but canonical identity differs from the reference for the requested molecule
- **THEN** structural validation remains successful and identity is `mismatch`.

### Requirement: Conservative name resolution
Only explicit simple name requests matching supported draw/generate forms SHALL be resolved. Open-ended structure-generation instructions SHALL be `not_applicable`. Verification SHALL use an injectable existing Name→Structure resolver and SHALL NOT ask another LLM to judge identity.

#### Scenario: Open-ended structure request
- **WHEN** a prompt asks for a molecule with abstract properties rather than naming one molecule
- **THEN** identity is `not_applicable` and no reference lookup occurs.

### Requirement: Strong isolated chemical identity comparison
The system SHALL compare canonical, stereochemistry-aware chemical identity using isolated ChemIO/RDKit-backed InChI conversion. It SHALL NOT compare literal SMILES strings or 2D drawings, and SHALL NOT import RDKit directly in application/UI code.

#### Scenario: Equivalent structures use different SMILES
- **WHEN** proposal and reference represent the same stereochemical identity with distinct SMILES text
- **THEN** identity is `verified`.

### Requirement: Mismatch is explicit and guarded
When a trusted reference differs from the proposal, the UI SHALL show structural validity and identity mismatch separately, SHALL NOT insert automatically, and SHALL require both an explicit override choice and a confirmation before insertion. A mismatch SHALL remain identifiable in the operation/report.

#### Scenario: User declines mismatch override
- **WHEN** a mismatch preview is shown and the user does not confirm Insert Anyway
- **THEN** no graph is inserted and editor state is unchanged.

### Requirement: Unavailable references do not imply correctness
Missing references SHALL produce `unverified`; resolver or canonicalization failures SHALL produce `reference_error`. Neither state SHALL be labeled verified.

#### Scenario: Reference resolver fails
- **WHEN** the resolver is unavailable or reports a network/reference error
- **THEN** identity is `reference_error` and ChemIO success is retained separately.

### Requirement: Offline semantic regression coverage
Tests SHALL use fake resolver references and no Internet. Coverage SHALL include equivalent-but-textually-different SMILES, valid mismatch, missing/failing resolver, open-ended prompt, independence from model response text, guarded mismatch insertion, explicit override, and `cholesterol-semantic-mismatch-01`.

#### Scenario: Cholesterol mismatch regression
- **WHEN** the fake reference resolver supplies the cholesterol identity and the fake model proposes another valid structure
- **THEN** `cholesterol-semantic-mismatch-01` reports `mismatch` without network access.
