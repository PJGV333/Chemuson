# architecture-boundaries Delta Specification

## MODIFIED Requirements

### Requirement: Temporary Exceptions Explicitly Documented

Every violation of a forbidden boundary SHALL be recorded exactly while it
exists. Eliminated violations SHALL be removed from the catalog and frozen
baseline; no replacement exception may be introduced for the same concern.

#### Scenario: M15 exceptions removed
- **GIVEN** autosave receives its collaborators through structural contracts
- **WHEN** the catalog is inspected
- **THEN** M15 has no temporary exceptions

## ADDED Requirements

### Requirement: Autosave Dependency Inversion

`chemuson.utils.autosave` SHALL not import ChemIO, GUI or PyQt6. It SHALL
receive document serialization and timer behavior through small typed
structural contracts or callbacks supplied by an authorized composition module.

#### Scenario: Autosave saves through injected collaborator
- **GIVEN** a fake document, undo stack, serializer and timer factory
- **WHEN** autosave writes a dirty document
- **THEN** it invokes the serializer with that document and writes its normal JSON payload

### Requirement: Autosave Behavioral Compatibility

Autosave SHALL preserve its observable directory, filename, metadata, backup
rotation, timing, error propagation and recovery-compatible payload behavior.

#### Scenario: Rotating backup behavior remains stable
- **GIVEN** an autosave manager with a backup limit of two
- **WHEN** it writes four dirty snapshots
- **THEN** exactly the two newest snapshots remain

### Requirement: M15 Debt Reduction

M15 SHALL have no current dependencies, target dependencies or temporary
exceptions, and no substitute exception may be added.

#### Scenario: Catalog observes a clean M15
- **GIVEN** the catalog and AST boundary analyzer
- **WHEN** they inspect M15
- **THEN** no M15-to-M01 runtime or M15-to-M09 TYPE_CHECKING crossing exists

### Requirement: Autosave Import Isolation

Importing `chemuson.utils.autosave` SHALL not load `chemuson.gui`,
`chemuson.chemio`, PyQt6 or RDKit transitively.

#### Scenario: Isolated subprocess import
- **GIVEN** a fresh Python subprocess
- **WHEN** it imports `chemuson.utils.autosave`
- **THEN** none of those modules are present in `sys.modules`
