## ADDED Requirements

### Requirement: ChemIO Is GUI-Independent

M01 SHALL not import M09 at runtime, under `TYPE_CHECKING`, locally, lazily,
through aliases, or through reexports.

#### Scenario: GUI import is inspected in every scope
- **WHEN** the architecture analyzer examines every M01 Python file
- **THEN** it finds no M01-to-M09 import

### Requirement: Persistence Uses Structural Contracts

`PersistenceManager` SHALL accept a GUI-neutral structural contract and SHALL
not name a GUI class or invoke a GUI-private method.

#### Scenario: Pure document is persisted
- **WHEN** a non-Qt fake provides the persistence document operations
- **THEN** save and load complete without importing GUI

### Requirement: No Temporary Exceptions

The module catalog SHALL contain zero `temporary_exceptions` in every module.

#### Scenario: Catalog is inventoried
- **WHEN** architecture tests inspect all module entries
- **THEN** the accumulated temporary-exception list is empty

### Requirement: Eliminated Persistence Exception Cannot Reappear

The former M01-to-M09 `ChemusonCanvas` TYPE_CHECKING identity SHALL be rejected
if it reappears, including normalized file-path or whitespace variants.

#### Scenario: Historical identity is reintroduced
- **WHEN** the frozen exception audit receives the historical identity
- **THEN** it reports the identity as unexpected and M01 as growing from zero
