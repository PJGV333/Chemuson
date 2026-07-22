## ADDED Requirements

### Requirement: Imported Depiction Ownership Is Clean2D

The module catalog SHALL list `depiction_quality` and `imported_depiction` as
M02 internals and SHALL not list `depiction_candidates` as an M01 internal.

#### Scenario: Ownership inventory is inspected
- **WHEN** the catalog is validated
- **THEN** imported-depiction quality and orchestration belong only to M02

### Requirement: M01 and M02 Dependency State Is Cycle-Free

M01 current and target dependencies SHALL be exactly M00. M02 current and
target dependencies SHALL be exactly M00 and M01.

#### Scenario: Dependency state is checked
- **WHEN** the catalog is parsed after the migration
- **THEN** M01 has no M02 dependency and M02 retains its M01 dependency

### Requirement: Persistence Is the Only Temporary Exception

M01 SHALL retain exactly the M01-to-M09 TYPE_CHECKING exception in
`chemio/persistence.py`, while M02 SHALL have no temporary exceptions.

#### Scenario: Exception baseline is counted
- **WHEN** temporary exceptions are inventoried
- **THEN** there is one TYPE_CHECKING exception and zero runtime exceptions
