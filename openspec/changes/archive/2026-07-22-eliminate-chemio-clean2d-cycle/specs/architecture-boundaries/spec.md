## ADDED Requirements

### Requirement: ChemIO Does Not Depend on Clean2D

No M01 file SHALL import M02, including a local, lazy, aliased, reexported, or
TYPE_CHECKING import.

#### Scenario: Lazy import is analyzed
- **WHEN** the architecture analyzer inspects every M01 Python file
- **THEN** it finds no import targeting M02 regardless of the import's scope

### Requirement: Clean2D May Consume ChemIO Import Services

M02 SHALL be permitted to consume M01 format parsing and isolated worker
services while remaining independent of GUI modules.

#### Scenario: Imported depiction calls ChemIO
- **WHEN** Clean2D orchestrates SMILES depiction candidates
- **THEN** it may invoke M01 parsing and worker services without creating an M01-to-M02 edge

### Requirement: M01/M02 Cycle Removed

The catalog SHALL contain no M01/M02 circular dependency or any circular
dependency entry.

#### Scenario: Cycle inventory is empty
- **WHEN** the module catalog is inspected after the migration
- **THEN** its circular dependency inventory is empty

### Requirement: Runtime Exception Reduction

The five historical M01-to-M02 runtime exception identities SHALL not
reappear, and the only active exception SHALL be the M01-to-M09 persistence
TYPE_CHECKING crossing.

#### Scenario: Historical runtime edge is reintroduced
- **WHEN** any former M01-to-M02 exception identity is added to the catalog
- **THEN** the frozen exception test fails
