## ADDED Requirements

### Requirement: M01 Has Only Core Dependencies

M01 current and target dependencies SHALL each be exactly `M00`; its temporary
exceptions and circular dependencies SHALL be empty. Persistence SHALL remain
an internal M01 API.

#### Scenario: M01 catalog entry is inspected
- **WHEN** the module catalog test reads M01
- **THEN** it finds only M00 dependencies and no GUI debt

### Requirement: Zero-Debt Catalog

Every module entry SHALL retain a `temporary_exceptions` field whose value is
an empty list, and all circular-dependency lists SHALL be empty.

#### Scenario: Full catalog is inspected
- **WHEN** all module entries are accumulated
- **THEN** no temporary exception or cycle is present
