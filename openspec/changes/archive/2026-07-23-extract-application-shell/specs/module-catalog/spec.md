## ADDED Requirements

### Requirement: M08 Catalogs the Internal Application Shell

M08 SHALL own `src/chemuson/gui/shell/` as an internal implementation path and
SHALL retain `ChemusonWindow` as its public API. The extraction SHALL NOT add a
module dependency, temporary exception or cycle.

#### Scenario: M08 inventory is inspected
- **GIVEN** the module catalog after shell extraction
- **WHEN** paths, APIs and dependencies are audited
- **THEN** the shell has exclusive M08 ownership, `ChemusonWindow` remains public, and dependency sets are unchanged
