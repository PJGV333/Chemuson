## MODIFIED Requirements

### Requirement: Module Identification with Stable IDs

Each module entry SHALL have a unique `id` field matching the pattern `M\d\d`. The catalog SHALL contain exactly 24 modules (M00-M23). IDs SHALL be persistent: once assigned, an ID is never reused even if the module is removed.

#### Scenario: ID uniqueness
- **GIVEN** the catalog contains 24 module entries
- **WHEN** a test extracts all `id` values
- **THEN** all 24 IDs are distinct

#### Scenario: ID format
- **GIVEN** a module entry with id `M05`
- **WHEN** the id is validated against pattern `^M\d\d$`
- **THEN** the format is valid

## ADDED Requirements

### Requirement: M23 Catalogs the Molecular Assistant Boundary

M23 SHALL own `src/chemuson/molecular_assistant/`, SHALL list M00 and M01 as its current and target dependencies, and SHALL forbid dependencies on M02, M04, M08-M13, M16 and M19. M00, M01 and M02 SHALL forbid a dependency on M23. M23 SHALL have no temporary exceptions or circular dependencies.

#### Scenario: M23 catalog ownership is inspected
- **WHEN** the module catalog and source tree are audited
- **THEN** M23 exclusively owns the molecular-assistant package and tests, depends only on M00/M01, and records the required forbidden dependencies

#### Scenario: Clean2D remains AI-independent
- **WHEN** M02 dependencies and forbidden dependencies are inspected
- **THEN** M02 has no M23 current or target dependency and explicitly forbids importing M23
