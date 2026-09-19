# resilience-boundary Specification

## Purpose
TBD - created by archiving change establish-resilience-boundary. Update Purpose after archive.
## Requirements
### Requirement: Resilience Has Canonical Ownership

M22 SHALL own the canonical autosave, recovery filesystem-policy and
crash-reporting implementations. M15 SHALL retain only import compatibility
shims for those historical paths.

#### Scenario: Canonical imports

- **WHEN** a runtime consumer imports autosave or crash reporting
- **THEN** the canonical import resolves under `chemuson.resilience`

### Requirement: Recovery Filesystem Policy Is GUI-Independent

M22 SHALL provide `read_autosave_metadata`, `list_autosave_entries` and
`archive_autosave` from `resilience/recovery.py` without importing
`chemuson.gui` or `PersistenceManager`. RecoveryController SHALL retain only
thin delegation plus Qt/document orchestration.

#### Scenario: Recovery ownership

- **WHEN** the recovery implementation and controller are inspected
- **THEN** filesystem policy has one M22 owner and the controller contains no
  duplicate policy implementation


#### Scenario: Import isolation

- **WHEN** the architecture contract inspects M22 source files
- **THEN** no ChemUSON GUI import is present

### Requirement: Runtime Behavior Is Preserved

The move SHALL preserve autosave metadata, rotation, recovery directory policy,
crash-log content and historical imports.

#### Scenario: Historical consumers

- **WHEN** a consumer imports `chemuson.utils.autosave` or
  `chemuson.utils.crash_reporter`
- **THEN** it receives the canonical M22 symbols without duplicate logic

