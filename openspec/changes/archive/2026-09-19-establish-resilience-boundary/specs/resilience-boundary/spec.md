## ADDED Requirements

### Requirement: Resilience Has Canonical Ownership

M22 SHALL own the canonical autosave and crash-reporting implementations. M15
SHALL retain only import compatibility shims for those historical paths.

#### Scenario: Canonical imports

- **WHEN** a runtime consumer imports autosave or crash reporting
- **THEN** the canonical import resolves under `chemuson.resilience`

### Requirement: Resilience Is GUI-Independent

M22 SHALL NOT import `chemuson.gui` or depend on ChemUSON widgets, controllers,
bootstrap or persistence modules.

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
