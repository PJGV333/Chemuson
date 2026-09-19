# platform-settings-boundary Specification

## Purpose
TBD - created by archiving change establish-platform-settings-boundary. Update Purpose after archive.
## Requirements
### Requirement: Platform Settings Is GUI-Independent

M21 SHALL own application settings and packaged-resource resolution without
importing `chemuson.gui`, `QtWidgets`, a main window, a canvas or a controller.

#### Scenario: Platform imports are isolated

- **WHEN** the platform architecture test inspects M21 source imports
- **THEN** it finds no GUI or widget dependency

### Requirement: Preference Behavior Is Preserved

The platform settings policy SHALL preserve the historical naming and numbering
keys, defaults and boolean coercion used by the GUI.

#### Scenario: Preference round-trip

- **WHEN** legacy values are loaded and typed preferences are saved
- **THEN** the normalized values and historical keys match the prior behavior

### Requirement: Resource Ownership Has One Canonical Implementation

M21 SHALL own the implementation of packaged-resource resolution, while
`chemuson.utils.resources` SHALL remain an import-only compatibility shim.

#### Scenario: Historical resource import

- **WHEN** an existing consumer imports `chemuson.utils.resources`
- **THEN** it resolves the canonical M21 helper without a duplicate definition

