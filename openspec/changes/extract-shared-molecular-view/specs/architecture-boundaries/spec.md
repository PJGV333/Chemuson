## ADDED Requirements

### Requirement: ChemCalc depends only on Core
M03 SHALL have M00 as its only current ChemUSON runtime dependency and SHALL
not import any M04 module at runtime.

#### Scenario: ChemCalc imports are analyzed
- **WHEN** the architecture test analyzes M03 source files
- **THEN** no M03-to-M04 import or M03 temporary exception is present
