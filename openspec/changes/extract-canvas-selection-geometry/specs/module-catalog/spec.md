## ADDED Requirements

### Requirement: M09 Catalogs Selection Geometry

M09 SHALL inventory `gui/canvas/selection_geometry.py` as internal
implementation and its architectural and numeric tests as M09 coverage without
adding a dependency, temporary exception, cycle or public API.

#### Scenario: M09 inventory is inspected
- **GIVEN** the module catalog after extraction
- **WHEN** paths, tests, APIs and dependencies are audited
- **THEN** selection geometry and its tests belong to M09 while dependency sets remain unchanged
