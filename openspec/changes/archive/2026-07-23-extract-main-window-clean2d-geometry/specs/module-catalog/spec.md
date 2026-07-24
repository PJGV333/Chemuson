## ADDED Requirements

### Requirement: M08 Catalogs Clean2D Geometry Integration

M08 SHALL inventory `gui/clean2d_geometry.py` as internal implementation and
its architectural characterization as an M08 test without adding a dependency,
temporary exception, cycle or public API.

#### Scenario: M08 inventory is inspected
- **GIVEN** the module catalog after extraction
- **WHEN** paths, tests and APIs are audited
- **THEN** the geometry module and test are owned by M08 and dependency sets remain unchanged
