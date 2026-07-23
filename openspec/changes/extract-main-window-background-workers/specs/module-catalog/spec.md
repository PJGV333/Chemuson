## ADDED Requirements

### Requirement: M08 Catalogs GUI Background Workers

M08 SHALL inventory `gui/background_workers.py` as internal implementation and
its architectural characterization as an M08 test without adding a dependency,
temporary exception, cycle or public API.

#### Scenario: M08 inventory is inspected
- **GIVEN** the module catalog after extraction
- **WHEN** paths, tests and APIs are audited
- **THEN** the worker module and test are owned by M08 and dependency sets remain unchanged
