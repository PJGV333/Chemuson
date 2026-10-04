## MODIFIED Requirements

### Requirement: Reserved Slot Is Not a Placeholder

M24 SHALL remain reserved until a future change demonstrates a cohesive boundary. M23 is assigned by this change to the independent `molecular_assistant` boundary and SHALL NOT be treated as part of M14 `update`.

#### Scenario: No forced extraction

- **WHEN** the current catalog is inspected after the audit
- **THEN** M14 retains its existing ownership, M23 owns only molecular assistance, and no M24 package or placeholder module exists
