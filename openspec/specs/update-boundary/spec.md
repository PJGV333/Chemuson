# update-boundary Specification

## Purpose
TBD - created by archiving change audit-update-boundary. Update Purpose after archive.
## Requirements
### Requirement: Update Boundary Remains Cohesive

M14 SHALL remain the canonical owner of `src/chemuson/update/`, including
policy, provider, semantic-version, security, portable, Windows, rollback and
telemetry responsibilities.

#### Scenario: Update audit

- **WHEN** the architecture audit evaluates M14
- **THEN** it records `audited / no structural change required`

### Requirement: Reserved Slot Is Not a Placeholder

M24 SHALL remain reserved until a future change demonstrates a cohesive boundary. M23 is assigned by this change to the independent `molecular_assistant` boundary and SHALL NOT be treated as part of M14 `update`.

#### Scenario: No forced extraction

- **WHEN** the current catalog is inspected after the audit
- **THEN** M14 retains its existing ownership, M23 owns only molecular assistance, and no M24 package or placeholder module exists

