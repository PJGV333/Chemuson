## ADDED Requirements

### Requirement: Existing Domain Module Audit Is Explicit

The repository SHALL document an audit of M05, M06, M07, M14, M16, M17, M18,
M19 and M20 without renumbering or duplicating existing module ownership.

#### Scenario: Audit coverage is complete

- **WHEN** the audit report and catalog are inspected
- **THEN** every listed module has an explicit audited status and no module is
  structurally changed without a separate boundary decision

### Requirement: Audited Modules Have No Unrecorded Debt

Each audited module SHALL have equal current and target dependencies, no
`temporary_exceptions`, and no `circular_dependencies`.

#### Scenario: Catalog debt is checked

- **WHEN** the architecture audit test reads the catalog
- **THEN** all audited modules satisfy the zero-debt condition
