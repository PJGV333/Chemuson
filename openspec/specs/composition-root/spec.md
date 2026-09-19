# composition-root Specification

## Purpose
TBD - created by archiving change audit-composition-root. Update Purpose after archive.
## Requirements
### Requirement: Composition Root Remains Consolidated

M19 SHALL remain the single owner of CLI parsing and application composition
under `src/chemuson/__main__.py` and `src/chemuson/app/`.

#### Scenario: Composition audit

- **WHEN** the architecture audit evaluates M19
- **THEN** it records `consolidated / no structural change required`

### Requirement: Reserved Slot Is Not a Placeholder

M24 SHALL remain reserved until a future boundary is justified by behavior and
ownership evidence.

#### Scenario: No forced composition split

- **WHEN** the current catalog is inspected after the audit
- **THEN** no M24 package or placeholder module exists

