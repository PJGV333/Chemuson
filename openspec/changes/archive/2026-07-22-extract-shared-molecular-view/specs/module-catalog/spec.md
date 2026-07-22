## MODIFIED Requirements

### Requirement: Circular Dependencies Documented

Known circular dependencies SHALL be registered in `circular_dependencies` with
`modules`, `edges`, `severity`, and `resolution_plan`. The current codebase has
one high severity chemio-to-clean2d cycle. The resolved ChemCalc-to-ChemName
cycle SHALL NOT be registered.

#### Scenario: Resolved ChemCalc-to-ChemName cycle is absent
- **WHEN** the catalog is inspected after ChemCalc uses the Core molecular view
- **THEN** no circular dependency entry contains both M03 and M04
