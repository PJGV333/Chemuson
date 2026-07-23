## MODIFIED Requirements

### Requirement: M19 Catalogs the Composition Root

M19 SHALL own `src/chemuson/__main__.py` and `src/chemuson/app/`. Its current
and target dependencies SHALL be M08, M15 and M18, matching its direct imports.
Its public API SHALL continue to contain `main`, associated with the installed
entry point.

#### Scenario: M19 inventory is inspected
- **GIVEN** the module catalog after establishing the composition root
- **WHEN** M19 paths, dependencies and public API are audited
- **THEN** both bootstrap paths exist, dependencies are exactly M08/M15/M18, and `main` resolves statically

#### Scenario: Bootstrap paths have exclusive ownership
- **GIVEN** all catalog paths
- **WHEN** `src/chemuson/__main__.py` and files below `src/chemuson/app/` are resolved
- **THEN** they belong to M19 and no other module
