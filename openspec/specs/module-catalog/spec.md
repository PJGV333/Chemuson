# module-catalog Specification

## Purpose
Define `architecture/modules.yml` as the authoritative inventory of ChemUSON
module ownership, responsibilities, APIs, dependencies, risk, tests and
temporary exceptions.
## Requirements
### Requirement: Module Catalog File Exists and is Valid YAML

The module catalog SHALL exist at `architecture/modules.yml` and SHALL be parseable as valid YAML containing a `modules` list.

#### Scenario: Catalog is readable
- **GIVEN** the file `architecture/modules.yml` exists
- **WHEN** a test parses it as YAML
- **THEN** the top-level key `modules` contains a list of module entries

### Requirement: Module Identification with Stable IDs

Each module entry SHALL have a unique `id` field matching the pattern `M\d\d`. The catalog SHALL contain exactly 20 modules (M00-M19). IDs SHALL be persistent: once assigned, an ID is never reused even if a module is removed.

#### Scenario: ID uniqueness
- **GIVEN** the catalog contains 20 module entries
- **WHEN** a test extracts all `id` values
- **THEN** all 20 IDs are distinct

#### Scenario: ID format
- **GIVEN** a module entry with id `M05`
- **WHEN** the id is validated against pattern `^M\d\d$`
- **THEN** the format is valid

### Requirement: Mandatory Fields Present

Each module entry SHALL contain: `id`, `name`, `title`, `responsibility`, `paths`, `status`, `risk_level`, `current_dependencies`, `target_dependencies`. The `paths` list SHALL reference directories or files that exist on disk relative to the repository root. Each Python file in `src/chemuson/` SHALL belong to at most one module. Declared paths SHALL NOT overlap. No module SHALL depend on itself.

#### Scenario: Path existence
- **GIVEN** a module entry with `paths: ["src/chemuson/core/"]`
- **WHEN** a test checks each path against the filesystem
- **THEN** all paths resolve to existing entries

#### Scenario: Status validity
- **GIVEN** a module entry with `status: stable`
- **WHEN** the status is checked against allowed values
- **THEN** it is one of: `stable`, `evolving`, `legacy`

### Requirement: Dependency Fields Distinct and Consistent

- `current_dependencies` SHALL list real runtime imports; imports within functions or lazy are still runtime; imports under `TYPE_CHECKING` do NOT enter `current_dependencies`.
- `target_dependencies` SHALL list dependencies permitted by the target architecture.
- `forbidden_dependencies` SHALL list dependencies that must never be imported.
- A module ID SHALL NOT appear in both `target_dependencies` and `forbidden_dependencies` of the same module.
- A runtime exception represents a current runtime dependency that is either absent from `target_dependencies` or present in `forbidden_dependencies`. A TYPE_CHECKING-only crossing may be documented with `type_checking_only: true`, but it must not appear in `current_dependencies`.

#### Scenario: No contradictory dependency rules
- **GIVEN** module M02 has `target_dependencies: [M00]` and `forbidden_dependencies: [M08]`
- **WHEN** a test checks for overlap between the two lists
- **THEN** the intersection is empty

### Requirement: Public API Accuracy

`public_api` SHALL contain deliberately public contracts: symbols listed in `__all__`; deliberate re-exports from `__init__.py`; functions associated with `[project.scripts]`; for a single-file module, a deliberately cataloged function like `M19.main`. A symbol imported directly from an internal module by another package does NOT automatically become public API. Validation SHALL use AST analysis, not runtime import.

#### Scenario: Public API symbol exists
- **GIVEN** module M00 declares `public_api` containing `MolGraph`
- **WHEN** a test parses `src/chemuson/core/__init__.py` with AST
- **THEN** `MolGraph` is found as a name definition or import target

### Requirement: Circular Dependencies Documented

Known circular dependencies SHALL be registered in `circular_dependencies` with
`modules`, `edges`, `severity`, and `resolution_plan`. The current codebase has
one high severity chemio-to-clean2d cycle. The resolved ChemCalc-to-ChemName
cycle SHALL NOT be registered.

#### Scenario: Resolved ChemCalc-to-ChemName cycle is absent
- **WHEN** the catalog is inspected after ChemCalc uses the Core molecular view
- **THEN** no circular dependency entry contains both M03 and M04

### Requirement: Empty Directories Not Cataloged

Directories that exist but contain no functional code (e.g., `templates/` with 0 entries) SHALL NOT receive a module ID.

#### Scenario: Empty directory excluded
- **GIVEN** `src/chemuson/templates/` exists with no Python files
- **WHEN** the catalog is inspected for module entries
- **THEN** no entry maps to `templates`

### Requirement: Imported Depiction Ownership Is Clean2D

The module catalog SHALL list `depiction_quality` and `imported_depiction` as
M02 internals and SHALL not list `depiction_candidates` as an M01 internal.

#### Scenario: Ownership inventory is inspected
- **WHEN** the catalog is validated
- **THEN** imported-depiction quality and orchestration belong only to M02

### Requirement: M01 and M02 Dependency State Is Cycle-Free

M01 current and target dependencies SHALL be exactly M00. M02 current and
target dependencies SHALL be exactly M00 and M01.

#### Scenario: Dependency state is checked
- **WHEN** the catalog is parsed after the migration
- **THEN** M01 has no M02 dependency and M02 retains its M01 dependency

### Requirement: Persistence Is the Only Temporary Exception

M01 SHALL retain exactly the M01-to-M09 TYPE_CHECKING exception in
`chemio/persistence.py`, while M02 SHALL have no temporary exceptions.

#### Scenario: Exception baseline is counted
- **WHEN** temporary exceptions are inventoried
- **THEN** there is one TYPE_CHECKING exception and zero runtime exceptions

### Requirement: M01 Has Only Core Dependencies

M01 current and target dependencies SHALL each be exactly `M00`; its temporary
exceptions and circular dependencies SHALL be empty. Persistence SHALL remain
an internal M01 API.

#### Scenario: M01 catalog entry is inspected
- **WHEN** the module catalog test reads M01
- **THEN** it finds only M00 dependencies and no GUI debt

### Requirement: Zero-Debt Catalog

Every module entry SHALL retain a `temporary_exceptions` field whose value is
an empty list, and all circular-dependency lists SHALL be empty.

#### Scenario: Full catalog is inspected
- **WHEN** all module entries are accumulated
- **THEN** no temporary exception or cycle is present

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

### Requirement: M08 Catalogs the Internal Application Shell

M08 SHALL own `src/chemuson/gui/shell/` as an internal implementation path and
SHALL retain `ChemusonWindow` as its public API. The extraction SHALL NOT add a
module dependency, temporary exception or cycle.

#### Scenario: M08 inventory is inspected
- **GIVEN** the module catalog after shell extraction
- **WHEN** paths, APIs and dependencies are audited
- **THEN** the shell has exclusive M08 ownership, `ChemusonWindow` remains public, and dependency sets are unchanged

### Requirement: M08 Catalogs GUI Background Workers

M08 SHALL inventory `gui/background_workers.py` as internal implementation and
its architectural characterization as an M08 test without adding a dependency,
temporary exception, cycle or public API.

#### Scenario: M08 inventory is inspected
- **GIVEN** the module catalog after extraction
- **WHEN** paths, tests and APIs are audited
- **THEN** the worker module and test are owned by M08 and dependency sets remain unchanged
