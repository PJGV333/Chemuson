# Spec: Module Catalog

## ADDED Requirements

### Requirement: Module Catalog File Exists and is Valid YAML

The module catalog SHALL exist at `architecture/modules.yml` and SHALL be parseable as valid YAML containing a `modules` list.

#### Scenario: Catalog is readable
- **GIVEN** the file `architecture/modules.yml` exists
- **WHEN** a test parses it as YAML
- **THEN** the top-level key `modules` contains a list of module entries

### Requirement: Module Identification with Stable IDs

Each module entry SHALL have a unique `id` field matching the pattern `M\d\d`. The catalog SHALL contain exactly 19 modules (M00-M18). IDs SHALL be persistent: once assigned, an ID is never reused even if a module is removed.

#### Scenario: ID uniqueness
- **GIVEN** the catalog contains 19 module entries
- **WHEN** a test extracts all `id` values
- **THEN** all 19 IDs are distinct

#### Scenario: ID format
- **GIVEN** a module entry with id `M05`
- **WHEN** the id is validated against pattern `^M\d\d$`
- **THEN** the format is valid

### Requirement: Mandatory Fields Present

Each module entry SHALL contain: `id`, `name`, `title`, `responsibility`, `paths`, `status`, `risk_level`, `current_dependencies`, `target_dependencies`. The `paths` list SHALL reference directories or files that exist on disk relative to the repository root.

#### Scenario: Path existence
- **GIVEN** a module entry with `paths: ["src/chemuson/core/"]`
- **WHEN** a test checks each path against the filesystem
- **THEN** all paths resolve to existing entries

#### Scenario: Status validity
- **GIVEN** a module entry with `status: stable`
- **WHEN** the status is checked against allowed values
- **THEN** it is one of: `stable`, `evolving`, `legacy`

### Requirement: Dependency Fields Distinct and Consistent

- `current_dependencies` SHALL list module IDs actually imported at runtime (verified by AST analysis).
- `target_dependencies` SHALL list module IDs permitted in the target architecture.
- `forbidden_dependencies` SHALL list module IDs that must never be imported.
- A module ID SHALL NOT appear in both `target_dependencies` and `forbidden_dependencies` of the same module.
- Dependencies guarded by `TYPE_CHECKING` SHALL NOT appear in `current_dependencies`; they SHALL be recorded in `temporary_exceptions` with `type_checking_only: true`.

#### Scenario: No contradictory dependency rules
- **GIVEN** module M02 has `target_dependencies: [M00]` and `forbidden_dependencies: [M08]`
- **WHEN** a test checks for overlap between the two lists
- **THEN** the intersection is empty

### Requirement: Public API Accuracy

When `public_api` is populated, each symbol SHALL exist in the module's `__init__.py` (declarable via `__all__` or direct assignment). Validation SHALL use AST analysis, not runtime import.

#### Scenario: Public API symbol exists
- **GIVEN** module M00 declares `public_api` containing `MolGraph`
- **WHEN** a test parses `src/chemuson/core/__init__.py` with AST
- **THEN** `MolGraph` is found as a name definition or import target

### Requirement: Circular Dependencies Documented

Known circular dependencies SHALL be registered in `circular_dependencies` with `modules`, `edges`, `severity`, and `resolution_plan`. The current codebase has two high/medium severity cycles (chemio↔clean2d, chemcalc↔chemname) that MUST be documented.

#### Scenario: Circular dependency registered
- **GIVEN** the catalog's `circular_dependencies` section
- **WHEN** a test checks for the chemio↔clean2d cycle
- **THEN** an entry exists with both module IDs and resolution plan

### Requirement: Empty Directories Not Cataloged

Directories that exist but contain no functional code (e.g., `templates/` with 0 entries) SHALL NOT receive a module ID.

#### Scenario: Empty directory excluded
- **GIVEN** `src/chemuson/templates/` exists with no Python files
- **WHEN** the catalog is inspected for module entries
- **THEN** no entry maps to `templates`
