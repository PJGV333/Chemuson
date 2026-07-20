# Spec: Architecture Boundaries

## ADDED Requirements

### Requirement: Path Ownership and Exclusivity

Each Python file in `src/chemuson/` SHALL belong to at most one module. Declared paths SHALL NOT overlap. All declared paths MUST exist. No module SHALL depend on itself. These rules are verified by the module catalog tests.

#### Scenario: Path exclusivity verified
- **GIVEN** two module entries with paths `["src/chemuson/core/"]` and `["src/chemuson/chemio/"]`
- **WHEN** the catalog is validated for path overlap
- **THEN** no file belongs to more than one module

#### Scenario: Self-dependency prohibited
- **GIVEN** a module entry with id `M00` and paths `["src/chemuson/core/"]`
- **WHEN** dependencies are checked against paths
- **THEN** M00 does not depend on itself

### Requirement: Forbidden Imports Detected by Static Analysis

The boundary test suite SHALL analyze all `.py` files in `src/chemuson/` using Python's `ast` module. Any import statement that crosses a `forbidden_dependencies` boundary SHALL produce a test failure with file path, line number, import expression, and violated rule.

#### Scenario: Forbidden import detected
- **GIVEN** module M00 (core) has `forbidden_dependencies: [M08]`
- **WHEN** `src/chemuson/core/model.py` contains `from chemuson.gui import Something`
- **THEN** the test fails with a message identifying the file, line, import, and rule

### Requirement: TYPE_CHECKING Imports Excluded from Runtime Dependencies

Imports wrapped in `if TYPE_CHECKING:` blocks SHALL NOT be treated as runtime dependencies. A runtime exception represents an actual not permitted or forbidden dependency; a crossing under `TYPE_CHECKING` MAY be recorded as exception with `type_checking_only: true`. The AST analysis SHALL detect the `TYPE_CHECKING` guard and skip those imports when checking boundary rules.

#### Scenario: TYPE_CHECKING import allowed
- **GIVEN** `src/chemuson/chemio/persistence.py` imports `chemuson.gui.canvas` under `TYPE_CHECKING`
- **WHEN** the boundary test analyzes the file
- **THEN** the import does not trigger a forbidden dependency violation

### Requirement: Temporary Exceptions Explicitly Documented

Every violation of a `forbidden_dependencies` rule that currently exists in the codebase SHALL be recorded as a `temporary_exception` in `architecture/modules.yml` with all mandatory fields: `source_id`, `target_id`, `file`, `import_path`, `reason`, `debt_ref`, `elimination_condition`, `type_checking_only`.

#### Scenario: Exception prevents test failure
- **GIVEN** an exception entry for `utils/autosave.py` importing `gui.canvas` with `type_checking_only: true`
- **WHEN** the boundary test encounters that import
- **THEN** the test allows it and does not fail

### Requirement: Exception List Cannot Grow Silently

The test suite SHALL verify that the set of active violations matches the documented `temporary_exceptions` exactly. A new violation not covered by an exception entry SHALL fail the test, preventing silent growth of technical debt.

#### Scenario: New undocumented violation fails
- **GIVEN** the catalog has 3 temporary exceptions
- **WHEN** a developer adds a 4th forbidden import without updating the catalog
- **THEN** the boundary test fails identifying the new violation

### Requirement: No Wildcards in Exceptions

Exception entries SHALL specify exact file paths and exact import expressions. Patterns like `gui/*`, `utils/*`, or `any import from gui` SHALL NOT be accepted as valid exceptions.

#### Scenario: Wildcard exception rejected
- **GIVEN** an exception entry with `file: "gui/*"`
- **WHEN** the test validates exception entries
- **THEN** the entry is rejected as invalid

### Requirement: Tools Isolation Enforced

No file in `src/chemuson/` SHALL import from `tools/` or `chemuson.tools`. This rule has no exceptions.

#### Scenario: Tools import blocked
- **GIVEN** a file in `src/chemuson/core/` imports `tools.chemname_acceptance`
- **WHEN** the isolation test runs
- **THEN** the test fails with the specific file and import identified

### Requirement: Exception Entry Schema Validated

Every `temporary_exception` entry SHALL be validated for completeness: all 8 mandatory fields present, `source_id` and `target_id` are valid module IDs, `file` is an existing path, `type_checking_only` is boolean.

#### Scenario: Incomplete exception rejected
- **GIVEN** an exception entry missing `elimination_condition`
- **WHEN** the catalog validation test runs
- **THEN** the test fails identifying the incomplete entry
