# architecture-boundaries Specification

## Purpose
Define enforceable module ownership, dependency, exception, public-surface and
tool-isolation rules for Python code under `src/chemuson`.
## Requirements
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

Imports wrapped in `if TYPE_CHECKING:` blocks SHALL NOT be treated as runtime dependencies. A runtime exception represents a current runtime dependency that is either absent from `target_dependencies` or present in `forbidden_dependencies`. A TYPE_CHECKING-only crossing may be documented with `type_checking_only: true`, but it must not appear in `current_dependencies`. The AST analysis SHALL detect the `TYPE_CHECKING` guard and skip those imports when checking boundary rules.

#### Scenario: TYPE_CHECKING import allowed
- **GIVEN** `src/chemuson/chemio/persistence.py` imports `chemuson.gui.canvas` under `TYPE_CHECKING`
- **WHEN** the boundary test analyzes the file
- **THEN** the import does not trigger a forbidden dependency violation

### Requirement: Temporary Exceptions Explicitly Documented

Every violation of a forbidden boundary SHALL be recorded exactly while it
exists. Eliminated violations SHALL be removed from the catalog and frozen
baseline; no replacement exception may be introduced for the same concern.

#### Scenario: M15 exceptions removed
- **GIVEN** autosave receives its collaborators through structural contracts
- **WHEN** the catalog is inspected
- **THEN** M15 has no temporary exceptions

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

### Requirement: Autosave Dependency Inversion

`chemuson.utils.autosave` SHALL not import ChemIO, GUI or PyQt6. It SHALL
receive document serialization and timer behavior through small typed
structural contracts or callbacks supplied by an authorized composition module.

#### Scenario: Autosave saves through injected collaborator
- **GIVEN** a fake document, undo stack, serializer and timer factory
- **WHEN** autosave writes a dirty document
- **THEN** it invokes the serializer with that document and writes its normal JSON payload

### Requirement: Autosave Behavioral Compatibility

Autosave SHALL preserve its observable directory, filename, metadata, backup
rotation, timing, error propagation and recovery-compatible payload behavior.

#### Scenario: Rotating backup behavior remains stable
- **GIVEN** an autosave manager with a backup limit of two
- **WHEN** it writes four dirty snapshots
- **THEN** exactly the two newest snapshots remain

### Requirement: M15 Debt Reduction

M15 SHALL have no current dependencies, target dependencies or temporary
exceptions, and no substitute exception may be added.

#### Scenario: Catalog observes a clean M15
- **GIVEN** the catalog and AST boundary analyzer
- **WHEN** they inspect M15
- **THEN** no M15-to-M01 runtime or M15-to-M09 TYPE_CHECKING crossing exists

### Requirement: Autosave Import Isolation

Importing `chemuson.utils.autosave` SHALL not load `chemuson.gui`,
`chemuson.chemio`, PyQt6 or RDKit transitively.

#### Scenario: Isolated subprocess import
- **GIVEN** a fresh Python subprocess
- **WHEN** it imports `chemuson.utils.autosave`
- **THEN** none of those modules are present in `sys.modules`

#### Scenario: RDKit remains absent after isolated import
- **GIVEN** a fresh Python subprocess
- **WHEN** it imports `chemuson.utils.autosave`
- **THEN** no module named `rdkit` or prefixed with `rdkit.` is present in `sys.modules`

### Requirement: Type-Safe Autosave Composition

`AutosaveManager` SHALL be generic over its document type and its serializer
SHALL accept that same type. The Qt adapter SHALL implement the controller
contract explicitly, without `__getattr__`, and SHALL not claim to be the core
manager. The autosave factory SHALL return `AutosaveController`; it SHALL not
use `Callable[..., ...]`, and the adapter SHALL not use `object` as a stand-in
for its concrete document type.

#### Scenario: Canvas manager creates a typed controller
- **GIVEN** `CanvasTabManager` and its autosave factory
- **WHEN** it creates a document tab
- **THEN** it stores an `AutosaveController` and supplies a serializer accepting the same canvas type

### Requirement: Autosave Lifecycle Compatibility

Autosave SHALL preserve idempotent `start`, periodic timer startup, dirty-only
initial debounce, stopping of both timers, restart and cancellation of debounce,
forced cleanup snapshots, and the dirty/clean policy executed by timer
callbacks.

#### Scenario: Lifecycle transitions preserve timer policy
- **GIVEN** an autosave manager with observable timers and clean state
- **WHEN** it starts, stops, restarts debounce, cancels debounce and cleans up after save
- **THEN** each timer action and snapshot policy matches the existing dirty/clean behavior

### Requirement: Canonical Specification Purposes

Canonical specifications SHALL use concrete Purpose text that describes their
actual scope. Archived historical copies remain immutable; this documentation
maintenance requirement does not add a runtime architectural rule.

#### Scenario: Canonical specs are reviewed after archive
- **GIVEN** the canonical specifications produced by Fase 1
- **WHEN** their Purpose sections are inspected
- **THEN** none retains an archive-created placeholder

### Requirement: ChemCalc depends only on Core
M03 SHALL have M00 as its only current ChemUSON runtime dependency and SHALL
not import any M04 module at runtime.

#### Scenario: ChemCalc imports are analyzed
- **WHEN** the architecture test analyzes M03 source files
- **THEN** no M03-to-M04 import or M03 temporary exception is present

### Requirement: ChemIO Does Not Depend on Clean2D

No M01 file SHALL import M02, including a local, lazy, aliased, reexported, or
TYPE_CHECKING import.

#### Scenario: Lazy import is analyzed
- **WHEN** the architecture analyzer inspects every M01 Python file
- **THEN** it finds no import targeting M02 regardless of the import's scope

### Requirement: Clean2D May Consume ChemIO Import Services

M02 SHALL be permitted to consume M01 format parsing and isolated worker
services while remaining independent of GUI modules.

#### Scenario: Imported depiction calls ChemIO
- **WHEN** Clean2D orchestrates SMILES depiction candidates
- **THEN** it may invoke M01 parsing and worker services without creating an M01-to-M02 edge

### Requirement: M01/M02 Cycle Removed

The catalog SHALL contain no M01/M02 circular dependency or any circular
dependency entry.

#### Scenario: Cycle inventory is empty
- **WHEN** the module catalog is inspected after the migration
- **THEN** its circular dependency inventory is empty

### Requirement: Runtime Exception Reduction

The five historical M01-to-M02 runtime exception identities SHALL not
reappear, and the only active exception SHALL be the M01-to-M09 persistence
TYPE_CHECKING crossing.

#### Scenario: Historical runtime edge is reintroduced
- **WHEN** any former M01-to-M02 exception identity is added to the catalog
- **THEN** the frozen exception test fails

### Requirement: ChemIO Is GUI-Independent

M01 SHALL not import M09 at runtime, under `TYPE_CHECKING`, locally, lazily,
through aliases, or through reexports.

#### Scenario: GUI import is inspected in every scope
- **WHEN** the architecture analyzer examines every M01 Python file
- **THEN** it finds no M01-to-M09 import

### Requirement: Persistence Uses Structural Contracts

`PersistenceManager` SHALL accept a GUI-neutral structural contract and SHALL
not name a GUI class or invoke a GUI-private method.

#### Scenario: Pure document is persisted
- **WHEN** a non-Qt fake provides the persistence document operations
- **THEN** save and load complete without importing GUI

### Requirement: No Temporary Exceptions

The module catalog SHALL contain zero `temporary_exceptions` in every module.

#### Scenario: Catalog is inventoried
- **WHEN** architecture tests inspect all module entries
- **THEN** the accumulated temporary-exception list is empty

### Requirement: Eliminated Persistence Exception Cannot Reappear

The former M01-to-M09 `ChemusonCanvas` TYPE_CHECKING identity SHALL be rejected
if it reappears, including normalized file-path or whitespace variants.

#### Scenario: Historical identity is reintroduced
- **WHEN** the frozen exception audit receives the historical identity
- **THEN** it reports the identity as unexpected and M01 as growing from zero
