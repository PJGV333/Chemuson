## MODIFIED Requirements

### Requirement: M15 Owns Only Historical Compatibility Shims

M15 SHALL own only the historical `utils` compatibility modules. Canonical
settings/resources SHALL belong to M21, and canonical autosave/crash/recovery
implementations SHALL belong to M22. M15's direct dependencies SHALL be limited
to the M21 and M22 shims it contains.

#### Scenario: Utility ownership is inspected

- **WHEN** M15 paths, imports and internal APIs are audited
- **THEN** no canonical implementation is assigned to M15 and its compatibility
  dependencies are explicitly recorded

### Requirement: M19 Catalogs the Composition Root

M19 SHALL own `src/chemuson/__main__.py` and `src/chemuson/app/`. Its current and
target dependencies SHALL be exactly M08, M18 and M22, matching its direct
imports and crash-hook installation. Its public API SHALL continue to contain
`main`.

#### Scenario: Composition-root dependencies are inspected

- **WHEN** M19 paths, dependencies and public API are audited
- **THEN** both bootstrap paths are M19-owned, dependencies are exactly M08/M18/M22,
  and `main` resolves statically

### Requirement: M20 and M09 Distinguish Canonical Selection from Shims

M20 SHALL own only `src/chemuson/gui/editor2d/selection/` and its five canonical
selection helper modules. M09 SHALL own the historical `src/chemuson/gui/canvas/`
shim paths and depend on M20. The M20 parent namespace SHALL remain owned by M08.

#### Scenario: Selection ownership is inspected

- **WHEN** canonical and historical selection paths are audited
- **THEN** canonical implementations belong exclusively to M20, canvas shims
  remain import-only under M09, and the dependency direction is M09 to M20

### Requirement: M21 Owns Explicit Platform Files

M21 SHALL own `src/chemuson/platform/__init__.py`, `settings.py` and
`resources.py` explicitly. M21 SHALL have no ChemUSON dependencies and SHALL NOT
import GUI widgets or controllers.

#### Scenario: Platform paths are inspected

- **WHEN** M21 paths and imports are audited
- **THEN** the three platform files are listed explicitly, no package directory
  wildcard is used, and no GUI dependency is present

### Requirement: M22 Owns Recovery Filesystem Policy

M22 SHALL own the canonical `recovery.py` implementation and the existing
`autosave.py` and `crash_reporter.py` implementations. The recovery module SHALL
provide `read_autosave_metadata`, `list_autosave_entries` and
`archive_autosave` without importing `chemuson.gui` or `PersistenceManager`.

#### Scenario: Recovery ownership is inspected

- **WHEN** recovery source, controller consumers and catalog entries are audited
- **THEN** the three filesystem-policy functions have one canonical M22 owner,
  RecoveryController retains only delegation plus UI/document orchestration, and
  M22 has no ChemUSON module dependencies

## ADDED Requirements

### Requirement: Global OpenSpec Validation Is Clean

All active and archived specifications SHALL satisfy the normative requirement
syntax accepted by strict OpenSpec validation.

#### Scenario: Strict validation runs

- **WHEN** `openspec validate --all --strict` is executed
- **THEN** validation completes without failures
