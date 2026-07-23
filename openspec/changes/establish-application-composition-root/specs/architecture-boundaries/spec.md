## ADDED Requirements

### Requirement: Bootstrap Is a Top-Level Dependency

Modules M00 through M18 SHALL NOT list M19 in current or target dependencies
and SHALL NOT import `chemuson.app` or `chemuson.__main__` at runtime, locally,
lazily, under `TYPE_CHECKING`, through aliases, or through reexports.

#### Scenario: Reverse bootstrap dependency is introduced
- **GIVEN** a source file owned by M00 through M18
- **WHEN** static analysis finds an import from `chemuson.app`
- **THEN** the architecture test fails with the source file and import expression

### Requirement: Non-Graphical CLI Paths Preserve Import Isolation

Architecture tests SHALL verify in a fresh Python process that importing the
entry point and executing non-graphical CLI actions do not load M08, PyQt6 or
the M19 graphical bootstrap.

#### Scenario: CLI import boundary is audited
- **GIVEN** a clean subprocess
- **WHEN** it imports `chemuson.__main__`
- **THEN** `chemuson.gui`, `PyQt6`, and `chemuson.app.bootstrap` are absent from `sys.modules`
