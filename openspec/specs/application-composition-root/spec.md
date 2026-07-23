# application-composition-root Specification

## Purpose
TBD - created by archiving change establish-application-composition-root. Update Purpose after archive.
## Requirements
### Requirement: M19 Owns Application Composition

M19 SHALL own the creation and lifecycle of the desktop application in
`src/chemuson/app/`. M08 SHALL own `ChemusonWindow` and SHALL NOT own the
application bootstrap function.

#### Scenario: Bootstrap ownership is inspected
- **GIVEN** the module catalog and Python sources
- **WHEN** ownership and definitions are analyzed
- **THEN** the application bootstrap belongs to M19 and `gui/main_window.py` does not define `run_app`

### Requirement: CLI Remains Lightweight

Importing `chemuson.__main__`, requesting `--version`, or requesting CLI help
SHALL NOT import PyQt6, `chemuson.gui`, or the graphical bootstrap module.

#### Scenario: Version is requested without Qt
- **GIVEN** a fresh Python process where graphical imports are observable
- **WHEN** `chemuson --version` is executed
- **THEN** it prints the application version and no GUI or PyQt6 module is loaded

#### Scenario: Entry point module is imported
- **GIVEN** a fresh Python process
- **WHEN** `chemuson.__main__` is imported without invoking graphical startup
- **THEN** the import succeeds without loading the GUI stack

### Requirement: CLI Delegates Once to the Composition Root

When no non-graphical CLI action terminates execution, `main()` SHALL lazily
load the M19 bootstrap and invoke it exactly once.

#### Scenario: Default invocation starts the GUI
- **GIVEN** the graphical bootstrap is replaced by a recording fake
- **WHEN** the CLI runs without `--version`
- **THEN** the fake bootstrap is imported after argument parsing and called exactly once

### Requirement: Startup Lifecycle Is Behaviorally Preserved

The M19 bootstrap SHALL preserve the existing startup order: install crash
reporting, create `QApplication`, set application name and version, create
`ChemusonWindow`, check autosaves, show the window, run the event loop, and
terminate with its exit code.

#### Scenario: Successful startup order is observed
- **GIVEN** recording substitutes for Qt, the window, version and crash reporter
- **WHEN** the bootstrap runs successfully
- **THEN** every startup operation occurs once in the established order and the event-loop code is passed to process termination

### Requirement: Startup Failure Behavior Is Preserved

The M19 bootstrap SHALL preserve the existing crash-log creation and user
notification behavior when startup raises an exception.

#### Scenario: Failure occurs with a Qt application instance
- **GIVEN** startup raises after a Qt application instance exists
- **WHEN** the bootstrap handles the exception
- **THEN** it writes one crash log and shows the existing critical dialog text

#### Scenario: Failure occurs without a Qt application instance
- **GIVEN** startup raises before a Qt application instance exists
- **WHEN** the bootstrap handles the exception
- **THEN** it writes one crash log and writes the existing failure text and log path to stderr
