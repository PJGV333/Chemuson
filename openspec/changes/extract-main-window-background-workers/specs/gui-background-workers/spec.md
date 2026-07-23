## ADDED Requirements

### Requirement: M08 Owns Private Background Workers

M08 SHALL own internal Qt workers for descriptor and Name→Structure operations
outside `main_window.py`. These workers SHALL NOT become public GUI API.

#### Scenario: Worker ownership is inspected
- **GIVEN** the M08 source tree and catalog
- **WHEN** background worker definitions are resolved
- **THEN** both definitions belong to M08 and `ChemusonWindow` remains the public window API

### Requirement: Worker Extraction Preserves Execution Contracts

The extracted workers SHALL preserve their constructor normalization,
`finished` signal payloads, deferred backend imports, timeout values, success
emissions and exception-to-error emissions.

#### Scenario: Worker source is characterized
- **GIVEN** the AST of the worker module
- **WHEN** constructors, signals and `run()` methods are inspected
- **THEN** their observable execution contracts match the pre-extraction implementations

### Requirement: Main Window Retains Thread Orchestration

`ChemusonWindow` SHALL continue to own thread creation, signal connections,
job registries, cancellation, progress UI and result callbacks.

#### Scenario: Main window jobs are inspected
- **GIVEN** descriptor and Name→Structure start handlers
- **WHEN** their orchestration is audited
- **THEN** only worker implementation moved and UI/thread lifecycle remains in the window

### Requirement: Background Workers Do Not Import Main Window

The worker module SHALL NOT import `gui.main_window` or the public window
class.

#### Scenario: Dependency direction is inspected
- **GIVEN** imports below `gui/background_workers.py`
- **WHEN** module dependencies are resolved
- **THEN** no reverse dependency to `main_window` exists
