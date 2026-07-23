# application-shell Specification

## Purpose
TBD - created by archiving change extract-application-shell. Update Purpose after archive.
## Requirements
### Requirement: M08 Owns the Application Shell

M08 SHALL own an internal `gui/shell/` package that assembles the main window
regions and their global connections. The shell SHALL NOT replace the public
`ChemusonWindow` class or own desktop process lifecycle.

#### Scenario: Shell ownership is inspected
- **GIVEN** the module catalog and GUI source tree
- **WHEN** the shell package is resolved
- **THEN** it belongs exclusively to M08 and the public window remains at `chemuson.gui.main_window.ChemusonWindow`

### Requirement: Main Window Delegates Shell Assembly Once

`ChemusonWindow.__init__` SHALL perform base Qt initialization and delegate
shell assembly exactly once without retaining a duplicate composition block.

#### Scenario: Constructor seam is audited
- **GIVEN** the AST of `ChemusonWindow.__init__`
- **WHEN** shell delegation is counted
- **THEN** exactly one delegation exists and docks, toolbars, tabs and status widgets are assembled by the shell

### Requirement: Shell Preserves Observable UI Composition

The extracted shell SHALL preserve existing construction order, Qt parents,
dock and toolbar areas, signal connections, initial visibility, initial canvas,
theme application and delayed update check.

#### Scenario: Existing GUI tests instantiate the window
- **GIVEN** an environment with the supported Qt runtime
- **WHEN** `ChemusonWindow` is constructed
- **THEN** its regions, attributes, connections and initial state match the pre-extraction window

### Requirement: Shell Extraction Does Not Move Domain Behavior

The shell package SHALL NOT contain handlers for chemistry, persistence,
documents, templates, validation, compchem or canvas editing.

#### Scenario: Shell source is inspected
- **GIVEN** files below `gui/shell/`
- **WHEN** their definitions and imports are audited
- **THEN** they contain composition only and do not own domain operations or canvas implementation
