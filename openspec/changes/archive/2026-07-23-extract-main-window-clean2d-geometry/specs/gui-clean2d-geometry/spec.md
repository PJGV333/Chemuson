## ADDED Requirements

### Requirement: M08 Owns Pure Clean2D Geometry Integration

M08 SHALL own the GUI integration geometry used to normalize, align and
complete Clean2D coordinates outside `main_window.py`. The implementation
SHALL remain internal and independent from Qt.

#### Scenario: Geometry ownership is inspected
- **GIVEN** the M08 source tree
- **WHEN** Clean2D integration transformations are resolved
- **THEN** their implementations belong to `gui/clean2d_geometry.py` and are not defined by `ChemusonWindow`

### Requirement: Geometry Extraction Preserves Numeric Contracts

The extracted functions SHALL preserve geometric centers, mean bond lengths,
bounded scaling, rigid alignment with optional reflection, and missing
hydrogen projection results.

#### Scenario: Existing regressions execute
- **GIVEN** the scaling, alignment and hydrogen projection regression inputs
- **WHEN** the extracted implementation is invoked through the compatibility aliases
- **THEN** all existing numeric expectations remain unchanged

### Requirement: Main Window Preserves Private Compatibility Aliases

`ChemusonWindow` SHALL retain its five existing private helper names as static
aliases during this phase so current controllers and tests do not change
their call sites.

#### Scenario: Private consumers are inspected
- **GIVEN** current calls through `ChemusonWindow` instances or the class
- **WHEN** the five helper names are resolved
- **THEN** they invoke the extracted pure functions without binding a window argument

### Requirement: Geometry Module Has No GUI Dependencies

The geometry module SHALL NOT import PyQt, `gui.main_window`, canvas modules or
domain backends.

#### Scenario: Imports are inspected
- **GIVEN** the AST of `gui/clean2d_geometry.py`
- **WHEN** its imports are enumerated
- **THEN** only Python standard-library imports are present
