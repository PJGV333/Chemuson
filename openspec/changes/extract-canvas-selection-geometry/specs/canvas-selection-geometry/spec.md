## ADDED Requirements

### Requirement: M09 Owns Selection Value Geometry

M09 SHALL own deterministic selection geometry and normalization rules outside
`canvas_selection.py`. The implementation SHALL remain internal and free of
canvas, scene, model, command and graphics-item state.

#### Scenario: Geometry ownership is inspected
- **GIVEN** the M09 source tree
- **WHEN** the seven selection value rules are resolved
- **THEN** their implementations belong to `gui/canvas/selection_geometry.py`

### Requirement: Extraction Preserves Numeric Contracts

The extracted functions SHALL preserve signed angular deltas, point rotation,
anchor scaling, optional-float and point equality, label-scale normalization
and custom-stroke normalization exactly.

#### Scenario: Numeric regressions execute
- **GIVEN** representative values at each tolerance and normalization boundary
- **WHEN** the extracted functions are invoked
- **THEN** their results match the pre-extraction formulas

### Requirement: Canvas Preserves Private Compatibility Names

`CanvasSelectionMixin` SHALL retain the seven existing private helper names.
Six SHALL be static aliases and custom-stroke normalization SHALL retain its
one-argument instance signature through a minimal delegating wrapper.

#### Scenario: Mixin consumers are inspected
- **GIVEN** current calls from canvas selection, text and structure mixins
- **WHEN** the seven private names are resolved through `ChemusonCanvas`
- **THEN** no consumer call site or MRO contract changes

### Requirement: Geometry Module Avoids Stateful GUI Dependencies

The geometry module SHALL NOT import QtGui, QtWidgets, canvas modules, domain
models, commands or graphics items. It MAY use `QPointF` from QtCore as its
point value type.

#### Scenario: Imports are inspected
- **GIVEN** the AST of `gui/canvas/selection_geometry.py`
- **WHEN** its imports are enumerated
- **THEN** only standard-library imports and `PyQt6.QtCore.QPointF` are present
