# canvas-selection-overlay Specification

## Purpose
TBD - created by archiving change extract-canvas-selection-overlays-and-handles. Update Purpose after archive.
## Requirements
### Requirement: Overlay geometry is consultative

The selection overlay module SHALL calculate padded bounds, transformed handle
positions and screen-space hit metrics without creating or mutating scene items.

#### Scenario: padding and handles
- **GIVEN** a selection bounds rectangle and the existing offsets/radius
- **WHEN** overlay geometry is calculated
- **THEN** the padded rectangle and rotate/move/scale positions match current
  canvas geometry

### Requirement: Handle hit policy is preserved

The module SHALL retain the minimum hit radius, distance priority, stable tie
behavior, invisible/None filtering and existing RuntimeError handling.

#### Scenario: closest handle wins
- **GIVEN** visible handles with different screen-space distances
- **WHEN** handle hit policy runs
- **THEN** the closest eligible handle kind is returned

