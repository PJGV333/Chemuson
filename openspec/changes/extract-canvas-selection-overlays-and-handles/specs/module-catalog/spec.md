## ADDED Requirements

### Requirement: M09 Catalogs Selection Overlay Geometry

M09 SHALL inventory `gui/canvas/selection_overlay.py` and its tests as internal
selection coverage without adding dependencies, exceptions or cycles.

#### Scenario: Overlay ownership is inspected
- **GIVEN** the M09 catalog and source tree
- **WHEN** overlay geometry ownership is enumerated
- **THEN** the new module and tests belong to M09 while the interactive canvas
  remains M09-owned
