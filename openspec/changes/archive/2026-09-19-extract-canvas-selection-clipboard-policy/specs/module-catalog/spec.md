## ADDED Requirements

### Requirement: M09 Catalogs Selection Clipboard Policy

M09 SHALL inventory `gui/canvas/selection_clipboard.py` and its tests as
internal clipboard policy coverage without changing its interactive canvas
ownership or adding dependencies, exceptions or cycles.

#### Scenario: Clipboard ownership is inspected
- **GIVEN** the M09 catalog and source tree
- **WHEN** clipboard policy ownership is enumerated
- **THEN** the new module and tests belong to M09 while paste coordination
  remains in the existing canvas mixins
