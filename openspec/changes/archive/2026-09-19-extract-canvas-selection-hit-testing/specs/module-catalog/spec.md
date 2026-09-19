## ADDED Requirements

### Requirement: M09 Catalogs Selection Hit Testing

M09 SHALL inventory `gui/canvas/selection_hit_testing.py` and its functional and
architectural tests as internal consultative selection coverage without adding
dependencies, temporary exceptions or circular dependencies.

#### Scenario: Hit-testing ownership is inspected
- **GIVEN** the M09 catalog and source tree
- **WHEN** selection hit-testing ownership and tests are enumerated
- **THEN** the new query module and tests belong to M09 and the existing M09
  dependency sets remain unchanged
