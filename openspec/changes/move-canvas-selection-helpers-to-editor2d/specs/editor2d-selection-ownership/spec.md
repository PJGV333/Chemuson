## ADDED Requirements

### Requirement: M20 Owns Canonical Selection Helpers

M20 SHALL own the canonical implementations of `selection_geometry`,
`selection_bounds`, `selection_hit_testing`, `selection_overlay` and
`selection_clipboard` under `src/chemuson/gui/editor2d/`.

#### Scenario: Canonical helper paths exist
- **WHEN** the M20 source package is inspected
- **THEN** all five helper modules exist under `gui/editor2d/` and contain the
  implementations used by canvas mixins

### Requirement: Legacy Canvas Imports Remain Compatible

The five historical `gui.canvas` helper paths SHALL remain importable as thin
compatibility shims that re-export the corresponding M20 module without owning
function or class definitions.

#### Scenario: Existing consumers import through legacy paths
- **WHEN** a consumer imports a helper from `chemuson.gui.canvas`
- **THEN** the same public helper symbols are available from the M20
  implementation

### Requirement: M09 and M20 Have Exclusive Ownership

M09 SHALL own canvas mixins and depend on M20 for the extracted helper policy.
M20 SHALL own only the canonical editor2d helper package and SHALL declare no
ChemUSON module dependencies, temporary exceptions or circular dependencies.

#### Scenario: Catalog ownership is audited
- **WHEN** `architecture/modules.yml` is inspected
- **THEN** the five helpers are listed only in M20 internal API, M20 owns their
  tests, and M09 lists M20 in current and target dependencies
