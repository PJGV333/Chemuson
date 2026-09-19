# editor2d-selection-ownership Specification

## Purpose
TBD - created by archiving change move-canvas-selection-helpers-to-editor2d. Update Purpose after archive.
## Requirements
### Requirement: M20 Owns Canonical Selection Helpers

M20 SHALL own the canonical implementations of `selection_geometry`,
`selection_bounds`, `selection_hit_testing`, `selection_overlay` and
`selection_clipboard` under `src/chemuson/gui/editor2d/selection/`. The parent
`src/chemuson/gui/editor2d/__init__.py` SHALL remain a logic-free namespace
owned outside M20.

#### Scenario: Canonical helper paths exist
- **WHEN** the M20 source package is inspected
- **THEN** all five helper modules exist under `gui/editor2d/selection/` and
  contain the implementations used by canvas mixins

### Requirement: Legacy Canvas Imports Remain Compatible

The five historical `gui.canvas` helper paths SHALL remain importable as thin
compatibility shims that re-export the corresponding M20 module without owning
function or class definitions.

#### Scenario: Existing consumers import through legacy paths
- **WHEN** a consumer imports a helper from `chemuson.gui.canvas`
- **THEN** the same public helper symbols are available from the corresponding
  `chemuson.gui.editor2d.selection` implementation

### Requirement: M09 and M20 Have Exclusive Ownership

M09 SHALL own canvas mixins and depend on M20 for the extracted helper policy.
M20 SHALL own only `gui/editor2d/selection/`, SHALL declare no ChemUSON module
dependencies, and SHALL have no temporary exceptions or circular dependencies.
M20 SHALL NOT depend on M09.

#### Scenario: Catalog ownership is audited
- **WHEN** `architecture/modules.yml` is inspected
- **THEN** M20 is named `gui.editor2d.selection`, owns only the selection
  path and its tests, M09 lists M20 as current and target dependency, and the
  parent namespace is assigned without path overlap

