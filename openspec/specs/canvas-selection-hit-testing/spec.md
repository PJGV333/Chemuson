# canvas-selection-hit-testing Specification

## Purpose
TBD - created by archiving change extract-canvas-selection-hit-testing. Update Purpose after archive.
## Requirements
### Requirement: Semantic parent promotion

The query module SHALL promote an item or descendant to its nearest semantic
`CompositeDiagramItem` root and SHALL return `None` when no such root exists.

#### Scenario: descendant promotion
- **GIVEN** a child nested below a composite diagram root
- **WHEN** semantic parent resolution runs
- **THEN** the composite root is returned

### Requirement: Scene item query

The query module SHALL resolve selectable scene items in the existing scene
order, skip disposable orphan atoms through an explicit callback, and map a
text child of an atom to its atom parent.

#### Scenario: orphan exclusion
- **GIVEN** an orphan atom for which the callback returns true
- **WHEN** the scene item query runs
- **THEN** that atom is ignored

### Requirement: Click priority

Click resolution SHALL preserve annotation priority, then molecular atom pick,
then the scene item fallback, then molecular bond fallback.

#### Scenario: atom pick priority
- **GIVEN** a scene item and an atom ID returned by geometric picking
- **WHEN** click resolution runs
- **THEN** the atom item is returned

### Requirement: Selected annotation query

The selected annotation query SHALL consider only items in the scene that are
visible, contain the scene point in their scene bounds and local shape, and
return the item with the greatest z-value.

#### Scenario: removed annotation
- **GIVEN** a selected annotation whose `scene()` is not the active scene
- **WHEN** annotation resolution runs
- **THEN** that item is ignored

