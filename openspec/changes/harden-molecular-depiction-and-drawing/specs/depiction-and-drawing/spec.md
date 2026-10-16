# Spec Delta — depiction y dibujo molecular

## ADDED Requirements

### Requirement: Clean2D preserves manually drawn ring integrity

Clean2D SHALL preserve the molecular graph exactly and SHALL NOT accept coordinates that visually open, collapse, or degenerate a manually drawn six-membered ring. If no candidate passes both graph and geometry safety checks, Clean2D SHALL leave the original coordinates unchanged and report its existing controlled result state/reason.

#### Scenario: Manually drawn benzene cycle
- **GIVEN** six manually placed carbon atoms connected by six ring bonds, with aromatic and non-aromatic drawing variants
- **WHEN** quick or publication Clean2D evaluates the complete ring
- **THEN** atom/bond IDs, endpoints, orders/styles, aromatic flags, charges, stereo metadata, and selection SHALL be unchanged
- **AND** any accepted coordinates SHALL keep the six-membered cycle closed, non-degenerate, and free of newly introduced bond crossings

#### Scenario: No safe ring candidate
- **GIVEN** all generated candidates open/collapse the ring or violate invariants
- **WHEN** Clean2D completes
- **THEN** it SHALL reject those candidates and leave all original coordinates unchanged

### Requirement: Wedge endpoints respect atomic-label clearance

Solid and hashed wedges SHALL terminate near the visible boundary of atom labels with a bounded, non-overlapping clearance. Clearance SHALL be computed against the rendered label geometry, not a constant atom-center offset, and SHALL behave consistently for short and wide labels and at supported zoom scales.

#### Scenario: Wedge meets a short or wide label
- **GIVEN** a wedge endpoint next to F, Cl, Br, or a wider atom label
- **WHEN** the bond item is rendered
- **THEN** the wedge SHALL neither intersect the visible glyph shape nor leave an excessive gap outside the measured acceptance range
- **AND** wedge direction, stereo metadata, and hashed-mark count SHALL be preserved

### Requirement: Wedge-to-branch junctions remain continuous and legible

Wedge geometry SHALL integrate with two outgoing bond branches without producing an unsupported rigid base, a gap at the shared atom, an unbounded spike, or obscuring either branch/label. Existing solid-wedge and hashed-wedge stereo orientation SHALL remain unchanged.

#### Scenario: Bifurcating tetrahedral junction
- **GIVEN** one wedge terminates at an atom with two outgoing bonds at varied angles and widths
- **WHEN** the production bond renderer updates its items
- **THEN** local wedge contours SHALL meet the rendered node region naturally and keep both outgoing bond strokes legible
- **AND** the molecular graph and atom coordinates SHALL remain unchanged

### Requirement: Canvas drawing is not blocked by the visible edge

The canvas SHALL allow new atoms, bonds, and ring vertices to be placed at scene positions near or beyond the current visible viewport. It SHALL not clamp user coordinates to the viewport. The scene's navigable bounds SHALL expand or otherwise make committed out-of-viewport geometry reachable while preserving the user's visible focus where possible.

#### Scenario: Draw a bond past the viewport edge
- **GIVEN** the user starts a bond near a visible canvas edge and commits its endpoint beyond the currently visible area
- **WHEN** the drawing operation completes
- **THEN** the graph SHALL contain the requested atom/bond at the chosen scene coordinates
- **AND** the geometry SHALL remain reachable by viewport navigation without prior manual repositioning of the structure

#### Scenario: Undo edge expansion
- **GIVEN** a new structure extends beyond the prior scene bounds
- **WHEN** the user undoes the drawing operation
- **THEN** the molecular graph and canvas navigation bounds SHALL return to a consistent state

### Requirement: Ring and Haworth drawing remains undoable and chemically faithful

Ring/Haworth templates and any local snapping/normalization SHALL preserve intended ring connectivity, atom/bond semantics, stereochemical substituent orientation, and the existing undo/redo contract. Geometry assistance SHALL not mutate unrelated atoms or silently alter a user-placed substituent.

#### Scenario: Place or snap a ring/Haworth template
- **GIVEN** a user intentionally places a ring or Haworth template near an existing structure
- **WHEN** placement/snapping is committed
- **THEN** the template SHALL be geometrically coherent and connected only at the intended anchor(s)
- **AND** one undo SHALL restore the prior graph and coordinates
