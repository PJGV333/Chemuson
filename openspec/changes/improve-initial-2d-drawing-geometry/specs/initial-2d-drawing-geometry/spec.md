# Initial 2D Drawing Geometry Specification

## Purpose
Define the editorial geometry used when ChemUSON creates new skeletal structures in the 2D canvas. Keep it independent from physical 3D geometry and from raster appearance.

## ADDED Requirements

### Requirement: Editorial Acyclic Skeletal Angles
The canvas SHALL use a deterministic 120° interior-angle convention for default sp3 skeletal growth and acyclic zigzag chains. The 109.5° tetrahedral reference SHALL remain distinct and SHALL NOT be substituted into 2D layout defaults. The policy SHALL NOT change the global `angle_step_deg` setting.

#### Scenario: New sp3 bond by short click
- **GIVEN** an atom with one or more structural neighbors and no intentional pointer-directed drag
- **WHEN** the user creates a default bond by a short click
- **THEN** the selected local direction is approximately 120° from the adjacent skeletal direction when geometrically feasible
- **AND** existing atom coordinates remain unchanged

#### Scenario: Free acyclic chain
- **GIVEN** fixed-angle chain drawing is enabled
- **WHEN** the user previews and commits an acyclic chain
- **THEN** all new bonds have equal scene-space lengths and internal zigzag angles are 120° within numeric tolerance
- **AND** preview coordinates equal committed atom coordinates

#### Scenario: Explicit pointer direction
- **GIVEN** a bond is drawn by a directed drag
- **WHEN** the cursor or `Alt`/`Shift` modifiers express an orientation/length choice
- **THEN** the existing drag, snapping, fixed-length, and modifier semantics remain in effect
- **AND** the system does not silently realign previously drawn atoms

### Requirement: Geometry-Aware Automatic Placement
An automatic short-click suggestion SHALL be computed from local bond vectors, molecular context, and collision constraints in scene coordinates. It SHALL NOT quantize that geometry to an unrelated absolute screen/scene angle grid when no pointer direction was expressed.

#### Scenario: Substituent on a rotated benzene vertex
- **GIVEN** a regular benzene ring rotated to an arbitrary scene angle and a feasible isolated carbon vertex
- **WHEN** a new exocyclic single bond is placed by short click
- **THEN** the bond points radially outward and the two local angles are approximately 120°
- **AND** the existing six ring vertices are unchanged

#### Scenario: Collision-aware candidate
- **GIVEN** one or more automatically suggested directions intersect existing atoms or bonds
- **WHEN** a new bond direction is chosen
- **THEN** the deterministic valid local candidate with the best existing cursor/environment score is selected, or the least-colliding fallback is used without moving existing atoms

### Requirement: Regular Ring Construction
A regular n-member ring SHALL be generated as a regular polygon with uniform scene-space side lengths, exact polygon closure within floating-point tolerance, and internal angle `180° - 360°/n`. Rotation, translation, zoom, and device scale SHALL NOT alter model coordinates or polygon regularity.

#### Scenario: Isolated benzene
- **GIVEN** the ring tool creates an isolated six-member ring with side length L
- **WHEN** its model coordinates are measured before rendering
- **THEN** all six side lengths equal L within numeric tolerance
- **AND** all six internal angles equal 120° within numeric tolerance
- **AND** the cycle closes without duplicate vertices or bonds

#### Scenario: Other regular ring sizes
- **GIVEN** the ring tool creates a regular three-, four-, or five-member ring
- **WHEN** its model coordinates are measured
- **THEN** its internal angles are 60°, 90°, or 108° respectively and all side lengths are equal

#### Scenario: Ring anchored to an existing edge
- **GIVEN** an existing bond anchors a new ring
- **WHEN** the new ring is committed
- **THEN** the existing edge remains fixed, the new polygon is regular at that edge length, and atom/bond fusion occurs only at matching existing vertices/edges

### Requirement: Model Geometry and Rendering Are Separate
The canvas renderer SHALL render atom and bond items from model coordinates without changing those coordinates. Junction trimming, label clearance, double/aromatic offsets, antialiasing, and device scaling SHALL be evaluated as visual effects, not used as evidence of incorrect model geometry.

#### Scenario: Geometry audit
- **GIVEN** a generated molecule and the corresponding rendered image
- **WHEN** its quality is assessed
- **THEN** coordinate-based lengths, angles, closure, and topology are reported separately from pixel appearance
- **AND** existing editorial rendering styles and Clean2D remain unchanged unless a distinct visual defect is independently reproduced

### Requirement: Short Click Creates One New Bond Undoably
A short click beginning on an existing atom SHALL create the automatically suggested new bond/atom when no target atom is snapped. Release hit-testing SHALL NOT reinterpret the start atom's adjacent bond as the user's target bond. Direct clicks on an existing bond SHALL retain their current order/style behavior.

#### Scenario: Short click at a ring atom
- **GIVEN** an aromatic ring atom with two adjacent ring bonds and no free target atom at the suggested endpoint
- **WHEN** the user presses and releases within the short-click threshold
- **THEN** exactly one new atom and one bond are added at the suggested endpoint
- **AND** the ring topology, bond order, and aromaticity are unchanged
- **AND** one undo restores the exact initial model and redo restores the same IDs and coordinates
