# Amide/Urea 2D Depiction Specification

## Purpose
Specify graph-aware planar placement of amide/urea nitrogen substituents and accurate text clearance for visible atom labels, independently of chemical graph coordinates and global Clean2D.

## ADDED Requirements

### Requirement: Graph-Aware Nitrogen Geometry
The 2D canvas SHALL identify amide/urea nitrogen environments from structural atom/bond connectivity, bond orders, and aromaticity, not from display strings. A neutral amide/urea N attached to a carbonyl carbon by a structural single bond SHALL use a trigonal-planar editorial placement target near 120° when adding automatic substituents. A urea N is an amide N attached to a carbonyl carbon that is also structurally single-bonded to another N. Existing atom coordinates SHALL remain unchanged.

#### Scenario: Sequential urea drawing
- **GIVEN** the user constructs `R–NH–C(=O)–NH–R′` using canvas bond events and short clicks
- **WHEN** a new bond is automatically placed from an N conjugated with the carbonyl
- **THEN** the candidate respects local planar geometry near 120° when feasible
- **AND** the direction avoids unnecessary turns and avoidable crossings
- **AND** the graph's chemistry/connectivity and previously committed coordinates are unchanged

#### Scenario: Distinguish nitrogen environments
- **GIVEN** an amine, amide/urea N, aromatic N, charged/higher-coordinate N, or coordination center
- **WHEN** the canvas selects an automatic local direction
- **THEN** the graph determines which geometry context applies
- **AND** the amide/urea rule is not indiscriminately imposed on unrelated N environments

#### Scenario: Explicit drag
- **GIVEN** a user performs a directed drag from any nitrogen environment
- **WHEN** the bond is committed
- **THEN** the explicit cursor direction and existing fixed-length/modifier semantics are preserved

### Requirement: Preserve Reference Bond Lengths
Automatic placement SHALL preserve configured center-to-center bond length within strict floating-point tolerance. Atom-label size, NH/NH2 hydrogens, charges, substituents, or rendering collisions SHALL NOT shorten, lengthen, or move chemical graph coordinates as a visual workaround.

#### Scenario: Heteroatom-labeled substituent
- **GIVEN** a new bond ends at an atom rendered as N, NH, NH2, or a substituted N
- **WHEN** the gesture commits
- **THEN** its MolGraph center distance remains the configured reference length within `1e-7` units
- **AND** label clearance is handled only by rendering geometry

### Requirement: Label-Path-Aware Bond Trimming
Visible bond trimming SHALL be based on the actual visible label/charge geometry intersecting the bond path, including a bounded optical clearance and the rendered stroke width. Projection-only support of distant glyph points SHALL NOT truncate a bond when the stroke-clearance corridor does not intersect the occupied text shape. Existing styles, caps, double-bond profiles, wedges, hashes, junctions, and optical element anchoring SHALL remain stable outside a demonstrated case.

#### Scenario: NH beside an oblique bond
- **GIVEN** an atom has an `NH`/`NH2` label optically anchored near the N and a neighboring bond leaves at an oblique angle
- **WHEN** the bond's visible path is computed
- **THEN** the path ends near, but does not touch or overlap, the actual glyph shape
- **AND** its remaining visible length is not reduced by remote parts of the label that the path does not approach

#### Scenario: No label-path intersection
- **GIVEN** a label is offset away from a bond ray and the stroke-plus-clearance corridor misses the glyph path
- **WHEN** bond label shrink is computed
- **THEN** that label contributes no projection-only trim to that endpoint
- **AND** center coordinates and all other endpoint/junction trims remain unchanged

### Requirement: Preview, Persistence, and Depiction Invariants
Preview geometry SHALL match the committed MolGraph endpoint. Undo/redo and CMSN save/load SHALL preserve graph IDs, elements, bond orders, aromaticity, and coordinates. The updated trim SHALL render deterministically at zoom/DPR and in light/dark workspace themes without changing established Visual Lab depiction refinements.

#### Scenario: Draw, undo, redo, and reopen
- **GIVEN** a test amide/urea is created by Qt events
- **WHEN** its preview is committed, undone/redone, saved, and reopened
- **THEN** preview/commit endpoints match and the graph snapshot is identical after redo and reload
- **AND** the measured visible label gap remains within its configured optical bounds
