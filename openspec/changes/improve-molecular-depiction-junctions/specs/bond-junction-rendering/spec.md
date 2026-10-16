# Spec Delta — geometría visual de uniones de enlaces

## ADDED Requirements

### Requirement: Bond-junction rendering preserves molecular identity

Bond-junction rendering SHALL be a deterministic presentation-only transformation. It MUST NOT mutate atom identities, coordinates, connectivity, bond orders/styles, charges, aromaticity, stereochemistry, document data, or undo history.

#### Scenario: Rebuild items leaves the molecular graph unchanged
- **GIVEN** a molecular graph containing wedges, multiple bonds, labels, charges, and rings
- **WHEN** the canvas rebuilds or refreshes its bond graphics items
- **THEN** every graph field and molecular coordinate SHALL remain unchanged
- **AND** wedge/hash endpoint orientation and `bold` style semantics SHALL remain unchanged

### Requirement: Hashed-wedge marks use a centered, bounded profile

Hashed-wedge marks SHALL be centered at equal intervals along the existing bond segment. Their widths SHALL vary monotonically from narrow to wide, interpolated at each mark center. Geometry SHALL be finite, deterministic, and independent of molecular graph mutation or Qt rasterization.

#### Scenario: Hashed wedge meets a simple or bold bond
- **GIVEN** a hashed wedge meeting a simple or thick `bold` bond at a shared atom
- **WHEN** its transverse marks are rendered
- **THEN** mark centers SHALL have uniform longitudinal spacing and symmetric half-interval margins at the two endpoints
- **AND** mark widths SHALL grow monotonically without changing mark count or bond orientation
- **AND** no terminal mark SHALL be centered exactly on the shared atom
- **AND** `bold` SHALL remain a non-stereochemical stroke

### Requirement: Existing solid-wedge geometry remains bounded and label-safe

Solid wedge geometry SHALL remain finite, non-self-intersecting, and bounded by the local bond length and stroke profile. This change SHALL NOT alter wedge contours or enter the visible shape of neighboring atom labels, charges, or explicit hydrogens.

#### Scenario: Junction-angle regression matrix
- **GIVEN** solid-wedge junctions at 30°, 60°, 90°, 120°, 150°, and near-collinear angles, at multiple lengths and stroke widths
- **WHEN** their existing contours are generated
- **THEN** the narrow/wide endpoint ordering SHALL be preserved
- **AND** no polygon SHALL form an unbounded miter or self-intersection

#### Scenario: Labeled heteroatom junction
- **GIVEN** an N/O/S atom with charge and/or explicit hydrogen near a wedge junction
- **WHEN** the gallery is rendered
- **THEN** its existing label glyphs SHALL remain visible and the molecule's graph SHALL be unchanged

### Requirement: Junction comparisons are reproducible across rendering surfaces

The curated gallery SHALL be reproducible from source and render the same molecular examples in light and dark application contexts. It SHALL use the production Qt graphics items; theme changes SHALL affect page/surface colors without changing molecular geometry. The optional comparison sheet SHALL show matching before/after panels at equal dimensions.

#### Scenario: Theme and comparison rendering
- **GIVEN** the same molecule, style, and coordinates
- **WHEN** the gallery is rendered in light/dark context and composed with its baseline
- **THEN** atom positions and junction topology SHALL be unchanged
- **AND** only expected theme and antialiasing differences SHALL occur
