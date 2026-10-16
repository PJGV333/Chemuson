# Spec Delta — geometría visual de uniones de enlaces

## ADDED Requirements

### Requirement: Bond junctions preserve molecular identity

Bond-junction rendering SHALL be a deterministic presentation-only transformation. It MUST NOT mutate atom identities, coordinates, connectivity, bond orders/styles, charges, aromaticity, stereochemistry, document data, or undo history.

#### Scenario: Rebuild items leaves the molecular graph unchanged
- **GIVEN** a molecular graph containing wedges, multiple bonds, labels, charges, and rings
- **WHEN** the canvas rebuilds or refreshes its bond graphics items
- **THEN** every graph field and molecular coordinate SHALL remain byte-for-byte/equality unchanged
- **AND** wedge/hash endpoint orientation and `bold` style semantics SHALL remain unchanged

### Requirement: Junction profiles retain individual stroke geometry

A junction calculation that handles a multiple bond SHALL represent its sigma and pi strokes separately, including their actual offsets, widths, trims, and caps. It SHALL NOT replace a double bond with one centerline and an enlarged scalar width. Solid and hashed wedge contours SHALL retain their distinct visual profiles.

#### Scenario: Solid wedge meets a double bond
- **GIVEN** a solid wedge meeting a double bond at a shared atom
- **WHEN** the join is rendered
- **THEN** the sigma stroke SHALL connect without an unjustified peak or gap
- **AND** the pi stroke SHALL retain its independent offset and intended end trim
- **AND** neither line SHALL be displaced into an atom label

#### Scenario: Hashed wedge meets a bold bond
- **GIVEN** a hashed wedge meeting a thick `bold` bond at its wide end
- **WHEN** the join is rendered
- **THEN** the hash marks SHALL remain monotonically ordered in width with uniform longitudinal spacing
- **AND** the terminal contact SHALL have no isolated protrusion beyond the connected stroke envelope
- **AND** the hashed orientation and `bold` semantics SHALL remain unchanged

### Requirement: Wedge geometry remains bounded and label-safe

Solid wedge geometry SHALL remain finite, non-self-intersecting, and bounded by the local bond length and stroke profile. A junction adjustment SHALL NOT enter the visible shape of a neighboring atom label, charge, or explicit hydrogen.

#### Scenario: Junction-angle regression matrix
- **GIVEN** junctions at 30°, 60°, 90°, 120°, 150°, and near-collinear angles, at multiple scales and widths
- **WHEN** wedge corners are adapted to neighboring strokes
- **THEN** the polygon SHALL preserve the narrow/wide endpoint ordering
- **AND** no corner SHALL form an unbounded miter, detached spike, or self-intersection

#### Scenario: Labeled heteroatom junction
- **GIVEN** an N/O/S atom with charge and/or explicit hydrogen at a wedge junction
- **WHEN** the scene is rendered
- **THEN** the molecular label, charge, and hydrogen glyphs SHALL remain unobscured

### Requirement: Junction comparisons are reproducible across rendering surfaces

The curated junction gallery SHALL be reproducible from source and render the same molecular examples in light and dark application contexts. Geometry SHALL be deterministic across zoom and vector export; themes SHALL affect only the applicable color/surface tokens, not molecular geometry.

#### Scenario: Theme, zoom, and export comparison
- **GIVEN** the same molecule, style, and coordinates
- **WHEN** rendered at supported zoom scales, in light/dark context, and through vector export where available
- **THEN** atom positions and junction topology SHALL be unchanged
- **AND** only expected antialiasing/device-pixel differences SHALL occur
