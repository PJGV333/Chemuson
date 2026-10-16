# Spec Delta — uniones de wedges, dobles y etiquetas

## ADDED Requirements

### Requirement: Hashed wedges integrate legibly with adjacent double bonds

A hashed stereobond meeting an atom that also carries a double bond SHALL meet the shared atom region without an abrupt gap, uncontrolled crossing, or collision with a secondary pi stroke. The renderer SHALL distinguish the primary sigma stroke from the offset pi stroke; it SHALL preserve hashed mark count, direction, bond order, stereo metadata, and atom coordinates.

#### Scenario: Hashed wedge approaches a nearby sigma/pi double
- **GIVEN** a hashed wedge and a double-bond neighbor at varied angles, including near-parallel pi geometry
- **WHEN** production bond items are rendered
- **THEN** hash marks SHALL terminate/transition consistently at the shared-node contour without colliding with the visible pi stroke
- **AND** the sigma/pi depiction and graph identity SHALL remain unchanged

### Requirement: Labels are optically centered in double-bond contexts

Atom labels adjacent to double bonds SHALL be visually centered relative to their anchor atom and laid out with awareness of the rendered sigma/pi orientation and neighboring bond direction. Any adjustment SHALL be limited to demonstrated optical error and SHALL preserve readable spacing for common heteroatoms and halogens.

#### Scenario: Heteroatom label beside a double bond
- **GIVEN** an O, N, halogen, or short label in horizontal, vertical, or diagonal double-bond geometry
- **WHEN** the depiction is refreshed
- **THEN** the visible label centroid SHALL remain within a bounded tolerance of its intended optical anchor
- **AND** no sigma, pi, charge, or neighboring label glyph SHALL be obscured

### Requirement: Wedge-label clearance is consistent and bounded

Solid and hashed wedges SHALL use a coherent clearance policy based on the rendered wedge/mark stroke and visible label contour. The policy SHALL prevent overlap while avoiding unnecessary separation and SHALL be checked on both wedge tips and bases.

#### Scenario: Solid and hashed wedge meet a label
- **GIVEN** short and wide labels at either end of solid/hashed wedges
- **WHEN** the production items are rendered
- **THEN** the measured visible gap SHALL remain inside the accepted minimum/maximum range for the relevant stroke width
- **AND** no wedge fill, hash mark, or cap SHALL intersect the label/charge silhouette

### Requirement: Solid wedges join bold Haworth edges cleanly

Solid wedge contours SHALL meet bold edges in cyclic Haworth-like depictions without bulky corner blobs, unbounded spikes, gaps wider than the participating stroke widths, or polygon self-intersection.

#### Scenario: Wedge meets bold edge in pentagon and hexagon
- **GIVEN** a solid wedge incident to a bold ring edge in a pentagonal or hexagonal Haworth depiction
- **WHEN** the production canvas is rendered
- **THEN** the local contours SHALL appear continuous and bounded at matching scale in light and dark themes
- **AND** ring connectivity, bold semantics, wedge orientation, and all atom coordinates SHALL be preserved

### Requirement: Existing depiction controls do not regress

Refinements SHALL NOT degrade ordinary double bonds, unlabeled wedges, common heteroatom labels, simple rings, or depiction cases already accepted by the molecular-depiction-excellence campaign.

#### Scenario: Standard depiction controls
- **GIVEN** ordinary doubles, unlabelled solid/hashed wedges, simple rings, F/Cl/Br/O/N labels, and prior junction/Haworth cases
- **WHEN** rendering and geometry regression tests run
- **THEN** their established geometry, graph invariants, and bounded clearance SHALL remain intact
