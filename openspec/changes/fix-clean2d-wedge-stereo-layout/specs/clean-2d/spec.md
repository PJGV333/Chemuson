# Spec Delta — Clean2D con enlaces en cuña sólida y rayada

## Purpose

Fix the stereo-layout preservation check so depictions that contain solid or hashed wedge bonds can be cleaned, while keeping every real stereochemical inversion rejected.

## ADDED Requirements

### Requirement: Clean 2D preserves the depicted wedge and hash configuration
Clean 2D SHALL determine stereochemical risk from the cyclic order of the real structural neighbours of each stereo center together with each neighbour's bond style and explicit bond stereochemistry. A candidate that keeps that cyclic order and those styles SHALL NOT be treated as a stereochemical risk. A candidate that inverts it SHALL be rejected.

#### Scenario: Single solid wedge depiction is cleaned
- **GIVEN** an acyclic depiction whose only stereo feature is one `WEDGE` bond
- **WHEN** Clean 2D runs in `quick` or `publication` mode over the whole depiction
- **THEN** an accepted candidate SHALL be selected and at least one coordinate SHALL change
- **AND** the wedge bond, its endpoints, and its style SHALL remain unchanged
- **AND** the cyclic neighbour configuration of the stereo center SHALL be unchanged

#### Scenario: Single hashed wedge depiction is cleaned
- **GIVEN** an acyclic depiction whose only stereo feature is one `HASHED` bond
- **WHEN** Clean 2D runs in `quick` or `publication` mode over the whole depiction
- **THEN** an accepted candidate SHALL be selected and at least one coordinate SHALL change
- **AND** the hashed bond style and the cyclic neighbour configuration SHALL be unchanged

#### Scenario: Tetrahedral center with wedge and hash is cleaned
- **GIVEN** a center carrying one `WEDGE` and one `HASHED` bond to different neighbours
- **WHEN** Clean 2D runs over the whole depiction
- **THEN** an accepted candidate SHALL be selected
- **AND** both wedge styles and the cyclic neighbour configuration SHALL be unchanged

#### Scenario: Order-preserving relayout is not a stereo risk
- **GIVEN** a stereo center and a candidate that rotates the whole neighbourhood or changes the angular gaps between neighbours without changing their cyclic order
- **WHEN** Clean 2D evaluates stereo preservation
- **THEN** the stereo layout signature SHALL be equal before and after
- **AND** the candidate SHALL NOT be rejected as a stereo risk

#### Scenario: Rotation and scale do not change the signature
- **GIVEN** a depiction with a stereo center and a rigid rotation of its neighbourhood
- **WHEN** the stereo layout signature is computed for both layouts
- **THEN** both signatures SHALL be equal

#### Scenario: Inverted arrangements are rejected
- **GIVEN** a stereo center and a candidate that mirrors the neighbourhood, swaps two neighbours, or changes a `WEDGE` bond into a `HASHED` bond
- **WHEN** Clean 2D evaluates stereo preservation
- **THEN** the stereo layout signature SHALL differ
- **AND** the candidate SHALL be rejected as a stereo risk

### Requirement: Stereo signature covers selection boundary neighbours
For a partial selection, every real structural neighbour of a selected stereo center SHALL participate in the stereo layout signature, including neighbours outside the selection, using their graph coordinates.

#### Scenario: Boundary neighbour orders a partially selected center
- **GIVEN** a stereo center inside the selection with at least one structural neighbour outside it
- **WHEN** the stereo layout signature is computed
- **THEN** the outside neighbour SHALL appear in that center's cyclic order

#### Scenario: Crossing a boundary neighbour inverts the configuration
- **GIVEN** a selected neighbour of a partially selected stereo center
- **WHEN** a candidate moves that neighbour across the ray of a neighbour outside the selection
- **THEN** the stereo layout signature SHALL differ
- **AND** the candidate SHALL be rejected as a stereo risk

#### Scenario: Partial selection never applies an inversion
- **GIVEN** a selection that leaves a stereo center or a wedge bond on the selection boundary
- **WHEN** Clean 2D runs in `quick` mode
- **THEN** the result state SHALL be one of `applied`, `no-op`, `preserve-only`, or `failed-controlled`
- **AND** the depicted cyclic configuration and wedge styles SHALL be unchanged for every returned state

### Requirement: Wedge cleaning is deterministic
Cleaning a depiction that contains wedge bonds SHALL produce the same result state, candidate source, coordinates, and geometry hash for the same input graph and options.

#### Scenario: Same input twice
- **GIVEN** the same wedge depiction and the same Clean 2D mode and target bond length
- **WHEN** Clean 2D runs twice on equivalent graphs
- **THEN** the result state, selected coordinates, and geometry hash SHALL be identical

## MODIFIED Requirements

### Requirement: Clean 2D enforces chemical and canvas safety categories
Clean 2D SHALL reject candidates that risk chemical identity, stereochemical meaning, or canvas boundary integrity.

#### Scenario: Invariant violation rejected
- **GIVEN** a candidate changes anything other than allowed target coordinates
- **WHEN** Clean 2D evaluates invariants
- **THEN** the candidate SHALL be rejected with `invariant-violation`

#### Scenario: Stereo risk rejected
- **GIVEN** a candidate risks changing the visual or metadata interpretation of stereochemistry
- **WHEN** Clean 2D evaluates stereo preservation
- **THEN** the candidate SHALL be rejected with `stereo-risk`

#### Scenario: Geometric relayout is not stereo risk
- **GIVEN** a candidate that keeps the cyclic neighbour order and the wedge/hash styles of every stereo center
- **WHEN** Clean 2D evaluates stereo preservation
- **THEN** the candidate SHALL NOT be rejected with `stereo-risk` on geometric grounds alone

#### Scenario: Boundary bond risk rejected
- **GIVEN** Clean 2D targets a selection with bonds crossing the selection boundary
- **WHEN** a candidate would distort or disconnect those boundary bonds beyond accepted limits
- **THEN** the candidate SHALL be rejected with `boundary-bond-risk`
