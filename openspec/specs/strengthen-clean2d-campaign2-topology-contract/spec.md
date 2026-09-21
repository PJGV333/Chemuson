# strengthen-clean2d-campaign2-topology-contract Specification

## Purpose
This capability strengthens deterministic, JSON-safe topology evidence for
ring systems, block adjacency, connectors, articulation structure, branch
points, and reusable fused/spiro/bridged/macrocycle facts without changing
Clean2D geometry or candidate behavior.

## Requirements

### Requirement: Topology evidence covers reusable graph structure

The topology summary SHALL expose derived ring systems, block adjacency,
connectors, flexible connectors, attachment atoms, articulation atoms and bonds,
branch points, fused/spiro/bridged evidence, and macrocycle blocks when present.

#### Scenario: Ring and block topology is explicit

- **GIVEN** a graph with rings and existing multilayer blocks
- **WHEN** topology evidence is collected
- **THEN** ring systems SHALL identify their ring and atom membership
- **AND** rigid/semi-rigid block IDs SHALL reference existing block records
- **AND** block adjacency and connector records SHALL expose stable endpoints.

### Requirement: Derived topology is general and deterministic

Topology facts SHALL be derived from `MolGraph`, `MultilayerChemicalGraph`, and
`BlockGraph` without molecule names, fixture IDs, new block kinds, or a parallel
molecular model.

#### Scenario: Independent topology families

- **GIVEN** acyclic, branched, monocyclic, fused, spiro, bridged, linker,
  multiblock, or macrocyclic graph structures
- **WHEN** each structure is summarized
- **THEN** the facts SHALL reflect graph connectivity and existing motifs/blocks
- **AND** repeated summaries SHALL have identical canonical ordering.

### Requirement: Topology output is strict JSON

The summary SHALL serialize with `json.dumps(summary, allow_nan=False,
sort_keys=True)`. Non-finite numeric values SHALL become `None`, and unordered
collections SHALL be canonicalized.

#### Scenario: Non-finite metadata is safe

- **GIVEN** block or connector metadata contains a non-finite float
- **WHEN** topology evidence is serialized
- **THEN** serialization SHALL succeed without emitting `NaN`, `Infinity`, or
  `-Infinity`.

### Requirement: Evidence collection remains observational

Collecting strengthened topology evidence SHALL NOT change coordinates,
complexity classification, preserve-only policy, local/global repair flags,
candidate sources, ranking, debug snapshot outcome, block unwrap, or local graph
cleaner behavior.

#### Scenario: Before and after behavior is unchanged

- **GIVEN** a graph and its existing Clean2D execution path
- **WHEN** topology evidence is collected before and after execution
- **THEN** graph identity, coordinates, policy fields, candidate sources, and
  result state SHALL remain equal.
