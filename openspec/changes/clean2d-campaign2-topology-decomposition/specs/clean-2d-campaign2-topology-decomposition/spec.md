# Specification: Clean2D Campaign 2 topology/decomposition

## ADDED Requirements

### Requirement: Topology decomposition reuses the multilayer model

The Campaign 2 topology summary SHALL derive blocks, motifs, connectors, and
metadata from `MultilayerChemicalGraph` and its existing `BlockGraph` rather
than creating a parallel molecular graph.

#### Scenario: Existing multilayer decomposition is summarized

- **GIVEN** a medium or large `MolGraph` and its multilayer model
- **WHEN** the topology summary is requested
- **THEN** it SHALL expose the existing block nodes and block edges
- **AND** it SHALL preserve the selected atom IDs and existing block kinds
- **AND** it SHALL not mutate the graph or coordinates.

### Requirement: Decomposition evidence is JSON-safe and stable

The topology summary SHALL contain only JSON-safe primitive values, lists, and
mappings, and SHALL canonicalize IDs and collections so repeated summaries of
the same input are equal.

#### Scenario: Repeated summary is deterministic

- **GIVEN** the same graph and multilayer model
- **WHEN** the summary is requested twice
- **THEN** the mappings and canonical JSON representations SHALL be equal
- **AND** block, connector, motif, and component order SHALL not depend on hash
  iteration order.

### Requirement: Blocks and connectors remain auditable

The summary SHALL expose atom count, bond count, ring count, connected
components, block count, block-kind counts, block membership, anchors, motif
IDs, connector count, connector kinds, connector endpoints, and connector
weights when available.

#### Scenario: Block and connector evidence is complete

- **GIVEN** a structure containing multiple rigid blocks joined by linkers or
  bridges
- **WHEN** the summary is serialized
- **THEN** each reported block SHALL identify its stable ID, kind, atom IDs,
  anchors, and motif IDs
- **AND** each reported connector SHALL identify its stable edge ID, kind,
  participating block IDs, atom IDs, and weight.

### Requirement: Campaign 2 is observational initially

The topology summary SHALL NOT alter Clean2D geometry, candidate generation,
ranking, local graph cleaning, block unwrap, backend routing, or multilayer
constraint behavior.

#### Scenario: Geometry remains outside the campaign

- **GIVEN** a graph and an existing Clean2D execution path
- **WHEN** topology evidence is collected
- **THEN** the existing coordinates, candidate sources, result state, and
  quality decisions SHALL be unchanged by collecting the summary.
