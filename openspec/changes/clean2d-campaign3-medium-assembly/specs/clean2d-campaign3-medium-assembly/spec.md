# Specification: Campaign 3 medium molecule assembly

## ADDED Requirements

### Requirement: Medium assembly plan is topology-derived and deterministic

The Clean2D medium assembly path SHALL expose a JSON-safe plan containing an anchor block, ordered block IDs, ordered connector IDs, and flexible connector count. The plan MUST be derived from `describe_clean2d_topology` and MUST be identical for identical graph topology and selected atom IDs.

#### Scenario: Multi-block medium graph receives a global assembly order

- **WHEN** a selected graph contains multiple rigid/semi-rigid blocks connected by attachment, bridge, or linker edges
- **THEN** the plan selects a deterministic anchor, traverses block adjacency, orders structural connectors before flexible linkers, and records the result in candidate metadata

### Requirement: Global assembly precedes local polish

The block-aware candidate SHALL consume the topology-derived connector order before any local polish candidate is considered. Existing hard gates, candidate ranking, and local polish safety checks MUST remain authoritative.

#### Scenario: Simple structures retain existing behavior

- **WHEN** a simple acyclic or single-ring graph does not require hierarchical block assembly
- **THEN** no medium assembly operation is applied and its result state and candidate behavior remain unchanged

### Requirement: Medium evidence preserves invariants and excludes runtime noise

Campaign evidence SHALL cover at least three multi-block medium/large controls and at least two simple controls. Each result MUST preserve graph identity, chemistry, stereo metadata, selection boundaries, finite coordinates, and JSON serialization. Observable comparisons MUST ignore runtime-only values and MUST make individual regressions visible.

#### Scenario: Assembly candidate preserves the graph contract

- **WHEN** the medium assembly candidate is generated or rejected
- **THEN** atom IDs, bond IDs/endpoints/order, element identity, charges, stereo metadata, connected components, and selection boundary metadata remain unchanged
