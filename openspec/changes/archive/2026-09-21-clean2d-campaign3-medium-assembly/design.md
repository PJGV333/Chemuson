# Design: Campaign 3 medium molecule assembly

## Assembly plan

`plan_clean2d_block_assembly` derives a JSON-safe plan from `describe_clean2d_topology`:

- choose the largest rigid/semi-rigid block as the deterministic anchor;
- traverse the block adjacency graph from that anchor;
- order attachment/bridge edges before flexible linker edges, then by connector ID;
- append disconnected blocks and edges deterministically;
- expose the ordered block IDs, connector IDs, and flexible connector count.

The existing `block_constraints` candidate consumes the connector order. Geometry operations remain the existing rigid-edge transforms; this change makes their global traversal explicit and stable rather than introducing a second block model.

## Safety

The plan is observational and does not mutate `MolGraph`, `BlockGraph`, coordinates, or policy state. Candidate metadata records the plan. Existing hard gates and ranking remain authoritative. `local polish` remains downstream of block assembly.

## Evidence

The campaign covers biphenyl-like, diphenyl-ether-like, and triphenyl-like medium/large multi-block cases, plus simple aromatic and acyclic controls. Evidence compares result state, candidate sources, topology plan, and geometry metrics with runtime excluded.
