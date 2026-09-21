# Design: Campaign 3 global block placement

## Candidate construction

Use `plan_clean2d_block_assembly` and the topology summary as a block tree/forest. Select the largest rigid or semi-rigid block as the root. For each parent-child attachment:

1. keep the child block's internal coordinates;
2. calculate the parent attachment direction from the parent centroid;
3. allocate separated angular sectors using sibling count and occupied sectors;
4. rotate the child around its attachment atom toward the assigned sector;
5. translate it to the target attachment distance;
6. continue breadth-first through the block graph.

The operation is topology-derived and contains no fixture-name branches. It is not force-directed and does not implement Campaign 5 routing.

## Safety and ranking

The candidate is generated before the complex `preserve-only` return, but it does not disable that policy. It must pass existing graph invariants, stereo/selection checks, crossing and collision checks, ring degeneracy, finite-coordinate validation, and displacement limits. Rejected candidates remain visible with a stable reason; if no candidate passes, the existing preserve-only result is returned.

## Evidence

Run the same regression corpus at baseline commit `043fe96` and the corrected HEAD. Record operation result, selected source, candidate metrics, and a separate internal assembly record for direct-connected rings, one-atom connectors, flexible chains, three blocks, branched multiblocks, ring-chain-ring, and aromatic/aliphatic branches. Runtime is excluded from semantic comparison.
