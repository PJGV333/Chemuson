# Design: Strengthen Campaign 2 topology contract

## Derived evidence

Keep `MultilayerChemicalGraph` and `BlockGraph` as the only block model. Extend
the existing summary with:

- ring systems grouped from existing ring motifs by shared atoms;
- rigid and semi-rigid block IDs classified from existing `BlockKind` values;
- block adjacency and connector records from existing block edges;
- flexible connector records, attachment atom IDs, and rotatable bond IDs;
- articulation atom and bond IDs derived from the selected covalent MolGraph;
- branch-point atom IDs derived from covalent degree;
- fused, spiro, bridged, and macrocycle evidence linked to motif/block IDs.

No new block kind is needed. Spiro and articulation facts are observational
properties derived from graph structure and existing motifs/blocks.

## Canonicalization and JSON

All unordered collections are sorted by numeric or stable tuple keys. Floats are
normalized through one JSON-safe helper; non-finite values become `None`.
Connector weights and nested block metadata use the same sanitizer. Tests use
`json.dumps(summary, allow_nan=False, sort_keys=True)`.

## Matrix and protection

Use independent synthetic graph builders with no molecule-name routing. Cover
acyclic chain, branched tree, monocycle, fused bicyclic, spiro, bridged,
two-rings-plus-flexible-linker, multiblock, and macrocycle. Add disconnected and
selected-subgraph cases. Capture graph coordinates, policy fields, and candidate
source/state evidence before and after summary collection to prove the adapter
is observational.
