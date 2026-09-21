# Design: Campaign 2 topology/decomposition

## Reuse existing layers

`build_multilayer_chemical_graph` remains the only decomposition entry point.
The implementation reads `MultilayerChemicalGraph.motif_graph` and
`BlockGraph`; it does not reconstruct rings, blocks, or connectors in a second
model. `Clean2DComplexityProfile` remains the policy-facing complexity surface.

## Observational summary

Add one JSON-safe summary adapter beside the existing complexity policy. The
adapter canonicalizes atom IDs, block IDs, motif IDs, block kinds, block edges,
and connected components using sorted values. It reports the existing
`BlockGraph` nodes and edges as decomposition evidence; it does not alter the
mutable `MolGraph`, coordinates, constraints, candidate generation, ranking, or
local repair.

The summary contains:

- selected atom IDs and covalent connected components;
- atom, bond, ring, and component counts;
- block records with IDs, kinds, atoms, anchors, motifs, and safe metadata;
- connector records derived from `BlockGraph.edges` with edge kind, block IDs,
  atom IDs, and weight;
- stable block-kind counts.

Only primitive JSON-safe values, lists, and mappings are emitted. Unknown
metadata values are represented by a stable string form rather than leaking
object identity or memory addresses.

## Determinism

The adapter sorts all records by stable IDs and all ID collections
lexicographically/numerically. Repeated calls over the same graph and layer
model must produce equal mappings and equal canonical JSON.

## Scope boundary

This campaign initially changes the reporting surface and tests only. It does
not route Clean2D through the summary, change `block_unwrap`, change
`local_graph_cleaner`, alter geometry, or add a dependency. A later campaign may
use the proven decomposition evidence for placement or routing under a separate
OpenSpec.
