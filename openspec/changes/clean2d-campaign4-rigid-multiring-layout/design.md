# Design: Clean2D Campaign 4 rigid and multiring layout

## Existing infrastructure

- `MultilayerChemicalGraph`, `MotifGraph`, `BlockGraph`, `BlockKind`, `BlockEdgeKind`, and `build_multilayer_chemical_graph` already represent rings, fused blocks, linkers, terminal substituents, bridges, and cavities.
- `describe_clean2d_topology` already emits deterministic ring-system, block, connector, articulation, and attachment evidence.
- `_cycle_basis_ordered`, `ring_degeneracy_score`, `classify_clean2d_layout_quality`, `evaluate_clean2d_layout`, `stereo_layout_signature`, and `assert_clean2d_invariants` are existing geometry and safety surfaces.
- `fused_aromatic_template`, `scaffold_depiction`, `block_unwrap`, and Campaign 3 global placement are existing candidate sources with established fallback behavior.

## Reuse and extension

Extend `complex_policy.py` with a general rigid-system descriptor derived from the existing multilayer model and graph coordinates. For each canonical ring system it records member atoms/bonds, ring membership, shared atoms/bonds, external attachments/neighbors, centroid, principal orientation, attachment vectors, external substituent count, and local congestion. The descriptor is observational, JSON-safe, sorted, and does not create a second topological representation.

Add a bounded `rigid_multiring_layout` candidate in `engine.py`. It may rotate a rigid system around its existing centroid and choose immediate outward attachment directions. It SHALL NOT translate complete blocks, change Campaign 3 parent-child/root decisions, route flexible chains, or enumerate an unbounded global search. Candidate selection uses existing quality metrics and hard gates; local orientation cost is only an internal deterministic construction heuristic, not Campaign 7 global ranking.

The candidate is emitted for reusable topology signals: fused, spiro, bridged, polycyclic, congested attachment, or multiple rigid systems. Systems with explicit stereo or unsafe bridge geometry may be rejected and fall back to existing safe candidates/preserve-only.

## Candidate safety

Preserve atom/bond identity, endpoints, bond order, aromaticity, charges, stereo signature, finite coordinates, ring degeneracy, no-new-crossing, collision, bond-length, bounding-box, selection, and displacement contracts. The candidate's local displacement budget is finite and derived from target bond length plus rigid-system diameter; it is recorded with all hard-gate results. A safe candidate must improve an applicable ring/attachment/visual metric or remain rejected as no-improvement.

## Evidence

Generate a deterministic `evidence/baseline.json` from commit `346b229` and `evidence/current.json` from the implementation. Each of the twelve required fixture families records topology, baseline/current result and source, candidate metadata, metrics, hard gates, delta, and whether the internal candidate was selected. Runtime is excluded from deterministic comparison.

## Campaign 5 boundary

Flexible connector torsion, long-branch routing, branch competition, and global redistribution of child blocks remain unchanged and are reported as non-goals. Campaign 4 only evaluates the first bond/attachment direction adjacent to the rigid system.
