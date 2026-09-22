# Design: Complete Clean2D Campaign 4 target improvement

## Topology and ownership

Reuse `MultilayerChemicalGraph`, `MotifGraph`, `BlockGraph`, `BlockKind`, `BlockEdgeKind`, `describe_clean2d_rigid_systems`, and the existing safety/evidence surfaces. Do not create a second topology representation.

The descriptor must expose `multiple_rigid_blocks` and `rigid_system_count` for linked rigid systems. A topology-built fixture with at least three fused rings must be classified as `polycyclic`, not inferred from a fixture name.

## Candidate semantics

`_candidate_from_rigid_multiring_layout` may construct only topology-derived local geometry:

- For spiro systems, keep the shared spiro center fixed and orient each ring subspace independently into separate sectors.
- For fused systems with substituents, adjust only the first exocyclic substituent atom when terminal and safe; never route a chain or translate the global block.
- For congested attachments, preserve the attachment bond length approximately while changing only a safe local direction.
- For bridged systems where a safe reconstruction is unavailable, emit a rejected candidate with `preserve-only` fallback.

The candidate must never change Campaign 3 global block placement. It may consume the existing assembly plan as read-only context.

## State model

Candidate construction records only construction-time facts:

- `hard_gate_checks`: deterministic JSON-safe mapping from gate name to `bool`.
- `hard_gates_passed`: conjunction of that mapping.
- metrics before and after construction.

Post-construction reporting records separately:

- `accepted_by_engine`: survived general candidate evaluation/ranking.
- `selected`: final engine selection.
- `source`, rejection reason, and deltas.

Construction must not set `accepted_by_engine` or `selected` equal to hard-gate status.

## Evidence

Corrective evidence is deterministic JSON, excludes runtime from canonical comparison, records baseline/current commit references, and reports every target's construction, acceptance, selection, source, metrics, delta, and rejection reason. It includes protected Campaign 3, acyclic, aromatic, and stereo-sensitive controls.
