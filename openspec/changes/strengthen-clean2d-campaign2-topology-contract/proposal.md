# Proposal: Strengthen Clean2D Campaign 2 topology contract

## Why

Campaign 2 currently exposes generic block and connector records, but Campaign 3
needs reusable evidence for ring systems, articulation structure, branch points,
flexible connectors, and fused/spiro/bridged/macrocycle topology.

## What Changes

Extend `describe_clean2d_topology(...)` with deterministic derived evidence from
`MolGraph`, `MultilayerChemicalGraph`, and `BlockGraph`. Add an independent
nine-shape topology test matrix, disconnected/selected-subgraph coverage, strict
JSON checks, and an immutability test. No layout or candidate behavior changes.

## TARGET FAMILY

General medium/large graph topologies used by future block placement: acyclic,
branched, monocyclic, fused, spiro, bridged, linker-connected, multiblock, and
macrocyclic structures.

## NON-TARGET FAMILIES

Geometry, candidate generation, candidate ranking, local repair, block unwrap,
local graph cleaner, complexity policy decisions, GUI, persistence, and backend
routing.

## HARD INVARIANTS

- No second molecular graph or new `BlockKind` is introduced.
- Existing multilayer and block graph objects remain the source of block data.
- Evidence collection does not mutate graph data, coordinates, policies, or
  candidate behavior.
- JSON output passes `json.dumps(..., allow_nan=False, sort_keys=True)`.

## PROMOTION GATES

The focused matrix, existing topology/complexity tests, architecture tests,
compileall, Ruff, strict OpenSpec validation, full suite, and scope review must
pass with only the four documented baseline failures.

## ROLLBACK CONDITION

Revert if any new evidence is nondeterministic, non-JSON-safe, molecule-specific,
or changes coordinates, policy classification, candidate sources, ranking,
block unwrap, or local graph cleaning.
