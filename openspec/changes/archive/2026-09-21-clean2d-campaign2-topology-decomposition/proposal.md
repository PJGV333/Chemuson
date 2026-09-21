# Proposal: Clean2D Campaign 2 — Topology and Decomposition

## Why

Campaign 1 established reproducible evidence, but the existing multilayer
decomposition was not exposed as a stable, auditable summary for medium and
large structures.

## What Changes

Add an observational JSON-safe summary over the existing multilayer and block
graph objects, with deterministic block, connector, component, and topology
metadata. Geometry and candidate selection remain unchanged.

## TARGET FAMILY

Medium and large structures whose multilayer model contains rings, fused systems,
rigid/semi-rigid blocks, linkers, bridges, or terminal substituents.

## NON-TARGET FAMILIES

Simple acyclic molecules, local geometric polishing, global block placement,
connector routing, candidate ranking, backend selection, and GUI behavior.

## EXPECTED IMPROVEMENT

Make the existing multilayer decomposition auditable and deterministic by
exposing stable block, connector, component, ring, and topology metadata without
changing coordinates or Clean2D candidate selection.

## HARD INVARIANTS

- `MolGraph` remains the source of truth for atom and covalent-bond identity.
- Atom IDs, bond IDs, endpoints, orders, aromaticity, stereo metadata, and
  connected components remain unchanged.
- The summary is observational and JSON-serializable.
- Existing `MultilayerChemicalGraph`, `BlockGraph`, motifs,
  `Clean2DComplexityProfile`, `local_graph_cleaner`, `block_unwrap`, and
  multilayer constraint contracts are reused.
- No new parallel molecular graph or geometry algorithm is introduced.

## METRICS OBSERVED

- atom, bond, ring, and connected-component counts;
- deterministic block count and block-kind counts;
- deterministic connector/edge count and kinds;
- block atom IDs, anchors, motif IDs, and connector endpoints;
- stable serialization across repeated decomposition;
- existing Campaign 1 before/after evidence remains unchanged.

## BASELINE INPUT

- Campaign 1 archived evidence:
  `openspec/changes/archive/2026-09-21-clean2d-campaign1-benchmark-observability/evidence/baseline.json`
- Campaign 2 baseline command results are recorded in `baseline.md`.
- Existing decomposition tests and multilayer constraints are the behavioral
  baseline; no production geometry changes are permitted in this phase.

## PROMOTION GATES

1. The Campaign 2 OpenSpec validates strictly.
2. New topology contract tests pass.
3. Existing multilayer, block unwrap, complex policy, architecture, and full
   regression tests preserve their baseline outcomes.
4. The summary is deterministic, JSON-safe, and derived only from existing
   multilayer structures.
5. No changes occur in GUI, persistence, candidate ranking, or layout geometry
   behavior.

## ROLLBACK CONDITION

Revert the Campaign 2 commit if any existing geometry result, candidate source,
quality state, graph invariant, or architecture boundary changes unexpectedly,
or if topology summaries differ for identical graph/layer inputs.
