# Proposal: Campaign 3 medium molecule assembly

## Why

Medium structures with multiple rigid blocks and flexible connectors need an explicit global assembly order before local polishing. The current engine has block-aware operations, but their traversal order is implicit and is not exposed as reproducible evidence.

## What Changes

- Add a deterministic, topology-derived medium assembly plan with an anchor block, block order, connector order, and flexible connector count.
- Consume that plan in the existing block-constraint candidate before local polishing.
- Add medium multi-block regression evidence and a distributive comparison that requires no regression in simple cases.
- Preserve all chemical, stereo, selection, and graph invariants.

## Scope

In scope: `src/chemuson/clean2d/complex_policy.py`, the existing Clean2D block candidate in `engine.py`, regression tests and OpenSpec evidence.

Out of scope: new block models, changes to `core/layers.py`, RDKit policy, local polish heuristics, GUI behavior, or the existing Campaign 1/2 baselines.
