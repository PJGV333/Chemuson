# Proposal: Tighten Campaign 3 global placement safety

## Problem

Campaign 3 added `global_block_placement`, but its displacement exception can bypass the local safety helper without an explicit finite global budget, and the complex-preserve path can return global placement before comparing it with safe scaffold or unwrap alternatives.

## Scope

- Add an explicit topology- and target-dependent finite displacement budget for global placement.
- Evaluate all global-placement hard gates explicitly after candidate construction.
- Compare safe candidates in the complex-preserve path using existing quality fields and deterministic source tie-breaking.
- Extend the seven Campaign 3 fixture records and tests with safety and competition evidence.

## Out of scope

- No changes to `safety.py` global policy.
- No Campaign 7 ranking system.
- No changes to historical baselines or the four known global failures.
- No Campaign 4 work.
