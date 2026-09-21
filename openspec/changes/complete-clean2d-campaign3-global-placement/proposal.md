# Proposal: Complete Campaign 3 global block placement

## Why

Campaign 3 currently exposes deterministic connector traversal, but its evidence does not demonstrate that the internal assembly candidate improves medium multiblock layouts. Complex profiles can enter `preserve-only` before a topology-aware candidate can compete, while external candidates may already produce a good result.

## What Changes

- Add a topology-derived global block placement candidate before the complex preserve-only exit.
- Preserve each block's internal geometry while placing children from parent attachments, free angular sectors, sibling occupancy, and target separation.
- Apply existing invariants and hard gates without weakening preserve-only protection.
- Capture explicit baseline-vs-current and internal-candidate metrics for a broader medium family.

## Scope

In scope: the existing Clean2D engine and Campaign 3 evidence/tests, plus this corrective OpenSpec and the historical Campaign 2 promotion checkboxes.

Out of scope: Campaign 4, detailed connector routing, force-directed optimization, RDKit changes, GUI behavior, and changes to `core/layers.py`.
