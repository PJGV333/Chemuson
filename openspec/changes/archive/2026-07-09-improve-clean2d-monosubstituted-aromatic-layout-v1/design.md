## Context

The previous Clean 2D production change introduced `simple_aromatic_template` for isolated aromatic rings and intentionally scoped it to quick/publication modes. Monosubstituted aromatics such as toluene-like, phenol-like, and aniline-like corpus cases are the next smallest aromatic family because they contain one ring and one terminal substituent.

This change should reuse the same candidate-selection pipeline and safety checks. It should not alter propose mode, complex-policy routing, or broad candidate ranking.

## Goals / Non-Goals

**Goals:**
- Improve at least one monosubstituted aromatic corpus case with measurable geometry or visual metrics.
- Preserve chemical identity, aromaticity, stereo metadata, selection metadata, and existing safety checks.
- Keep the new candidate out of propose, fused systems, high-risk complex-policy, tetrandrine-like, and naphthalene/fused coverage.

**Non-Goals:**
- Support disubstituted or para-disubstituted aromatics.
- Support fused aromatics, macrocycles, cyclophanes, tetrandrine-like cases, or any complex-policy guard family.
- Change UI, public formats, RDKit/CoordGen integrations, geometry gates, or expected corpus states to hide regressions.

## Decisions

- Add a distinct `monosubstituted_aromatic_template` source instead of overloading `simple_aromatic_template`. This keeps baseline diffs and debug snapshots explicit.
- Limit detection to a single 5/6-member neutral aromatic ring with exactly one neutral terminal substituent atom for v1. A one-atom substituent is enough for the target corpus and avoids accidental chain/fused behavior.
- Place the substituent radially outward from the regularized ring center through the ring anchor. If the original substituent direction is already outward enough, the candidate still normalizes bond length while preserving orientation intent.
- Let existing invariant, safety, and ranking checks accept or reject the candidate; do not relax chemical or geometric safety.

## Risks / Trade-offs

- The candidate may duplicate existing `clean2d_v2` improvements. Mitigation: use explicit source/metadata and tests to prove measurable improvement and scope.
- Dedupe may hide equivalent candidates. Mitigation: allow the explicit monosubstituted candidate to survive dedupe like the simple aromatic template.
- Over-detection could affect protected families. Mitigation: terminal-substituent-only detection, no propose generation, and tests for fused/high-risk/tetrandrine exclusions.
