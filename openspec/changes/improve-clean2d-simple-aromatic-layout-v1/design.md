## Context

Clean 2D has regression corpus, metrics, deterministic baselines, baseline reports, diff review, and blocking-policy infrastructure. That makes it possible to attempt a small production change and inspect its effects across the corpus instead of relying on isolated snapshots.

The target family is simple isolated or substituted aromatic rings in the regression corpus. The implementation should use the smallest safe production hook discovered in the existing Clean 2D pipeline, favoring candidate ranking/acceptance over new geometry engines unless the code already exposes a narrowly reusable coordinate utility.

## Goals / Non-Goals

**Goals:**
- Produce at least one real, measurable Clean 2D improvement for a simple aromatic corpus case.
- Preserve atom/bond identity, endpoints, orders, aromaticity, stereo metadata, selection metadata, and document structure.
- Avoid regressions for high-risk complex-policy cases and protected tetrandrine-like coverage.
- Document before/after changes with baseline report compare/review output.

**Non-Goals:**
- Rewrite `engine.py` or introduce a new layout engine.
- Change UI, public formats, RDKit/CoordGen integrations, selection behavior, or chemical invariants.
- Repair macrocycles, cyclophanes, tetrandrine-like cases, or other high-risk complex structures.
- Convert geometry metrics or diff review results into gates.

## Decisions

- Start with production candidate selection/ranking or safety acceptance because it can select an already-generated safe aromatic candidate without changing chemical data structures.
- Restrict any aromatic-specific logic to simple, low-risk structures detected by existing graph/corpus characteristics where possible. Complex/high-risk cases must continue through the existing protected policy path.
- Use tests to compare pre-existing and improved behavior at the API/result/metric level rather than updating expected corpus states to hide differences.
- Use temporary baseline reports under `/tmp` for empirical before/after evidence and keep those JSON files out of git.

## Risks / Trade-offs

- A ranking adjustment could select a better aromatic candidate but affect unrelated molecules. Mitigation: constrain the scoring/acceptance change to simple aromatic candidates and run the requested regression suites.
- Safety relaxation could accidentally accept worse candidates. Mitigation: preserve existing invariant, crossing, collision, ring-degeneracy, and complex-policy checks.
- Baseline changes may be observationally significant. Mitigation: use compare/review and confirm no unexpected contract-risk outside the intended family.
