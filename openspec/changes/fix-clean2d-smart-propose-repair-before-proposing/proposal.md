# Proposal: Fix smart propose to repair distorted structures before proposing alternatives

## Context
The current implementation of `smart propose` fails when the input structure is heavily distorted. Instead of first attempting a "safe" repair to reach a baseline quality, it rejects the layout because the base quality is too low. This leads to the user seeing a "rebuild required" message and no proposed alternatives, even if a simple repair could have made the structure acceptable for proposal.

## Problem
The test `test_smart_clean_repairs_distorted_structure_before_proposing` reveals that when the structure is distorted (e.g., bond lengths < 21.47 instead of 40.0), `smart propose` rejects candidates immediately because the `baseline_bad` flag is set due to the poor quality of the initial layout.

## Proposed Solution
Modify the `rank_clean2d_candidates` flow to:
1. Detect if the geometry is "bad" but repairable.
2. Apply a preliminary repair (like a global normalization or local nudge) if the layout is too distorted to proceed with standard proposal.
3. Ensure that the resulting repaired geometry satisfies the minimum requirements (e.g., bond lengths >= 0.7 * target) before proceeding to propose alternatives.

This ensures that `smart propose` is helpful for distorted structures by first "cleaning" them to a workable state.
