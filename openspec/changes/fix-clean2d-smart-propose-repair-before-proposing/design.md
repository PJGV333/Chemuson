# Design: Smart Propose Repair Flow

## Implementation Details

### Analysis of Current Flow
In `rank_clean2d_candidates`:
1. `baseline_score` and `quality` are calculated for the current layout.
2. `baseline_bad` is set if `quality.quality_class == "needs_rebuild"`.
3. If `mode == PROPOSE` and `baseline_bad` is true, the engine immediately returns a "rebuild required" message.

### Proposed Changes
We need to introduce a "pre-propose repair" step. If `mode == PROPOSE` and the geometry is distorted:
1. Check if the layout is "distorted" (e.g., `quality_class == "needs_rebuild"`).
2. If it is, attempt a "safe" repair first. This repair should be a deterministic, non-generative operation (e.g., global normalization or a specific "repair" candidate).
3. If the repair is successful and moves the layout out of "needs_rebuild" into "needs_polish" or "good", then proceed to propose alternatives.
4. If no safe repair is possible, then return the "rebuild required" message.

### Technical Constraints
- Must preserve chemical invariants (atom/bond IDs, elements, charges, etc.).
- Must not modify fused/naphthalene behavior.
- Must be limited to the smart/propose flow.
- The repair must ensure bond lengths stay above `target * 0.70`.

### Impacted Files
- `src/chemuson/clean2d/engine.py`: Modify `rank_clean2d_candidates`.
