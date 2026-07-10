# Spec: Fix smart propose repair

## Description
Correct the behavior where `smart propose` rejects distorted structures without attempting a preliminary repair.

## Requirements
1. When `mode == PROPOSE` and the input geometry is distorted (quality class `needs_rebuild`), the engine should attempt a "safe" repair.
2. The "safe" repair must be a deterministic operation that preserves chemical invariants.
3. If the repair is successful (moving the layout to a better quality class), the engine should then proceed to propose alternatives.
4. If the repair is unsuccessful or not possible, the engine should return the standard "rebuild required" message.
5. Final bond lengths after repair must be >= 70% of the target length.

## Invariants
- Atom/Bond IDs and counts must be preserved.
- Chemical elements, charges, and stereo metadata must be preserved.
- Fused systems (naphthalene) must not be affected by this change.

## Success Criteria
- `test_smart_clean_repairs_distorted_structure_before_proposing` passes.
- Other existing tests pass.
- `openspec validate --all --strict` passes.
