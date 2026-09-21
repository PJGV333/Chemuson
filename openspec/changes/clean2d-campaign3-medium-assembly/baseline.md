# Campaign 3 baseline

- Base commit: `043fe96 Archive strengthened Clean2D topology contract`
- Worktree: clean before this change
- Existing global baseline after Campaign 2 strengthening: `1553 passed, 20 skipped, 4 failed`.
- The four failures are the established candidate/stereo failures:
  - `test_generate_candidates_attempts_rdkit_for_cyclic_graphs`
  - `test_chiral_smiles_import_creates_wedge_or_hash`
  - `test_amino_acid_chiral_smiles_import_creates_wedge_or_hash`
  - `test_tetrandrine_import_preserves_visual_stereo`
- Existing focused medium observations:
  - `multiblock_biphenyl_like`: applied; before `needs_rebuild`, after `good`.
  - `multiblock_diphenyl_ether_like`: applied; before `needs_rebuild`, after `good`.
  - `selection_boundary_biphenyl_linker`: failed-controlled with `worse-quality`; selection-boundary behavior is retained.
  - `fused_aromatic_current_baseline`: applied; before `needs_polish`, after `good`.
- No production or architecture files were changed before this OpenSpec.
