# Corrective Campaign 3 baseline

- Reference pre-Campaign-3 commit: `043fe96` (`Archive strengthened Clean2D topology contract`).
- Comparison commit: `0b1647b` (`Archive Clean2D Campaign 3 medium assembly`).
- Comparison uses seven fixtures: biphenyl-like, diphenyl-ether-like, triphenyl-like, fused-plus-sidechain, ring-chain-ring, branched multiblock, and aromatic-plus-aliphatic branch.
- Baseline engine has no `global_block_placement` candidate and keeps triphenyl-like and branched multiblock in `preserve-only`.
- Current engine emits a topology-derived candidate for **7/7** fixtures; **3/7** are accepted by the candidate hard gates, and triphenyl-like plus branched multiblock select `global_block_placement`.
- Accepted candidates reduce visual score and introduce no new crossings; rejected candidates retain explicit reasons and remain visible in evidence.
- The complete deterministic record is `evidence/baseline.json`; it intentionally omits `runtime_ms`.
- The repository baseline still has the four documented pre-existing global failures; this corrective change does not alter their ownership.
