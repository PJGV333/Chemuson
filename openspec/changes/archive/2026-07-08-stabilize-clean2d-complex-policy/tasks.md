## 1. Tests first

- [x] 1.1 Add or adjust a focused test for `test_tetrandrine_like_hierarchical_blocks_do_not_select_local_graph` that asserts policy classification and that `local_graph` is not selected.
- [x] 1.2 Add or adjust a focused test for `test_complex_engine_does_not_call_global_redraw_candidates` that asserts global redraw candidates are not called for high-risk complex quick clean.
- [x] 1.3 Add classification unit tests for hierarchical blocks, macrocycles, cyclophanes, internal cavities, and complex bridged structures.
- [x] 1.4 Add tests proving simple structures are not forced into preserve-only by the complex policy.
- [x] 1.5 Add diagnostic/snapshot tests that verify policy evidence is observable when available and that snapshots remain optional.

## 2. Complex-policy contract

- [x] 2.1 Identify the smallest policy data shape needed to represent preserve-only requirement, local-repair eligibility, global-redraw eligibility, internal-route eligibility, and stable reasons.
- [x] 2.2 Extend or adapt `classify_clean2d_complexity` to expose the required route decisions without generating candidates.
- [x] 2.3 Ensure hierarchical block-layout signals classify as high risk for global redraw.
- [x] 2.4 Ensure macrocycle, cyclophane, internal cavity, and complex bridged-system signals classify as high risk by default.
- [x] 2.5 Reuse existing stable reason vocabulary where possible; add only minimal metadata if a new policy detail is needed.

## 3. Engine routing integration

- [x] 3.1 Gate global redraw candidate generation using the complex-policy result before candidate generators are called.
- [x] 3.2 Gate `local_graph_cleaner` using explicit local-repair eligibility from the policy.
- [x] 3.3 Preserve existing safety gates for any route that remains allowed.
- [x] 3.4 Avoid changing general candidate ranking; only block or allow routes according to the policy contract.
- [x] 3.5 Stop and record a follow-up if satisfying the contract requires a broad `engine.py` redesign.

## 4. Diagnostics and snapshots

- [x] 4.1 Include JSON-serializable policy evidence in result or candidate metadata where available.
- [x] 4.2 Ensure `quality_diagnostic` can surface policy evidence without changing selection or ranking.
- [x] 4.3 Ensure opt-in debug snapshots capture available policy evidence when enabled.
- [x] 4.4 Verify disabled snapshots do not change complex-policy routing.

## 5. Out-of-scope guards

- [x] 5.1 Confirm this change does not attempt to fix `test_naphthalene_fused_system_does_not_collapse`.
- [x] 5.2 Confirm this change does not attempt to fix `test_smart_clean_repairs_distorted_structure_before_proposing`.
- [x] 5.3 Confirm no RDKit/CoordGen, MolGraph, canvas, or global layout-engine changes are introduced.

## 6. Verification

- [x] 6.1 Run the new complex-policy tests.
- [x] 6.2 Run relevant existing Clean 2D complex, local graph, quality reporting, and debug snapshot tests.
- [x] 6.3 Run `openspec validate stabilize-clean2d-complex-policy --strict`.
