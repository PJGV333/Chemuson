## Follow-ups

- Review conceptual duplication between `engine.py`, `safety.py`, `depiction_candidates.py`, and `Clean2DController` in a separate OpenSpec change if a broader migration becomes necessary.
- Consider whether imported depiction candidate diagnostics in `chemio/depiction_candidates.py` need their own stable adapter once a consumer requires that surface.
- Consider whether controller-level UI messaging should display stable diagnostic reasons in a future UI-focused change.
- Investigate unrelated Clean 2D geometry/complex-structure test failures observed while running the broad related suite: `test_naphthalene_fused_system_does_not_collapse`, `test_tetrandrine_like_hierarchical_blocks_do_not_select_local_graph`, `test_smart_clean_repairs_distorted_structure_before_proposing`, and `test_complex_engine_does_not_call_global_redraw_candidates`. These are outside this reporting-only change because fixing them would require algorithmic or geometry behavior changes.
