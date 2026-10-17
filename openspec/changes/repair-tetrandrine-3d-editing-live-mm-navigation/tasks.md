# Tasks

## 0. OpenSpec gate
- [x] Trace the local tetrandrina template and capture the pre-change baseline.
- [x] Validate this change with `openspec validate repair-tetrandrine-3d-editing-live-mm-navigation --strict`.
- [x] Commit the validated specification separately before production code.

## 1. Tetrandrina calculation and scene publication — commit 2
- [x] Add a regression using the exact local template graph and assert the original graph is unchanged.
- [x] Expand supported abbreviations only on a calculation copy; preserve source stereo exactly.
- [x] Normalize only contradictory aromatic bridge flags with explicit single order and include diagnostics.
- [x] Publish coordinates for source atoms, derived abbreviation atoms, and H into a renderable scene.
- [x] Verify ETKDG, scene construction, finite mesh, and the nonconverged-MMFF status are separately reported.

## 2. Direct 3D drawing — commit 3
- [ ] Add RED tests for world-space preview geometry, snapping, empty-space creation, and Escape cancellation.
- [ ] Implement a 3D sphere-and-bond preview excluded from picking/document state.
- [ ] Keep existing-atom drops bond-only and empty drops atom+bond as one undo action.
- [ ] Implement and test the double-click rule with no duplicate first-click mutation.
- [ ] Verify valence rejection, explicit/derived H behavior, and undo/redo.
- [ ] Verify ethane and butane drawing/editing paths.

## 3. Live MM playback and telemetry — commit 4
- [ ] Record baseline worker-frame timing for ethane and butane (captured in `baseline.md`).
- [ ] Add tests for compatible snapshot interpolation, incompatible snapshots, and visual-only energy semantics.
- [ ] Add timed playback so grouped real frames remain visible progressively; never persist interpolated coordinates.
- [ ] Show force field, measured E, measured sample ΔE, step/max steps, and fixed-atom constraints with units.
- [ ] Clear pending frames and cancel stale work immediately on edits/document switches.
- [ ] Test cancellation, final-result publication, and ethane/butane real MM sequences.
- [ ] Measure playback tick cost and report workload, frame count, and timing separately from visual acceptance.

## 4. Camera navigation affordances — commit 5
- [ ] Add camera-basis axis projection tests across camera orientations.
- [ ] Render a persistent lower-left X/Y/Z indicator.
- [ ] Render transient orbit guidance only for empty-space camera drags and distinguish it from atom motion.
- [ ] Test that overlays do not affect picking or molecular geometry.

## 5. Final verification
- [ ] Run focused chemistry, interaction, animation, cancellation, viewport, and architecture tests.
- [ ] Run compileall, scoped Ruff, OpenSpec strict validation, and `git diff --check`.
- [ ] Do not rerun the full pytest suite while the known SIGSEGV remains open; record this limitation.
- [ ] If GUI remains unauthorized, keep `manual-test.md` and report visual acceptance as pending.
- [ ] Confirm all requested work is committed locally on `feature/chemuson-3d`, no push/merge/PR occurred, and Git status is clean.
