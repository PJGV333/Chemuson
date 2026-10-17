# Tasks

## 0. OpenSpec gate
- [x] Trace the local tetrandrina template and capture the pre-change baseline.
- [x] Validate this change with `openspec validate repair-tetrandrine-3d-editing-live-mm-navigation --strict`.
- [x] Commit the validated specification separately before production code.

## 1. Tetrandrina calculation and scene publication — commit 2
- [x] Add a regression using the exact local template graph and assert the original graph is unchanged.
- [x] Expand supported abbreviations only on a calculation copy; preserve source stereo exactly.
- [x] Normalize only contradictory aromatic bridge flags with explicit single order and include diagnostics; treat this as RDKit sanitization only, not identity repair.
- [x] Publish coordinates for the imported graph, derived abbreviation atoms, and H into a scene; do not equate scene publication with correct chemical identity.
- [x] Verify ETKDG, scene construction, finite mesh, and the nonconverged-MMFF status are separately reported.
- [x] Audit formula, composition, components, connectivity, expansion, aromaticity, and valences against the reference; record the local-template defect in `identity-audit.md` without changing the user's model.
- [ ] Apply a user-authorized correction to the local template and revalidate its identity; the current 93-atom scene is not tetrandrina acceptance.

## 2. Direct 3D drawing — commit 3
- [x] Add RED tests for world-space preview geometry, snapping, empty-space creation, and Escape cancellation.
- [x] Implement a 3D sphere-and-bond preview excluded from picking/document state.
- [x] Keep existing-atom drops bond-only and empty drops atom+bond as one undo action.
- [x] Implement and test the double-click rule with no duplicate first-click mutation.
- [x] Verify valence rejection, explicit/derived H behavior, and undo/redo.
- [x] Verify ethane and butane drawing/editing paths.

## 3. Live MM playback and telemetry — commit 4
- [x] Record baseline worker-frame timing for ethane and butane (captured in `baseline.md`).
- [x] Add tests for compatible snapshot interpolation, incompatible snapshots, and visual-only energy semantics.
- [x] Add timed playback so grouped real frames remain visible progressively; never persist interpolated coordinates.
- [x] Show force field, measured E, measured sample ΔE, step/max steps, and fixed-atom constraints with units.
- [x] Clear pending frames and cancel stale work immediately on edits/document switches.
- [x] Test cancellation, final-result publication, and ethane/butane real MM sequences.
- [x] Measure playback tick cost and report workload, frame count, and timing separately from visual acceptance.

## 4. Camera navigation affordances — commit 5
- [x] Add camera-basis axis projection tests across camera orientations.
- [x] Render a persistent lower-left X/Y/Z indicator.
- [x] Render transient orbit guidance only for empty-space camera drags and distinguish it from atom motion.
- [x] Test that overlays do not affect picking or molecular geometry.

## 5. Final verification
- [x] Run focused chemistry, interaction, animation, cancellation, viewport, and architecture tests.
- [x] Run compileall, scoped Ruff, OpenSpec strict validation, and `git diff --check`.
- [x] Do not rerun the full pytest suite while the known SIGSEGV remains open; record this limitation.
- [x] GUI/OpenGL acceptance was not run; keep `manual-test.md` and report visual acceptance as pending.
- [ ] Obtain explicit authorization before changing the defective user-local tetrandrina template; identity acceptance remains blocked.
- [x] Confirm all implemented work is committed locally on `feature/chemuson-3d`, no push/merge/PR occurred, and Git status is clean after the final documentation commit.
