# Depth-occlusion framebuffer evidence

These are real `QOpenGLWidget.grabFramebuffer()` captures made in the current logged-in Wayland desktop with Qt/OpenGL. They are component evidence, not an integrated ChemUSON application acceptance.

## Reproduce

From the repository root, with the project's virtual environment and an authorized graphical session:

```bash
QT_QPA_PLATFORM=wayland ./.venv/bin/python openspec/changes/repair-tetrandrine-3d-editing-live-mm-navigation/evidence/depth-occlusion/component_depth_occlusion_probe.py before
QT_QPA_PLATFORM=wayland ./.venv/bin/python openspec/changes/repair-tetrandrine-3d-editing-live-mm-navigation/evidence/depth-occlusion/component_depth_occlusion_probe.py after
```

`before` dynamically loads the `MolecularViewportWidget` and mesh implementation from baseline revision `74ebc9a`; it does not alter the worktree. `after` uses the current implementation. Both runs use the same deterministic scenes: ethane and benzene from fixed SMILES/conformer seed 29, and the unchanged repaired tetrandrina fixture with seed 173. The tetrandrina scene is reported as unoptimized where minimization did not converge; no chemistry or template was changed for these images.

Each run captures all three molecules in four states: docked, after camera orbit, floating, and re-docked. The adjacent `before/results.json` and `after/results.json` record scene composition, camera state, context depth-buffer size, and depth-test/write-mask/function values at molecular draw entry. The Qt surface reported a 24-bit depth buffer; `GL_DEPTH_FUNC` was `GL_LESS` (`513`), and `GL_DEPTH_WRITEMASK` was true.

The before captures show atom/bond overlap from lost depth testing after `QPainter`; draw-time samples enter with `GL_DEPTH_TEST` false. After the fix, draw-time samples show depth testing enabled across frames, including after orbit and dock recreation. Visual inspection of the real before/after framebuffers confirms the foreground spheres occlude rear bonds/atoms and ball-and-stick rods meet the sphere surfaces. `GL error=0` is recorded but is not used alone as evidence of visual correctness.

## Files

- `captures/before/`: 12 pre-fix framebuffer PNGs and machine-readable manifest.
- `captures/after/`: 12 post-fix framebuffer PNGs and machine-readable manifest.
- `captures/axial/before/`: real pre-fix bond-through-sphere and sphere draw-order control captures.
- `captures/axial/after/`: real post-fix ball-and-stick, space-filling, rotated and re-docked pixel-regression captures. Recreate them with:

  ```bash
  QT_QPA_PLATFORM=wayland CHEMUSON_DEPTH_EVIDENCE_DIR=openspec/changes/repair-tetrandrine-3d-editing-live-mm-navigation/evidence/depth-occlusion/captures/axial/after ./.venv/bin/python -m pytest -q tests/test_viewport3d_depth_occlusion.py::test_real_viewport_occludes_bonds_with_depth_after_qpainter_overlays
  ```
- `component_depth_occlusion_probe.py`: reproducible component capture and GL-state probe.
- `tests/test_viewport3d_depth_occlusion.py::test_real_viewport_occludes_bonds_with_depth_after_qpainter_overlays`: real pixel/state regression, including the known axial overlap in ball-and-stick and space-filling.

The full ChemUSON main-window workflow, manual picking/editing acceptance, and user-run CachyOS integrated acceptance remain outside this component probe.
