# Qt/OpenGL component probe

These artifacts preserve a real Qt OpenGL framebuffer run of `MolecularViewportWidget` for the repaired local tetrandrina fixture. The probe hosts the widget in a minimal `QMainWindow`/`QDockWidget`, then captures docked, floating, and re-docked states. It is **not** the ChemUSON main window and does not verify template loading or integrated desktop acceptance.

## Reproduction

Run from the repository root in an authorized graphical session with the project virtual environment and an isolated configuration directory:

```sh
CHEMUSON_CONFIG_HOME="$TMPDIR/chemuson-render-test/config-home" \
XDG_CONFIG_HOME="$TMPDIR/chemuson-render-test/config-home" \
./.venv/bin/python \
  openspec/changes/repair-tetrandrine-3d-editing-live-mm-navigation/evidence/component-opengl/component_probe.py
```

The probe writes its three PNGs and `tetrandrina_scene_gl.json` under `$TMPDIR/chemuson-render-test/baseline/`. The committed JSON and images are the captured run, not a synthetic rendering or a CPU mesh preview.

## Captured result

- Corrected fixture; RDKit conformer has 88 scene atoms and 94 bonds.
- Docked, floating, and re-docked framebuffer images each contain molecular pixels outside the lower-left XYZ-axis region; all three were visually inspected.
- Qt reports a current OpenGL context matching the renderer, OpenGL 4.6 Core / Mesa 26.2.4, and `gl_error=0` for all captures.
- The conformer is `unoptimized`, with no energy because MMFF94 did not converge.
- The probe's `QOpenGLBuffer.size()` observation is `0`. This value is preserved as observed; it is **not** used to claim an uploaded byte count. The independent framebuffer readback and visual inspection establish visible component pixels.

Integrated ChemUSON acceptance in the user's CachyOS desktop remains pending and is intentionally not claimed here.
