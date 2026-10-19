# Validation report

## Baseline and focused suites

The earlier integrated viewport/interaction suite completed with 458 passed, 2 skipped in 191.72 s before the local-template repair. After the chemistry change, the final focused regression run was:

```text
48 passed, 1 skipped in 7.55s
```

This run covered aromatic bridge refusal across conformer/optimization/streaming paths, the preserved malformed fixture, corrected tetrandrina identity, RDKit generation/scene output, geometry service behavior, and the local Clean2D tetrandrina regression. The historical full-suite Qt/SIP SIGSEGV remains open, so the full pytest suite was not rerun. Compileall passed; scoped Ruff, strict OpenSpec validation, and `git diff --check` were run after final documentation changes.

## Chemical identity audit

`identity-audit.md` preserves the original defect evidence and records the authorized repair. The saved entry was reread through `graph_from_template_payload` without loading the write-on-load library: the local graph has 42 atoms/48 bonds and the independent corrected fixture matches its MolBlock. Supported abbreviation expansion yields 46 heavy atoms/52 bonds; the expanded formula is `C38H42N2O6`, with one component, cycle rank 7, 24 aromatic atoms/bonds, and no specified stereochemistry. This matches PubChem CID 73078 connectivity without stereo comparison. The stored SMILES parses to the same formula. A verified full backup was created before atomic JSON replacement; a byte-level comparison confirmed only the selected entry's MolBlock and SMILES string values changed; parsed library data outside those values is identical.

The repaired worker probe used seed 173, MMFF94, and 80 maximum iterations. ETKDG completed with 46 finite heavy-atom coordinates; MMFF94 did not converge (`method=unoptimized`, `energy=None`, `minimization_not_converged`). Scene construction yielded 88 atoms and 94 bonds plus a finite CPU mesh. These results verify chemical/CPU pipeline output only—not minimized geometry or OpenGL pixels. The old 93-atom output and 120,384-vertex CPU mesh remain historical evidence for the malformed graph and are not evidence for the repaired identity or visible rendering.

## Live MM sequence and telemetry

The real RDKit `optimize_stream` run used MMFF94, seed 29, `max_iters=240`, `steps_per_update=8`, and one fixed atom-position constraint (`fixed_atom_ids=(1,)`). Energies below are RDKit snapshot values in kcal/mol. The player interpolates coordinates for 80 ms per transition; interpolation carries no scientific energy and no gradient is calculated.

| Molecule | Worker time | Frames / steps | Snapshot energies E in emitted order (kcal/mol) | Gap min / median / max | Gaps <16 ms |
|---|---:|---|---|---:|---:|
| Ethane | 0.153254 s | 8 / 8–64 by 8 | `[-0.839715, -4.381175, -4.722530, -4.733974, -4.734346, -4.734365, -4.734365, -4.734365]` | 0.026 / 0.039 / 0.363 ms | 7/7 |
| Butane | 0.155223 s | 30 / 8–240 by 8 | `[1555.874633, 260.330341, 229.826548, 99.136901, 39.432425, 5.342332, -1.830166, -3.005857, -3.957392, -4.326724, -4.624488, -4.900675, -4.965823, -5.019794, -5.034051, -5.062962, -5.067789, -5.073646, -5.075038, -5.075363, -5.075562, -5.075839, -5.075884, -5.075928, -5.075944, -5.075955, -5.075961, -5.075965, -5.075969, -5.075970]` | 0.119 / 0.180 / 0.249 ms | 29/29 |

With an 80 ms visual transition and deterministic 60 Hz player sampling, ethane used 42 ticks (0.700 s) and butane 174 ticks (2.900 s), with zero dropped frames. CPU `scene_with_coordinates` + `build_scene_mesh` cost was 4.973 ms mean / 5.223 ms p95 for ethane (10,104 vertices) and 8.046 ms / 8.735 ms for butane (17,736 vertices). Butane reached 240 steps with `converged=False`, final E `-5.075970386544 kcal/mol`; the UI must not imply convergence. These CPU timings do not include OpenGL presentation or establish visual FPS.

The final focused tests verify frame order, visual-only interpolation, sample energy/HUD semantics, cancellation on edits and stale-frame rejection, final snapshot publication, draw/snap/double-click/undo/redo/Escape, orientation projections, orbit-vs-atom movement, and overlay/picking isolation. They do not verify desktop pixels.

## Final checks and limitations

| Check | Result |
|---|---|
| Earlier integrated viewport/interaction suite (before local repair) | 458 passed, 2 skipped (191.72 s) |
| Post-repair focused chemistry suite | 48 passed, 1 skipped (7.55 s) |
| `./.venv/bin/python -m compileall -q src tests tools packaging` | pass |
| Scoped Ruff on `rdkit_backend.py` and `test_geometry3d_aromaticity.py` (`F401,F811,F821,E722,E741`) | pass |
| Broad Ruff selection over `src tests tools packaging` | exit 1 only for pre-existing `F401` unused `math` at `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3` |
| `openspec validate repair-tetrandrine-3d-editing-live-mm-navigation --strict` | pass |
| `git diff --check` | pass after final documentation edits |
| Full test suite | not rerun; historical Qt/SIP SIGSEGV remains open |
| Qt OpenGL component framebuffer (docked/floating/re-docked) | pass; see `evidence/component-opengl/`; not an integrated ChemUSON acceptance |
| Integrated ChemUSON GUI/OpenGL acceptance on CachyOS | pending user-run manual acceptance; no integrated GUI was launched |
| Local tetrandrina identity | repaired and read back; InChI/formula/connectivity verified, stereo not inferred |
| User-library protection | full backup verified; only selected MolBlock/SMILES values changed; backup and config excluded from Git |

Component pixel success is based on actual `QOpenGLWidget.grabFramebuffer()` readback plus visual inspection, not mocks, scene-mesh counts, `render_requested`, or upload status. It verifies the viewport component's dock lifecycle only. Integrated ChemUSON application acceptance remains pending. The local template has no explicit stereochemistry; none is inferred. No push, merge, or PR was performed.

## Compatibility closure — expanded coordinate IDs (2026-10-18)

`_geometry3d_position_ids_compatible` now requires every source MolGraph ID. Additional position IDs must match the exact declared calculation expansions, and `template_normalization.expanded_atom_ids` must agree with non-hydrogen `generated_atoms` entries named `calc:<id>` when both sources are present. Missing source positions, missing declared expansion positions, undeclared extras, duplicate declarations, and mismatched metadata are rejected. The repaired tetrandrina calculation declares `43–46` beyond the 42 authored graph IDs; those coordinates remain compatible and publish the 88-atom/94-bond scene.

The focused lifecycle regression verifies publication, geometry-command undo/redo, switching away from and back to the document tab, viewport reopening without a redundant generation job, and a CMSN save/load round trip. CMSN schema v1 persists only authored graph positions and its existing metadata whitelist. The reloaded graph-only projection remains compatible and is not falsely marked stale; the expanded calculation-only scene metadata is not claimed to round-trip.

## Component OpenGL framebuffer evidence

`evidence/component-opengl/` contains the reproducer, machine-readable JSON, and docked/floating/re-docked framebuffer captures. The probe used a real PyQt6 `MolecularViewportWidget` inside a `QMainWindow`/`QDockWidget` on Wayland, generated from the corrected local tetrandrina fixture. It is a component harness, not the ChemUSON main window.

| State | Scene | GL context / version | GL error | Framebuffer molecular pixels outside XYZ ROI |
|---|---:|---|---:|---:|
| Docked | 88 atoms / 94 bonds; 114,096 mesh vertices | current renderer context; OpenGL 4.6 Core, Mesa 26.2.4 | 0 | 43,548 |
| Floating | 88 atoms / 94 bonds; 114,096 mesh vertices | current renderer context; OpenGL 4.6 Core, Mesa 26.2.4 | 0 | 41,691 |
| Re-docked | 88 atoms / 94 bonds; 114,096 mesh vertices | current renderer context; OpenGL 4.6 Core, Mesa 26.2.4 | 0 | 41,907 |

Visual inspection confirms the molecular structure is present apart from the lower-left XYZ overlay in all three captures. The conformer remains `unoptimized` (`energy=None`, `minimization_not_converged`); this is a rendering check, not a minimization claim. The probe's `QOpenGLBuffer.size()` observation was `0`, so this report does not claim a verified VBO byte size. The framebuffer pixels and context/error observations establish component rendering only; they do not establish integrated application acceptance, fit/orbit behavior, or the user's desktop configuration.

## Closeout checks

- Baseline at `282cd1fdfc157f992030c1e377eabd81e04d3712`: compileall passed; collection passed; four focused compatibility/tab/undo/regen tests passed; scoped Ruff passed.
- RED check before the compatibility fix: missing declared expansion position failed as expected (`1 failed, 3 passed`).
- Final focused suites: `tests/test_compchem3d_dock.py`, `tests/test_geometry3d_persistence.py`, and `tests/test_persistence_contract.py` — `69 passed, 1 skipped in 52.33s`.
- `./.venv/bin/python -m compileall -q src tests tools packaging` — exit 0, no output.
- Scoped Ruff on `main_window.py`, `test_compchem3d_dock.py`, and the preserved component probe with `F401,F811,F821,E722,E741` — `All checks passed!`.
- `openspec validate repair-tetrandrine-3d-editing-live-mm-navigation --strict` — valid.
- `git diff --check` — pass.
- The preserved `component_probe.py` was executed from the repository; process exit 0 and all three real framebuffer captures were byte-verified against their saved evidence copies.
- Full `pytest -q` was not rerun, consistent with the documented Qt/SIP SIGSEGV constraint and the user's focused scope.

## Depth occlusion and bond surface integration (2026-10-18)

### Reproduction and diagnosis
- Used a real Qt/OpenGL component on Wayland, not a CPU-only mesh or a `GL error=0` proxy.
- The pre-fix framebuffer probe loaded `widget.py` and `mesh.py` from baseline revision `74ebc9a`; current code was not modified to produce the baseline.

- On each baseline frame, `GL_DEPTH_TEST` was false before molecular drawing; `GL_DEPTH_WRITEMASK` was true, `GL_DEPTH_FUNC` was `GL_LESS` (513), and the Qt context reported a 24-bit depth buffer. This isolates the lost enable state rather than a missing buffer or write mask.
- The controlled axial framebuffer visibly shows the center-to-center bond overwriting the foreground sphere. The no-bond control shows the rear sphere overwriting the front sphere. Both images are under `evidence/depth-occlusion/captures/axial/before/`.

### Implementation and regressions
- `MolecularViewportWidget.paintGL()` restores depth testing, write mask, `GL_LESS`, and clear-depth 1.0 before clearing both color and depth on every frame.
- Ball-and-stick cylinder lanes are shortened individually to meet each endpoint sphere using the lane offset and bond radius; tests cover single, double, triple, and aromatic lanes. `space_filling` retains its spheres and is verified by real framebuffer pixels.
- Added fake-state and real Qt/OpenGL regressions. The real test runs repeated frames with orientation axes, Live MM HUD, orbit guide, and geometry handle; it checks state at the molecular draw boundary, front-atom pixel color, rotation, and dock floating/re-docking.

### Final validation
- Focused functional suites (viewport depth, widget, camera, geometry, picking, geometry service, Live MM stream, CompChem dock, and persistence): `197 passed, 3 skipped in 74.49s`.
- Real Wayland framebuffer regression with axial occlusion, overlays, camera rotation, and dock cycle: `1 passed in 2.21s`; optional image output was saved under `evidence/depth-occlusion/captures/axial/after/`.
- `./.venv/bin/python -m compileall -q src tests tools packaging` — exit 0.
- Scoped Ruff (`F401,F811,F821,E722,E741`) on modified Python sources/tests and the evidence probe — `All checks passed!`.
- `./.venv/bin/python -m pytest --collect-only -q` — `2281 tests collected in 0.57s`.
- `openspec validate repair-tetrandrine-3d-editing-live-mm-navigation --strict` — valid after all documentation updates.
- Real component captures: 12 before and 12 after across ethane, benzene, and tetrandrina (docked, rotated, floating, re-docked). All 36 pre-fix draw samples had depth testing off; all 36 post-fix samples had it on. Both runs reported write mask true, `GL_LESS`, 24 depth bits, and `GL error=0`.
- Coordinate SHA-256 values, molecular identities, and energies match exactly between before/after runs; tetrandrina was rendered from the unchanged repaired fixture with `method=unoptimized`, without changing chemistry, ETKDG, or its template.
- Visual inspection of the axial regression and ethane, benzene, and tetrandrina captures confirms opaque atom spheres and no bond overdraw through front surfaces. `GL error=0` alone is not treated as acceptance evidence.
- Full `pytest -q` remains unrun due the documented Qt/SIP SIGSEGV. The integrated ChemUSON GUI was not launched; user manual acceptance remains pending.
