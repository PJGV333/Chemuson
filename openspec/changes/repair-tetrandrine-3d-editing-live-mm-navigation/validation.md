# Validation report

## Baseline and final focused suite

See `baseline.md` for the captured baseline. The historical full-suite Qt/SIP SIGSEGV remains open, so the full pytest suite was not rerun. The final scoped integrated run was:

```text
458 passed, 2 skipped in 191.72s
```

It covered local-template/geometry chemistry, service and RDKit live streaming, 3D dock/edit interactions, viewport, camera/geometry/picking, and architecture tests. Compileall passed. OpenSpec strict validation, scoped Ruff on changed Python files, and `git diff --check` are recorded after final documentation changes.

## Chemical identity audit

`identity-audit.md` records the local MolBlock/fixture match, composition, graph and aromaticity comparison. The source graph is malformed before import; import and abbreviation expansion are faithful. The generated 47-heavy-atom / 93-total scene sanitizes only after demoting bridge bond 9 but remains `C39H46N2O6`, not reference tetrandrina `C38H42N2O6`. Atom 11 duplicates atom 6's 2D coordinates, is terminal, and is attached to atom 5 while aromatic closure 5–6 is absent. No source template or fixture was changed. A graph-equivalent candidate was tested only on a disposable RDKit copy and was not applied.

The original conformer probe reported `method=unoptimized`, `energy=None`, `minimization_not_converged`; generation took 0.264745 s and produced a scene of 93 atoms / 98 bonds and a finite 120,384-vertex CPU mesh. This verifies pipeline output for the malformed graph only—not tetrandrina identity and not OpenGL pixels.

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
| `pytest` focused integrated run | 458 passed, 2 skipped (191.72 s) |
| `python -m compileall -q src tests tools packaging` | pass |
| Ruff on changed Python files | pass |
| `openspec validate repair-tetrandrine-3d-editing-live-mm-navigation --strict` | pass |
| `git diff --check` and staged diff check | pass |
| Full test suite | not rerun; historical Qt/SIP SIGSEGV remains open |
| Real GUI / hardware OpenGL acceptance | not run; manual steps in `manual-test.md` |
| Local tetrandrina identity | **blocked** pending explicit authorization to correct the user-local template; current scene is not the reference molecule |

No visual success is inferred from mocks, offscreen widgets, mesh counts, `render_requested`, or OpenGL upload. The local template has no explicit stereochemistry; none is inferred. No push, merge, or PR was performed.
