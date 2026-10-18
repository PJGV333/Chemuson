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
| Real GUI / hardware OpenGL acceptance | not run; remains pending per `manual-test.md` |
| Local tetrandrina identity | repaired and read back; InChI/formula/connectivity verified, stereo not inferred |
| User-library protection | full backup verified; only selected MolBlock/SMILES values changed; backup and config excluded from Git |

No visual success is inferred from mocks, offscreen widgets, mesh counts, `render_requested`, or OpenGL upload. The local template has no explicit stereochemistry; none is inferred. No push, merge, or PR was performed.
