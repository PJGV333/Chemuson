# Validation report

## Baseline
See `baseline.md`. Baseline focused tests passed (164 passed, 2 skipped); compileall and collection passed. Repository-wide Ruff remains red for one pre-existing F401. Full pytest remains unrun because of the known Qt/SIP SIGSEGV.

## Per-commit validation

- Exact local template: `method=unoptimized`, `energy=None`, diagnostic `minimization_not_converged`; generation 0.264745 s, scene 93 atoms / 98 bonds, CPU mesh 120,384 vertices, all finite. This verifies coordinate/scene/mesh construction only, not OpenGL pixels.
- Template graph fingerprint by structure: source 43 atoms / 48 bonds; calculation copy 47 / 52; source graph unchanged; aromatic bridge bond 9 normalized only on the calculation copy; no stereo specified or inferred.

| Commit track | Chemistry / interaction / animation tests | Architecture | Compileall | Scoped Ruff | OpenSpec / diff check |
|---|---|---|---|---|---|
| OpenSpec gate | baseline captured | pending | baseline pass | baseline pass scoped | `openspec validate ... --strict` valid; `git diff --check` pass |
| Local template | 340 passed, 1 skipped (chemistry + architecture group) | pass | pass | pass (Python scoped) | OpenSpec strict valid; `git diff --check` pass |
| Direct 3D drawing | pending | pending | pending | pending | pending |
| Live MM playback | pending | pending | pending | pending | pending |
| Camera navigation | pending | pending | pending | pending | pending |

## Acceptance caveats

- No visual success may be inferred from a mock, offscreen widget, mesh count, `render_requested`, or a successful OpenGL upload alone.
- The user-local tetrandrina template has no explicit stereo. The implementation must not assign stereochemistry; nonconverged minimization must remain reported as such.
- The full pytest suite is not rerun while the historical SIGSEGV remains open.
- Real GUI/OpenGL validation is pending unless explicitly authorized and run. See `manual-test.md`.
- No push, merge, or PR is authorized.
