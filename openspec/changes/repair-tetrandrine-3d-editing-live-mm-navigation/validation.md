# Validation report

## Baseline
See `baseline.md`. Baseline focused tests passed (164 passed, 2 skipped); compileall and collection passed. Repository-wide Ruff remains red for one pre-existing F401. Full pytest remains unrun because of the known Qt/SIP SIGSEGV.

## Per-commit validation

| Commit track | Chemistry / interaction / animation tests | Architecture | Compileall | Scoped Ruff | OpenSpec / diff check |
|---|---|---|---|---|---|
| OpenSpec gate | baseline captured | pending | baseline pass | baseline pass scoped | `openspec validate ... --strict` valid; `git diff --check` pass |
| Local template | pending | pending | pending | pending | pending |
| Direct 3D drawing | pending | pending | pending | pending | pending |
| Live MM playback | pending | pending | pending | pending | pending |
| Camera navigation | pending | pending | pending | pending | pending |

## Acceptance caveats

- No visual success may be inferred from a mock, offscreen widget, mesh count, `render_requested`, or a successful OpenGL upload alone.
- The user-local tetrandrina template has no explicit stereo. The implementation must not assign stereochemistry; nonconverged minimization must remain reported as such.
- The full pytest suite is not rerun while the historical SIGSEGV remains open.
- Real GUI/OpenGL validation is pending unless explicitly authorized and run. See `manual-test.md`.
- No push, merge, or PR is authorized.
