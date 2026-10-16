# Baseline — Molecular Depiction Visual Lab

## Git

- Requested base branch: `feat/molecular-depiction-excellence`.
- `git fetch origin feat/molecular-depiction-excellence`: completed successfully.
- Remote branch HEAD fetched: `ee5294800c5779f7451e625d8089025e8e6632ab` (`docs(openspec): record depiction delivery status`, 2026-10-16).
- Before branch creation, local branch was `release/v0.3.0-beta.1-prep`, clean; no files were carried over.
- Created `feat/molecular-depiction-visual-lab` directly from that fetched SHA; no merge/rebase/cherry-pick. No protected branch was changed.
- Baseline `git status --short`: empty.

## Exact command transcript

`baseline-command-output.txt` contains the exact output and exit status of the required commands run with repository `.venv/bin` first on `PATH` (Python 3.14.8, pytest 9.1.1, Ruff 0.16.10). The initial unactivated system-Python attempt lacked pytest/Ruff; it is not used as the environment baseline and is preserved in the session transcript, not this file.

| Required command | Baseline |
|---|---|
| `git status --short` | Exit 0, empty. |
| `python -m compileall src tests tools packaging` | Exit 0. |
| `pytest --collect-only -q` | Exit 0, 2,285 tests collected in 1.74 s. |
| `pytest -q` | Externally bounded at 600 s; terminated with timeout (exit 124) at approximately 34% progress. No assertion traceback was captured before termination. This follows known historical Qt suite duration; no blind retry. |
| `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` | Exit 1 for one pre-existing F401: unused `math` in `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; unrelated and not modified. |

## Relevant pre-existing visual/CI context

- Existing focused gallery tools: `tools/visual/render_bond_junction_gallery.py` and `tools/visual/render_wedge_double_label_gallery.py`, both use actual Qt production items/canvas and generate evidence.
- No separate reusable orbital calibration application exists. `tools/visual/verify_orbitals_rendering.py` is an orbital-render verification utility only; it is out of scope.
- The parent branch's `openspec/changes/refine-wedge-double-bond-junctions/validation.md` documents shard 4's order-dependent failure in `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend` and says shards 5–8 were not run in that campaign. It is historical, not attributed to this change.
- Baseline and verification will use bounded focused groups and the repository's 8-shard CI utility rather than repeatedly running the monolithic suite.
