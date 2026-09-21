# Strengthened Campaign 2 baseline

- Branch: `clean2d/campaign-implementation`
- Base commit: `a348baf Archive Clean2D Campaign 2 specification`
- Worktree: clean before this change
- `pytest -q tests/test_clean2d_topology_decomposition.py`: 2 passed in 0.06s
- Existing topology/complexity focus: 21 passed, 1 skipped in 2.28s
- `pytest -q tests/architecture`: 276 passed in 18.19s
- `python -m compileall -q src tests tools packaging`: exit 0
- `git diff --check`: exit 0
- Last full-suite baseline: 1,541 passed, 20 skipped, 4 failed; the four
  failures are the documented candidate/stereo-import failures.
