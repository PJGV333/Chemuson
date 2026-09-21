# Campaign 2 baseline

## Repository identity

- Branch: `clean2d/campaign-implementation`
- Base commit: `68835fb Archive Clean2D Campaign 1 specification`
- Baseline worktree: clean before Campaign 2 changes

## Commands and exact results

| Command | Exit | Result |
|---|---:|---|
| `git status --short --untracked-files=all` | 0 | clean |
| `pytest -q tests/test_clean2d_multilayer_constraints.py tests/test_block_unwrap_depiction.py tests/test_clean2d_complex_policy.py tests/test_clean2d_complex_preserve.py` | 0 | 21 passed, 1 skipped in 2.33s |
| `pytest -q tests/architecture` | 0 | 276 passed in 18.47s |
| `python -m compileall -q src tests tools packaging` | 0 | passed |
| `git diff --check` | 0 | passed |

## Baseline contract

The Campaign 2 implementation must preserve existing decomposition, multilayer
constraint, block unwrap, complex-policy, architecture, and geometry behavior.
Campaign 1 evidence remains the before/after corpus input:

`openspec/changes/archive/2026-09-21-clean2d-campaign1-benchmark-observability/evidence/baseline.json`
