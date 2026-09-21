# Campaign 1 baseline

## Repository identity

- Branch: `clean2d/campaign-implementation`
- Base commit: `31c565c8dfe15f2b4c2287eb4ece56e2b67d9732` (`Refine Clean2D campaign policy`)
- Baseline worktree: clean before Campaign 1 changes

## Commands and exact baseline results

| Command | Exit | Result |
|---|---:|---|
| `git status --short` | 0 | clean |
| `python -m compileall src tests tools packaging` | 0 | passed; verbose `Listing ...` output intentionally omitted here |
| `pytest --collect-only -q` | 0 | 1,559 tests collected |
| `pytest -q` | 1 | 1,535 passed, 20 skipped, 4 failed in 289.59s |
| `pytest -q tests/architecture` | 0 | 276 passed in 18.23s |
| `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` | 1 | one pre-existing `F401` in `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3` (`math`) |
| `openspec validate --all --strict` | 127 | `openspec: command not found` in the original baseline environment |
| `git diff --check` | 0 | passed |

## Exact baseline failures

1. `tests/test_clean2d_engine_candidates.py::test_generate_candidates_attempts_rdkit_for_cyclic_graphs`
2. `tests/test_smiles_stereo_import.py::test_chiral_smiles_import_creates_wedge_or_hash`
3. `tests/test_smiles_stereo_import.py::test_amino_acid_chiral_smiles_import_creates_wedge_or_hash`
4. `tests/test_smiles_stereo_import.py::test_tetrandrine_import_preserves_visual_stereo`

These failures are pre-existing and remain outside Campaign 1 scope.

## Evidence

Detailed machine-readable evidence is preserved at:

`openspec/changes/clean2d-campaign1-benchmark-observability/evidence/baseline.json`

The JSON report is the reproducible corpus evidence; this Markdown file records the human baseline summary and command outcomes only.
