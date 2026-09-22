# Campaign 4 baseline

- Branch: `clean2d/campaign-implementation`.
- Baseline commit: `346b229 Archive Clean2D Campaign 3 safety closure`.
- Pre-existing untracked files, outside scope and untouched: `docs/ui-modernization/PLAN.md`, `docs/ui-modernization/mockup-ui.html`.
- Exact baseline output: `/tmp/chemuson_campaign4_baseline.txt`.
- `pytest --collect-only -q`: 1607 tests collected.
- `pytest -q tests/architecture`: 276 passed.
- Full `pytest -q`: 1583 passed, 20 skipped, 4 failed.
- The four permitted historical failures are the candidate-generation cyclic-RDKit test and the three SMILES stereo-import tests documented by Campaign 3.
- `python -m compileall -q src tests tools packaging`: passed.
- Required Ruff selector: only historical F401 in `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py`.
- `git diff --check`: passed.
- `openspec validate --all --strict`: 16 passed, 25 historical Purpose-placeholder failures in unrelated specs; active Campaign 3 safety spec passed.
- Campaign 4 baseline evidence SHALL be generated from this commit before production code changes and SHALL exclude `runtime_ms` from deterministic comparison.
