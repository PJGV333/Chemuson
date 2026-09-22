# Baseline: Complete Clean2D Campaign 4 target improvement

- Branch: `clean2d/campaign-implementation`
- Baseline commit: `9a60968adb67ec47c3877ea36d01aebedab73457` (`Archive Clean2D Campaign 4 rigid multiring layout`)
- Previous implementation commit retained: `23b7046`
- Historical Campaign 3 safety reference: `346b229`
- Initial `git status --short --untracked-files=all`: clean

## Baseline commands and results

Environment used for pytest:

```text
PYTHONPATH=src:. uv run --no-project --with pytest --with PyQt6 --with numpy --with Pillow --with rdkit --with certifi --with PyYAML
```

| Command | Result |
|---|---|
| `python -m compileall -q src tests tools packaging` | exit 0 |
| `PYTHONPATH=src:. ... pytest --collect-only -q` | `1614 tests collected in 0.57s` |
| `PYTHONPATH=src:. ... pytest -q tests/architecture` | `276 passed in 19.47s` |
| `PYTHONPATH=src:. ... pytest -q` | `1590 passed, 20 skipped, 4 failed in 344.42s` |
| `uv run --no-project --with ruff ruff check src tests tools packaging --select F401,F811,F821,E722,E741` | one historical `F401`: unused `math` in `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3` |
| `git diff --check` | exit 0 |

The four full-suite failures are pre-existing and remain unchanged:

1. `tests/test_clean2d_engine_candidates.py::test_generate_candidates_attempts_rdkit_for_cyclic_graphs`
2. `tests/test_smiles_stereo_import.py::test_chiral_smiles_import_creates_wedge_or_hash`
3. `tests/test_smiles_stereo_import.py::test_amino_acid_chiral_smiles_import_creates_wedge_or_hash`
4. `tests/test_smiles_stereo_import.py::test_tetrandrine_import_preserves_visual_stereo`

## Corrective RED

Before production changes, the new corrective test file produced the expected RED result: `4 failed in 0.72s`, covering boolean gate-map semantics, true polycyclic classification, spiro promotion selection, and unsafe fused-substitution rejection.
