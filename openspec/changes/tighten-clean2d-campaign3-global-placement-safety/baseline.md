# Safety-closure baseline

- Branch: `clean2d/campaign-implementation`.
- HEAD: `9facaea Archive completed Clean2D Campaign 3 placement`.
- Pre-existing untracked files: `docs/ui-modernization/PLAN.md`, `docs/ui-modernization/mockup-ui.html`; not in scope.
- `pytest --collect-only -q`: 1604 tests collected.
- Full `pytest -q`: 1580 passed, 20 skipped, 4 failed.
- Allowed pre-existing failures:
  - `tests/test_clean2d_engine_candidates.py::test_generate_candidates_attempts_rdkit_for_cyclic_graphs`
  - `tests/test_smiles_stereo_import.py::test_chiral_smiles_import_creates_wedge_or_hash`
  - `tests/test_smiles_stereo_import.py::test_amino_acid_chiral_smiles_import_creates_wedge_or_hash`
  - `tests/test_smiles_stereo_import.py::test_tetrandrine_import_preserves_visual_stereo`
- `pytest -q tests/architecture`: 276 passed.
- `python -m compileall -q src tests tools packaging`: passed.
- Required Ruff selector: only the documented historical F401 in `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py`.
- `git diff --check`: passed.
- Exact command output is preserved in `/tmp/chemuson_campaign3_safety_baseline.txt`.
