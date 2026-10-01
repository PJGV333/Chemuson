# Baseline — 2026-10-01-modernize-ui-command-palette

Estado registrado antes de modificar el código de esta fase y, de forma
independiente, antes de la **corrección de UX de atajos** (Ctrl+K→Clean2D,
Ctrl+P→paleta). Entorno: `.venv` (Python 3.14), `QT_QPA_PLATFORM=offscreen`,
pytest 9.1.1, ruff 0.16.8. RAM/DOMINIOS de química no tocados.

## Referencia: HEAD `3c8d992` (fase 6 inicial, antes de la corrección)

Suite completa (`pytest -q`):

```
5 failed, 1813 passed, 20 skipped in 1324.81s (0:22:04)
EXIT=1
```

Los 5 fallos son de la clase RDKit química/async preexistente (documentada en
la baseline pre-Fase 6 y en `AGENT_REPORT.md`); ninguno importa nada de la
paleta, el AppBar ni la migración de atajos:

```
FAILED tests/test_clean2d_engine_candidates.py::test_generate_candidates_attempts_rdkit_for_cyclic_graphs
FAILED tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend
FAILED tests/test_smiles_stereo_import.py::test_chiral_smiles_import_creates_wedge_or_hash
FAILED tests/test_smiles_stereo_import.py::test_amino_acid_chiral_smiles_import_creates_wedge_or_hash
FAILED tests/test_smiles_stereo_import.py::test_tetrandrine_import_preserves_visual_stereo
```

`test_compchem3d_dock` (async) **pasa al re-ejecutarse en aislamiento** → flaky.

## Referencia: baseline pre-Fase 6 (HEAD `02c8ea8`, antes de crear la paleta)

Suite completa:

```
4 failed, 1787 passed, 20 skipped in 1222.29s (0:20:22)
EXIT=1
```

Mismos 4 fallos RDKit de la clase preexistente (stereo import ×3 + clean2d
candidates). La diferencia `1813 − 1787 = 26` es la mayoría de los tests nuevos
de `tests/test_command_palette.py`; el 5º fallo de `3c8d992` es el flaky
`compchem3d`.

## Comandos de baseline (protocolo AGENTS.md §1.2)

- `git status --short`
- `python -m compileall src tests tools packaging`
- `pytest --collect-only -q`
- `pytest -q`
- `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`

Ruff (scoped) sobre los archivos de la fase: `All checks passed!` (sin F401/F811/
F821/E722/E741 nuevos). El único F401 del repo completo es preexistente
(`tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py`, fuera de
alcance).
