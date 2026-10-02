# Baseline — 2026-10-01-modernize-ui-polish

Estado registrado **antes de modificar el código de esta fase**.
Rama `ui/modernization`, HEAD local/remoto `50a8241`
("Restore Clean2D shortcut and finalize command palette"), árbol limpio.
Entorno: `.venv` (Python 3.12.12), `QT_QPA_PLATFORM=offscreen` (vía conftest),
pytest 9.x, ruff 0.x. Logs completos en `/tmp/f7_compile.log`,
`/tmp/f7_collect.log`, `/tmp/f7_pytest.log`, `/tmp/f7_ruff.log`.

> Comparación solo baseline/final de **esta PC** (no contra otras máquinas).

## Referencia histórica (otra PC, orientativa)

```
1813 passed / 20 skipped / 5 failed
```

## Esta PC — antes de Fase 7

### `git status --short`
```
(vacío)
```

### `python -m compileall src tests tools packaging`
```
OK (sin errores; 88 ficheros escritos)
```

### `pytest --collect-only -q`
```
1838 tests collected in 0.50s
```

### `pytest -q`
```
5 failed, 1813 passed, 20 skipped in 883.95s (0:14:43)
EXIT=1
```

Los 5 fallos son la clase RDKit/async preexistente ya documentada
(`AGENT_REPORT.md`); ninguno importa la UI, la paleta ni los tokens:

```
FAILED tests/test_clean2d_engine_candidates.py::test_generate_candidates_attempts_rdkit_for_cyclic_graphs
FAILED tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend
FAILED tests/test_smiles_stereo_import.py::test_chiral_smiles_import_creates_wedge_or_hash
FAILED tests/test_smiles_stereo_import.py::test_amino_acid_chiral_smiles_import_creates_wedge_or_hash
FAILED tests/test_smiles_stereo_import.py::test_tetrandrine_import_preserves_visual_stereo
```

### `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`
```
Found 1 error.
F401 `math` imported but unused
 --> tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3:8
```

Este F401 es **preexistente** (no introducido por Fase 7) y vive en un test de
Clean2D fuera del alcance de esta fase; se deja intacto y se documenta aquí.
No se modifica (AGENTS.md §2.1: sin refactor oportunista).

## Criterio de parada de esta fase

Al finalizar, la suite SHALL coincidir con este baseline **en la misma PC**:
- mismos 5 fallos preexistentes (idénticos);
- `1813 passed + <tests nuevos>` passed, `20 skipped`;
- Ruff scoped: el único error preexistente (F401 en test de Clean2D) sin
  errores nuevos.
