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

## Baseline de la corrección post-push (tras `ca39a3a`, gate manual Fase 7)

Estado registrado **antes** de tocar código para las 3 correcciones
(QSS `opacity`, `iconSize` de Plantillas, onboarding "No volver a mostrar").
Rama `ui/modernization`, HEAD `ca39a3a`, árbol limpio. Logs en
`/tmp/f7fix/baseline_*.log`.

- `git status --short`: vacío.
- `python -m compileall src tests tools packaging`: OK (sin errores).
- `pytest --collect-only -q`: `1855 tests collected`.
- `pytest -q`: `5 failed, 1830 passed, 20 skipped` (los 5 fallos son la clase
  RDKit/async preexistente, idénticos a los del baseline original; `1830` =
  `1813` + `17` tests nuevos de la Fase 7 ya comprometidos en `ca39a3a`).
- `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`:
  `Found 1 error` — F401 `math` en
  `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3` (preexistente).

Criterio de parada de **esta corrección**: misma PC, mismos 5 fallos
preexistentes (idénticos); `1830 passed + <tests nuevos>` passed, `20
skipped`; Ruff scoped solo el F401 preexistente. Los tests nuevos añadidos
por la corrección: `test_qss_disabled_states_use_no_opacity`,
`test_qss_disabled_states_have_supported_visual_change` (sustituyen a
`test_qss_disabled_states_have_opacity`), `test_templates_tree_sets_thumbnail_icon_size`
y 3 de onboarding (`..._close_with_no_more_persists`,
`..._close_without_no_more_not_persisted`, `..._reappears_when_not_completed`);
neto `+5` tests.

## Baseline de la corrección HiDPI del render de thumbnails (post-push `f0d2a72`)

Estado registrado **antes** de tocar código para el único bug pendiente del
gate manual: en `template_preview_icon()` el DPR se aplicaba dos veces
(backing `logical × dpr` + `setDevicePixelRatio(dpr)` antes de pintar +
`painter.scale(dpr, dpr)` → `dpr²`), lo que sobredimensiona y recorta la
estructura a 200 % (evidencia `templates_after_200.png`).

Rama `ui/modernization`, HEAD `f0d2a72`, árbol limpio. Logs en
`/tmp/f7fix2/baseline_*.log`.

- `git status --short`: vacío.
- `python -m compileall src tests tools packaging`: OK (sin errores).
- `pytest --collect-only -q`: `1860 tests collected in 0.61s`.
- `pytest tests/test_ui_polish.py -q`: `22 passed in 5.18s`.
- `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`:
  `Found 1 error` — F401 `math` en
  `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`
  (preexistente, fuera de alcance).
- Suite completa en `f0d2a72` (referencia del commit): `5 failed` (clase
  RDKit/async preexistente), `1835 passed`, `20 skipped`.

Criterio de parada de esta corrección: mismos 5 fallos preexistentes
(idénticos), `1835 + 1` passed, 20 skipped, Ruff solo el F401 preexistente,
`git diff --check` OK y OpenSpec strict válido. Solo se toca el orden del DPR
en el render del thumbnail; `iconSize` (88×56), grafo, átomos, enlaces,
molblocks y geometría química permanecen intactos.

## Baseline de la intervención post-push (tras `fd0c342`, gate manual KDE/Wayland)

Estado registrado **antes** de tocar código para los tres puntos de esta
intervención: (1) render de la máscara del onboarding, (2) presentación
theme-aware de la tarjeta y (3) activación de Plantillas con un solo clic.

Rama `ui/modernization`, HEAD `fd0c342`, árbol limpio, local == remoto.
Logs en `/tmp/f7fix3/baseline_*.log`.

- `git status --short`: vacío.
- `python -m compileall src tests tools packaging`: OK (sin errores).
- `pytest --collect-only -q`: `1861 tests collected in 0.46s`.
- `pytest tests/test_ui_polish.py tests/test_template_dock.py -q`:
  `26 passed in 5.47s`.
- `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`:
  `Found 1 error` — F401 `math` en
  `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`
  (preexistente, fuera de alcance).

Defecto visual de partida (evidencia `onboarding_step1.png` antes de la
corrección): **1668 píxeles casi negros** (franjas/bordes de la máscara con
`CompositionMode_Clear`); tras la corrección por resta de caminos: **0**.

Criterio de parada de esta intervención: mismos 5 fallos preexistentes
(idénticos), `1861 + 10` tests (`1871 collected`), `26 + 10 = 36` en los
archivos dirigidos, Ruff solo el F401 preexistente, `git diff --check` OK y
OpenSpec strict válido. No se toca química, canvas, escena, hit-testing,
grafo molecular ni molblocks de las plantillas; la limpieza de
química/geometría de plantillas queda como backlog separado (design.md D11).
