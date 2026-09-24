# Baseline — 2026-09-24-modernize-ui-app-bar-document-tabs

Capturada el 2026-09-24 en la rama `ui/modernization` (HEAD `b649af4`,
cambio limpio salvo `M .gitignore` documentado), **antes** de tocar la
implementación de la Fase 3.

## Entorno de validación

El worktree usa su propio `.venv` (Python 3.14.7, `chemuson` desde `src/`).
Para la suite se instaló en ese venv `pytest==9.1.1`, `ruff==0.16.8` y
`PyYAML` (faltante; lo requieren los tests de arquitectura). Sin dependencias
de producción nuevas.

```
.venv/bin/python -m pytest -q
.venv/bin/python -m ruff check src tests tools packaging --select F401,F811,F821,E722,E741
.venv/bin/python -m compileall src tests tools packaging
```

## git status --short

```
 M .gitignore
```

(`.gitignore`: línea `*.egg-info/` añadida para el metapackage del venv;
sin efecto en código.)

## python -m compileall src tests tools packaging

```
compileall rc=0
```

## pytest --collect-only -q

```
1677 tests collected in 0.81s
```

## pytest -q (árbol pre-cambio)

Referencia estable pre-cambio (mismo HEAD `b649af4`, resultados de la Fase 2):
**1653 passed / 4 failed / 20 skipped**. Los 4 fallos son los preexistentes
(documentados desde la Fase 1) y se mantienen sin tocar:

```
FAILED tests/test_clean2d_engine_candidates.py::test_generate_candidates_attempts_rdkit_for_cyclic_graphs
FAILED tests/test_smiles_stereo_import.py::test_chiral_smiles_import_creates_wedge_or_hash
FAILED tests/test_smiles_stereo_import.py::test_amino_acid_chiral_smiles_import_creates_wedge_or_hash
FAILED tests/test_smiles_stereo_import.py::test_tetrandrine_import_preserves_visual_stereo
```

Nota metodológica: la captura `pytest -q` de este directorio comenzó sobre el
árbol limpio y finalizó (457 s) ya con parte de la Fase 3 aplicada, por lo que
registró `5 failed, 1652 passed, 20 skipped`: los mismos 4 fallos
preexistentes + 1 fallo espurio de inventario
(`tests/test_ui_svg_icons.py::TestSvgInventory::test_all_static_files_parse_with_valid_viewbox`,
que aún esperaba 55 SVG; se actualizó a 63 con los 8 SVG nuevos). No se
detectó ninguna otra regresión; la suite definitiva post-cambio se ejecuta de
nuevo completa y se registra en `tasks.md` (verificación 16).

## ruff check src tests tools packaging --select F401,F811,F821,E722,E741

```
F401 `math` imported but unused
 --> tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3:8
Found 1 error. (preexistente, no se toca)
```

## Resultado post-cambio (tras la suite completa)

Suite completa post-cambio (mismo entorno):

```
4 failed, 1676 passed, 20 skipped in 446.89s
```

Los 4 fallos son **exactamente** los 4 preexistentes listados arriba; nada
nuevo. 1676 = 1653 (Fase 2) + 23 tests nuevos
(`tests/test_ui_app_bar_tabs.py`). Sin regresiones.

Checks estáticos post-cambio: `compileall` rc=0; ruff scoped → solo el F401
preexistente; `git diff --check` limpio; `openspec validate
2026-09-24-modernize-ui-app-bar-document-tabs --strict` → válido.
