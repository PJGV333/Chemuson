# Validación — Clean2D con cuñas

## Identidad y alcance

- Rama: `clean2d/wedge-bond-campaign`, creada desde `origin/release/v0.3.0-beta.1-prep` en `9bc8192` (`docs(release): record beta preview acceptance gates`).
- Cambio OpenSpec: `fix-clean2d-wedge-stereo-layout`.
- Código de producción tocado: **un único archivo**, `src/chemuson/clean2d/local_graph_cleaner.py` (+99 / −21). No se tocaron ChemName, ChemIO, `core`, GUI, menús, persistencia `.cmsn`, CompChem, packaging, workflows ni dependencias.
- No se modificó `architecture/modules.yml`: no hay módulos nuevos ni cambios de import entre paquetes; `local_graph_cleaner` sigue siendo el módulo ya catalogado y la suite de arquitectura sigue en verde.
- La Fase A (`fix/beta-save-shortcut-menu-access`, remoto `47c588d`) no es ancestro de esta rama y no se transportó nada de ella.

### Archivos del cambio

- `src/chemuson/clean2d/local_graph_cleaner.py`
- `tests/test_clean2d_wedge_stereo.py` (nuevo, 20 casos)
- `tests/ci/expected_pytest_nodeids.txt` (2154 → **2174** node IDs, regenerado con `tools/ci_pytest_shards.py update-expected --confirm`)
- `openspec/changes/fix-clean2d-wedge-stereo-layout/{.openspec.yaml,README.md,proposal.md,design.md,tasks.md,baseline.md,validation.md,specs/clean-2d/spec.md}`

## Decisiones implementadas

- `stereo_layout_signature()` pasa de signos de producto vectorial por par de vecinos al **orden cíclico canónico** (mínimo por rotación) de los vecinos estructurales reales de cada centro estereogénico, con el estilo de enlace de cada vecino.
- Los vecinos y enlaces **de frontera** entran en la firma usando las coordenadas del grafo cuando la propuesta no los incluye (`_structural_bonds_touching_selection`, `_layout_point`, `_structural_neighbour_styles`, `_cyclic_neighbour_order`).
- Se conserva el bloque de enlaces con estereoquímica explícita (ID, extremos, estilo, `stereo`), ahora también para enlaces de frontera.
- No se cambian el ranking de candidatos, la política de complejidad (`multiple_stereo_centers` → `preserve_only`), el vocabulario estable de motivos de rechazo ni los umbrales geométricos.

## Red / green

- **Sin el parche** (`git stash` del único archivo de producción), `tests/test_clean2d_wedge_stereo.py`: **9 failed, 11 passed**. Fallan los dos casos de cuña simple en `quick` y `publication`, los dos de determinismo de cuña simple, la invariancia de rotación/re-encuadre, la cobertura de vecinos de frontera y el camino GUI.
- **Con el parche**: **20 passed in 4.24s**.

## Matriz observada en el motor (sonda desechable, misma definición de casos)

| caso | antes | después |
|---|---|---|
| 1 control acíclico sin estereo | `applied` / `rdkit_isolated` | `applied` / `rdkit_isolated` |
| 2 una cuña sólida | `failed-controlled` / ninguna | **`applied` / `rdkit_isolated`** |
| 3 una cuña rayada | `failed-controlled` / ninguna | **`applied` / `rdkit_isolated`** |
| 4 centro tetraédrico cuña + hash | `applied` / `rdkit_isolated` | `applied` / `rdkit_isolated` |
| 5 varios centros quirales | `applied` / `rdkit_isolated` | `applied` / `rdkit_isolated` |
| 6 selección con centro quiral | `failed-controlled` / ninguna | `failed-controlled` / ninguna |
| 7 selección con cuña en frontera | `failed-controlled` / ninguna | `failed-controlled` / ninguna |

Lectura: los casos 2 y 3 (el defecto reportado) pasan de no hacer nada a limpiarse. Los casos 6 y 7 siguen sin candidato aplicable porque la selección deja estereoquímica en la frontera; el contrato que se fija es que **nunca** se aplica una inversión, no que siempre exista re-encuadre. La sonda confirma que en los siete casos la configuración dibujada (orden cíclico + estilos) queda intacta.

## Matriz observada en la GUI (ventana real, `QT_QPA_PLATFORM=offscreen`)

Molécula de cuatro átomos con una cuña sólida y, por separado, una cuña rayada:

| acceso | antes | después |
|---|---|---|
| `_on_clean_2d_full` | 0 átomos movidos, `Limpieza 2D omitida: geometria_actual_necesita_reconstruccion` | 4 átomos movidos, `Estructura 2D optimizada`, deshacible |
| `_on_clean_2d` | 0 átomos movidos, mismo mensaje | 4 átomos movidos, `Estructura 2D optimizada`, deshacible |
| `_on_clean_2d_publication` | 0 átomos movidos, mismo mensaje | 4 átomos movidos, `Estructura 2D optimizada`, deshacible |
| `_on_clean_2d_propose` | `La estructura debe optimizarse antes de proponer conformeros 2D` | tras una limpieza previa: `Estructura 2D alternativa propuesta`, deshacible |

En todos los casos el orden cíclico alrededor del centro permanece igual y los estilos `wedge`/`hashed` no cambian. Undo restaura las coordenadas originales (verificado en la sonda y fijado en `test_gui_clean_2d_full_moves_wedge_molecule_and_undo_restores_it`).

## Evidencia de validación (comandos acotados externamente)

- `timeout 180 python -m compileall -q src tests tools packaging`: exit 0.
- `timeout 300 python -m pytest --collect-only -q`: exit 0; **2174 tests collected**.
- `timeout 120 ruff check src/chemuson/clean2d/local_graph_cleaner.py tests/test_clean2d_wedge_stereo.py --select F401,F811,F821,E722,E741`: **All checks passed!**
- `timeout 180 ruff check src tests tools packaging --select F401,F811,F821,E722,E741`: exit 1 por el único hallazgo preexistente fuera de alcance (`tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`, F401 `math`). No se modificó.
- `timeout 580 python -m pytest -q tests/architecture`: **280 passed in 11.90s**.
- Grupo Clean2D completo con el parche (`tests/test_clean2d_*.py`, `tests/test_clean_2d_*.py`, `tests/test_gui_clean2d_geometry.py`, `tests/test_wedge_geometry.py`): **397 passed in 191.27s** (377 en la base + 20 nuevos).
- Suite completa por shards, cada shard con timeout externo de 300 s:
  - `timeout 300 python tools/ci_pytest_shards.py plan --shards 8 --expected tests/ci/expected_pytest_nodeids.txt --output /tmp/shardB/plan.json` → distribución `[272, 272, 272, 272, 272, 272, 271, 271]`, `sha256=fc4037e59e1e65ecb05596dceb027201f563e39000a81f837d7a59f836de9aa0`.
  - Ocho ejecuciones `run --shard-index 0..7 --timeout-seconds 300`: exit 0 en las ocho, `timed_out=false`, `junit_matches_assignment=true`.
  - `timeout 300 python tools/ci_pytest_shards.py verify --artifacts /tmp/shardB --expected tests/ci/expected_pytest_nodeids.txt --shard-count 8 --plan-status success --shards-status success`: **`pytest shard coverage verified: 2174 unique node IDs, 8 complete successful shards`**.
  - Recuento JUnit agregado: **2174 tests, 2154 passed, 20 skipped, 0 failures, 0 errors**.
- `timeout 180 openspec validate fix-clean2d-wedge-stereo-layout --strict --no-interactive`: **valid**.
- `git diff --check`: sin salida.

## Fallos intermitentes previos

El fallo intermitente de CompChem observado en campañas anteriores no apareció en las ocho shards de esta validación. No se tocó CompChem ni se relajó ningún contrato.

## Límites y trabajo no realizado

- No se persigue que las selecciones parciales con estereoquímica en la frontera (casos 6 y 7) obtengan un candidato: se documentan como resultado seguro admisible. Mejorar el **mensaje** de ese caso (hoy `geometria_actual_necesita_reconstruccion`, que no menciona la causa estereoquímica) queda como seguimiento explícito, porque implicaría ampliar el vocabulario estable de motivos de rechazo.
- No se toca la política que deja `preserve_only` ante dos o más centros quirales.
- Las sondas desechables (`/tmp/wedge_probe.py`, `/tmp/wedge_analysis.py`, `/tmp/wedge_gui_probe.py`) no se versionan: su contenido útil quedó convertido en `tests/test_clean2d_wedge_stereo.py`.

## Estado Git

- Commits locales en la rama: `946de07` (`fix(clean2d): preserve wedge configuration by cyclic neighbour order`) y `4fd9dd6` (`docs(openspec): record clean2d wedge stereo campaign and evidence`).
- Publicación: sólo `origin/clean2d/wedge-bond-campaign`, SHA remoto `4fd9dd64bed21bc51ce6ce0341c8e807b1d496b7`. Sin merge a `release/v0.3.0-beta.1-prep` ni `main`, sin tags, sin releases.

## CI en la rama publicada

- Run de `test.yml` para `4fd9dd6`: id **38819097932**, conclusión **success**.
- Jobs en verde: `pytest-plan`, los ocho `pytest shard 0..7`, `pytest-summary`, `windows-smoke` y `flatpak-smoke`.
- `pytest-summary` registró: `pytest shard coverage verified: 2174 unique node IDs, 8 complete successful shards`.
