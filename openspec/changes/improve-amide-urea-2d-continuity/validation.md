# Validación — continuidad 2D de amidas y ureas

## Verificación focal

Entorno automatizado Qt: `QT_QPA_PLATFORM=offscreen`.

| Verificación | Resultado |
|---|---|
| `python -m compileall src tests tools packaging` | Exit 0. |
| `pytest --collect-only -q` | 2,592 tests recolectados; 2,592 IDs del manifiesto esperados, iguales sin faltantes ni sobrantes. |
| `pytest -q` en la matriz focal (geometría inicial, hidrógenos/labels, cadenas, canvas/bonds/junctions, cuñas/dobles, persistencia, Visual Lab, presets, arquitectura y contratos CI) | **216 passed en 8.99 s**. Incluye **14** regresiones nuevas de urea/amida/amina/labels/aromáticos. |
| Ruff scoped: `ruff check src/chemuson/gui/canvas/canvas_bond_geometry.py src/chemuson/gui/canvas/canvas_text.py tests/test_amide_urea_2d_continuity.py --select F401,F811,F821,E722,E741` | All checks passed. |
| `openspec validate improve-amide-urea-2d-continuity --strict --no-interactive` | Change valid. |
| `git diff --check` | Exit 0 (comprobado antes del commit). |

Las nuevas pruebas cubren `Ar–NH–C(=O)–NH–CH3` por eventos Qt, acetanilida, contextos primario/secundario/terciario, amida/urea por conectividad y orden, N aromático, coordinación, carga y valencia alta; preview/commit, drag dirigido, `NH2`/`NH`, no movimiento de centros, dimensiones de enlace, undo/redo, CMSN y zoom/DPR/temas. Las regresiones existentes comprueban el despeje de labels en estilos plain/bold/stereo y uniones.

## Comparación antes/después

La reproducción anterior se reejecutó en un worktree detached temporal del snapshot previo (`db996cef`); el JSON generado coincidió byte por byte con `evidence/before/urea-baseline-metrics.json`. El worktree temporal se eliminó. La reproducción posterior `evidence/after/reproduce_urea_after.py` genera un JSON idéntico al artefacto posterior committed; sus imágenes también se reproducen exactamente.

`Ar–NH–C(=O)–NH–CH3`: 11 átomos, 11 enlaces, referencia 40.0; las cinco distancias del fragmento siguen en 40.0. Los ángulos con los drags capturados siguen entre 119.485807° y 120.514193°. Los endpoints de preview y commit coinciden; snapshot del grafo idéntico tras undo/redo y save/load.

| Enlace | Trim antes → después | Path visible antes → después | Cambio visible |
|---|---:|---:|---:|
| Ar–N1H | 5.631875 → 5.631875 | 34.368125 → 34.368125 | 0.000000 |
| N1–C(=O) | 15.620745 → 8.200844 | 24.379255 → 31.799156 | +7.419902 |
| C=O | 8.204943 → 8.204943 | 31.795057 → 31.795057 | 0.000000 |
| C(=O)–N2H | 15.743313 → 8.372784 | 24.256687 → 31.627216 | +7.370529 |
| N2H–CH3 | 5.679473 → 5.644793 | 34.320527 → 34.355207 | +0.034679 |

La corrección recupera ~7.4 unidades de trazo en cada NH–carbonilo; no altera los centros, la carbonila, la convención de longitud ni la silueta aromática. Capturas/light/dark/1×/1.25× y métricas quedan en `evidence/before/` y `evidence/after/`; el test paramétrico también renderiza 1×/1.25×/1.5× con DPR 1/1.25/2. La hoja del lienzo sigue blanca y las capturas de papel son idénticas por hash entre temas claro/oscuro; los tests confirman estabilidad del clearance, no un cambio cosmético de papel por tema.

Las imágenes son capturas Qt/QImage offscreen y se revisaron visualmente en esta sesión. No equivalen a aceptación manual nativa en CachyOS; no se afirma esa revisión.

## Fallos globales basales, sin ocultar

- `timeout --signal=TERM 600s pytest -q` terminó con exit 124 a los 600 s, progreso ~30%, mostrando un único fallo y skips antes del timeout.
- `pytest -x -q` da **380 passed, 1 failed en 36.27 s**. El primer fallo es el mismo basal: `tests/test_beta_compact_shell_and_tabs.py::test_view_menu_submenus_reuse_actions_and_remove_obsolete_auxiliary_option`, que espera omitir `Visor molecular 3D` aunque el menú la ofrece. stderr: `QOpenGLWidget is not supported on this platform`. El baseline pre-cambio había llegado a 366 pasadas antes del mismo test; las 14 nuevas pasan antes de él.
- Ruff global requerido sigue fallando sólo por el F401 basal `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`. No se corrige una deuda ajena para enmascarar el resultado.

Por tanto, **la matriz focal está verde; la suite monolítica no se declara verde**. No se cambió Clean2D, química del grafo, serialización, renderer editorial general ni interfaz 3D.
