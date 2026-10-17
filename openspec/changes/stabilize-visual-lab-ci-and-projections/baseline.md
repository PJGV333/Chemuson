# Baseline — CI y Molecular Depiction Visual Lab

## Estado previo

- Rama: `feat/molecular-depiction-visual-lab`.
- HEAD: `87e298c46535e1dffb0477f36412b3dedc87bd98`; `origin/feat/molecular-depiction-visual-lab` apuntaba al mismo SHA.
- `git status --short`: vacío antes de iniciar esta propuesta.
- Transcripción obligatoria: `baseline-command-output.txt` (SHA-256 `f02be35e8ed22f32d92bbb1c6a021d2e3fd63c8b4d0941bde7b4e461ce6dd078`). Los comandos se ejecutaron antes de editar documentos/código.

## Baseline local

| Comando | Resultado |
|---|---|
| `git status --short` | PASS, limpio. |
| `python -m compileall src tests tools packaging` | PASS. |
| `pytest --collect-only -q` | PASS, 2,341 node IDs en 0.62 s. |
| `pytest -q` | 2,320 passed, 20 skipped, 1 fallo en 1,310.93 s (21:50). Falló `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend`: se observó `worker.finished` pero faltaron `thread.finished` y `job_finished`; queda como fallo basal de lifecycle Qt, distinto de los dos timeouts de CI. |
| `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` | Exit 1 por el F401 histórico `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; fuera de alcance. |

## GitHub Actions observado

Run: [38996089270](https://github.com/PJGV333/Chemuson/actions/runs/38996089270), workflow `test`, head `87e298c46535e1dffb0477f36412b3dedc87bd98`, Ubuntu 24.04, Python 3.11, `QT_QPA_PLATFORM=offscreen`. El plan SHA-256 `b274734b643c29e54c12b379375d2cf0570cba3cdc0e28c86e800b8eeea28c29` asignó 2,341 IDs como `[293,293,293,293,293,292,292,292]`; hash de colección `a7d09fbaff3bb32248792d66971323636799381e19dab736e52cef3d9e34744f`.

| Shard | IDs | Duración del runner | Resultado |
|---:|---:|---:|---|
| 0 | 293 | 34.466 s | PASS |
| 1 | 293 | 47.327 s | PASS |
| 2 | 293 | 35.954 s | PASS |
| 3 | 293 | 302.024 s | TIMEOUT 300 s; exit -15 / SIGTERM; último `tests/test_command_palette.py::test_palette_closes_after_execution` (~70%). |
| 4 | 293 | 108.171 s | PASS |
| 5 | 292 | 302.041 s | TIMEOUT 300 s; exit -15 / SIGTERM; último `tests/test_ui_polish.py::test_dark_text3_meets_aa` (~83%). |
| 6 | 292 | 73.698 s | PASS |
| 7 | 292 | 20.137 s | PASS |

Windows, Flatpak y plan de colección terminaron correctamente. `pytest-summary` falló correctamente: detectó que los dos shards no terminaron. Los JSON, logs, JUnit que sí se alcanzó a crear y el plan exacto están versionados en `baseline/ci-artifacts/`; `github-run-full.log.gz` conserva una compresión determinista del log completo descargado (el hash del original y del comprimido están en su README).

## Diagnóstico de los timeouts

La distribución original asignaba casi el mismo **número** de pruebas por shard y mantenía indivisibles los archivos, pero `_partition_by_file()` no usaba coste real. Los tiempos de ejecución (instrumento de este baseline: atributos `time` de JUnit en seis shards completos y diferencias entre las marcas UTC de finalización de pruebas en los dos logs incompletos) muestran concentración de trabajo:

- Shard 3: `test_clean2d_baseline_diff_review_policy.py` 11 pruebas/~72 s; `test_clean2d_baseline_report_script.py` 8/~40 s; `test_command_palette.py` completó 20 de 27/~67 s observados (proyección del archivo completo ~91 s); `test_chemname_acceptance.py` 76/~46 s; `test_clean2d_wedge_stereo.py` 35/~38 s; `test_beta_save_shortcut_menu_access.py` 20/~28 s. Los seis grupos ya suman alrededor de 291 s antes de todos los archivos ligeros restantes.
- Shard 5: `test_clean2d_baseline_reports.py` 15/~87 s; `test_ui_app_bar_tabs.py` 24/~55 s; `test_beta_compact_shell_and_tabs.py` 28/~55 s; `test_ui_polish.py` completó 4 de 38/~5 s hasta el timeout (proyección del archivo completo ~55 s); otros grupos de UI/assistant añaden más coste. Tres archivos principales ya suman ~198 s.
- Los seis shards completados variaron de 20.1 a 108.2 s con cuentas casi iguales (292–293), confirmando fuerte heterogeneidad del coste por archivo.

### Conclusión

La causa primaria confirmada es **balanceo por cantidad de node IDs, no por duración**, con carga acumulada de archivos Qt/Clean2D particularmente lentos. La ejecución imprimió pruebas progresivamente; el mayor intervalo observado entre resultados consecutivos dentro de los dos jobs incompletos fue ~14.2 s, no hubo pausa prolongada comparable a 120 s ni volcado `faulthandler_timeout=120`. No hay evidencia en estos logs de una única prueba colgada, lock o bloqueo de recursos no liberados. La lentitud se acumula en muchos casos de UI/Qt y pruebas de baseline que completan; no se identificó una regresión dependiente del orden para los shards que expiraron. Un timeout externo mató el proceso al llegar a 300 s, dejando el reporte/JUnit incompleto, no un fallo de aserción. El timeout no se sube en esta propuesta: primero se separan los archivos por pesos medidos y se vuelven a observar los ocho jobs.

El test Qt de CompChem que falla en el `pytest -q` local completo es una condición basal independiente; se reporta por separado y no explica los logs de timeout de Actions. Durante la validación de esta propuesta también se confirmó una falla sensible al orden en `tests/test_molecular_assistant_reference_resolution.py::test_ai_success_with_equivalent_nonidentical_smiles_is_verified`: el nodo pasa aislado y con la familia 3D/Molecular Assistant, pero falla después de ejecutar conjuntamente las suites UI `test_beta_compact_shell_and_tabs.py` y `test_beta_save_shortcut_menu_access.py`. El mismo prefijo ordenado en el worktree limpio del SHA basal `87e298c` reprodujo `69 passed, 1 failed` (79.38 s); no se modificaron el controller ni esas pruebas. Se documenta como falla basal de lifecycle Qt no causada por esta propuesta, y no se oculta ni se descuenta del resultado final.

## UI y proyecciones antes de cambio

- `ProductionCanvasAdapter.create()` crea `ChemusonCanvas`, fija `sceneRect(graph_scene_rect)` y centra en esa geometría; en el primer `showEvent`, `QTimer.singleShot(0, self._center_canvas_initial)` ejecuta el método heredado, que vuelve a llamar `center_on_paper()` sobre el tamaño/punto de centrado de producción. Esto confirma una ruta plausible para la escena vacía del laboratorio; se probará sobre viewport real después de procesar el timer.
- `ProjectionGeometry.to_molgraph()` produce únicamente `MolGraph` con coordenadas/enlaces y no asigna stereo; `Atom`/`Bond` y `.cmsn` no contienen descriptores α/β/D/L, convención Fischer o relación de silla. No se cambiará el modelo químico ni el formato de persistencia para esta compatibilidad.
- Haworth baseline y evidencia visual anteriores están registrados en el OpenSpec histórico. Este cambio generará una comparación nueva, de igual escala y renderer.
