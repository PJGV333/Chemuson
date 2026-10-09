# Propuesta — pytest y CI reproducibles

## Por qué

La ejecución GitHub Actions #343 sobre `f5f82a63c84c56c9dd3ee07e7c7f56e7eb700bef` falla en seis nodos pytest, aunque cinco pasan al ejecutarse solos. La reproducción local identifica un contrato Clean2D basado erróneamente en la presencia de un candidato tras deduplicación; dos verificaciones del smoke RDKit dependen del orden de importación global; y el aislamiento de preferencias Moleculares sólo cambia `XDG_CONFIG_HOME` después de que Qt puede haber cacheado la ruta de `QSettings`. CompChem y el panel lateral pasan localmente de forma aislada, por lo que sus diagnósticos deben distinguir una carrera de señales de Qt y una geometría transitoria de un fallo estable, no silenciar las aserciones.

El workflow ejecuta toda la suite en un único proceso y no ofrece un límite externo dedicado, resultados por grupo ni prueba mecánica de cobertura de colección.

## Alcance

- Corregir los contratos de los tests Clean2D, preferencias del Asistente Molecular y smoke worker empaquetado sin debilitar la validación.
- Instrumentar las señales `worker.finished`, `QThread.finished` y `job_finished`; medir las pestañas del panel sólo tras estabilizarse el layout, manteniendo intactos sus límites de legibilidad y recorte.
- Añadir plan determinista de shards, manifiesto explícito de IDs pytest, límites externos de cinco minutos, artefactos por shard y resumen que falle ante IDs faltantes, duplicados, señal de proceso, timeout o reporte incompleto.
- Mantener los jobs `windows-smoke` y `flatpak-smoke`.
- Crear una especificación OpenSpec independiente para esta campaña.

## Fuera de alcance

No se cambia código de producción, la lógica química, el diseño del panel, las preferencias de la aplicación, el protocolo del worker, las excepciones/skip de tests, el workflow de releases ni los smoke jobs de Windows/Flatpak. No se ejecuta el pytest monolítico. No se hace merge, PR, tag, release ni publicación de artefactos.

## Módulos probables

`tests/test_clean2d_engine_candidates.py`, `tests/test_molecular_assistant_ui.py`, `tests/test_rdkit_packaged_worker.py`, `tests/test_compchem3d_dock.py`, `tests/test_side_panel.py`, `tests/test_test_workflow_contract.py`, `tests/ci/expected_pytest_nodeids.txt`, `tools/ci_pytest_shards.py` y `.github/workflows/test.yml`.
