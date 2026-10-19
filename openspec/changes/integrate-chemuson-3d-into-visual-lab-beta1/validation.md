# Validación — integración 3D en visual lab Beta 1

## Estado

**Integración local validada de forma focal** en la rama temporal `integration/chemuson-3d-beta1`, con base Visual Lab `43a8a27c45cdabfd157bf3f2a52c9216b24a1552` y fuente 3D `ec486626333ccd2e0ae86671861708b30c508c5b`. No equivale a aceptación completa ni a readiness de release. El baseline anterior al merge está en `baseline.md`.

## Resolución y preservación

- Se resolvieron los seis conflictos manualmente. Se conservaron el `AGENT_REPORT.md` y las superficies de serialización/ChemIO de Beta 1; el catálogo conserva `projections`, `bond_junctions` y añade `viewport3d`, dependencias y pruebas 3D.
- `_rdkit_worker.py` y `rdkit_safe.py` combinan el worker RDKit existente, el modo congelado de una operación y el streaming requerido por 3D. `rdkit_io.py` conserva el importador/exportador de Beta 1. Se integró el reconocimiento MOL del código de enlace estéreo 3 como marca ambigua WAVY/EITHER (igual que el código 4), sin asignar E/Z. La regresión `test_molfile_import_keeps_ambiguous_double_bond_mark_explicit` falló antes con `PLAIN` y pasó después con `WAVY`.
- El builder mantiene la exportación MOL y los menús de Beta 1, y agrega el toggle del viewport. La auditoría de canvas, shell, docks, ventana principal y pruebas CompChem no encontró pérdida de los flujos 2D probados; los cambios 3D guardan geometría en el documento y mantienen el undo/redo probado.
- No se añadieron dependencias runtime. Se actualizaron los node IDs esperados del CI de 2.358 a 2.567 tras revisar los 209 tests nuevos; el manifiesto actualizado registra todas las pruebas, no oculta fallos.

## Verificaciones posteriores al merge

- `python -m compileall -q src tests tools packaging`: PASS.
- Colección completa (`pytest --collect-only -q`): **2.567 tests**.
- ChemIO/worker, packaged-worker, importación estéreo, geometría 3D, viewport y CompChem: **240 passed, 6 skipped**.
- Depiction, proyecciones y regresiones Clean2D acotadas: **159 passed, 3 skipped**.
- Persistencia CMSN, coordinación, roundtrip RDKit, tabs y exportación MOL: **52 passed, 4 skipped**.
- Arquitectura completa: **283 passed**. Contrato de shards: **6 passed**.
- Ruff scoped sobre módulos y pruebas afectadas (`F401,F811,F821,E722,E741`): PASS.
- `openspec validate integrate-chemuson-3d-into-visual-lab-beta1 --strict`: PASS.
- `git diff --check HEAD`: PASS. Escaneo del árbol sin marcadores de conflicto: PASS.
- Planificador CI: inicialmente detectó el manifiesto obsoleto (**2.358 esperados / 2.567 recolectados**); tras actualizar explícitamente `tests/ci/expected_pytest_nodeids.txt`, el plan de ocho shards pasó. Distribución de IDs: `[297,239,323,359,445,264,302,338]`; estimación LPT aproximadamente 151 s por shard. El SHA del plan en scratch es `d387a30bc83ee4fd61ac13e207fd090031398c163b3bcd3dab6efd7d5c30d763`.

## Baseline conocido y gates pendientes

- Antes del merge: colección 2.358; regresión focal 2D/depiction/arquitectura **450 passed, 3 skipped**. Ruff global tenía un F401 histórico (`math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`), fuera del alcance; no se corrigió.
- No se ejecutaron los ocho shards ni la suite pytest monolítica, respetando la restricción vigente por el SIGSEGV histórico de Qt/SIP. Por tanto, no se afirma que la suite completa esté verde.
- No se abrió la GUI integrada ni se generó evidencia de aceptación de la aplicación completa; la revisión manual visual/OpenGL corresponde al propietario en CachyOS.
- Las pruebas de contrato del worker empaquetado pasan, pero no se construyó ni ejecutó un PyInstaller/AppImage real. Live MM por streaming en builds congelados sigue sin aceptación: `_run_worker_stream` arranca `_rdkit_worker.py` directamente, mientras que `_run_worker` usa el protocolo de archivos del ejecutable congelado. Verificar ese camino antes de declarar readiness de release.
- `release/v0.3.0-beta.1-prep`, ramas remotas, tags y canales de publicación permanecen fuera de alcance y sin cambios.

## Decisión

Los conflictos están resueltos y las regresiones focales, de arquitectura, persistencia y planificación CI pasan. Es seguro conservar **un merge commit local** en esta rama aislada; esto no equivale a aceptación completa ni a readiness de Beta 1. Mantener bloqueada cualquier publicación hasta los shards/gates autorizados y la aceptación manual del propietario.
