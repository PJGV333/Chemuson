# Diseño — estabilización final beta.1

## Baseline y evidencia

La baseline está en `baseline.md`. El SHA inicial esperado coincide con HEAD. El run Actions #39679607676 corresponde al mismo SHA y ejecutó 2.597 IDs en ocho shards. Todos los IDs se ejecutaron y los JUnit coinciden con la asignación: 2.572 passed, 22 skipped, 3 failed, cero errors. Shards 0 y 6 terminaron exit 1, sin timeout ni señal fatal; el resumen falló correctamente.

- Falla 1: `test_compchem_controller_generates_async_with_fake_backend` recibió `worker.finished`, pero antes del límite preexistente de 3 s no observó `QThread.finished` ni `job_finished`; `active_jobs=(1,)` y `pending_results=(1,)`. El auto-connect original encolaba `thread.quit` hacia el event loop GUI. Una regresión controlada bloqueó el event loop principal después de conectar el worker: `worker.finished` sí llegó pero `thread.wait(1000)` falló; tras `QApplication.processEvents()` el thread terminó. Esta es una explicación reproducible del estado observado y no requiere ampliar el timeout. Se cambió sólo esa conexión a `Qt.ConnectionType.DirectConnection` (Qt documenta `QThread.quit()` como thread-safe), conservando `_record_worker_finished` y `thread.finished` en la cola GUI. El nuevo test primero falló con el código anterior y pasa con el cambio. No hay todavía un run remoto en el nuevo SHA; no se afirma que CI esté resuelto.
- La matriz local posterior añadió un bloqueo independiente: shard 6 abortó dos veces con `SIGABRT` en `tests/test_ui_polish.py::test_onboarding_close_releases_overlay`; la pila nativa incluye `QUndoStack::~QUndoStack`/`cleanChanged` durante destrucción QWidget. El archivo `tests/test_ui_polish.py` pasa aislado (38 tests). No se atribuye causalidad a la conexión CompChem ni se modifica fuera de alcance; requiere diagnóstico Qt separado y bloquea declarar verde la matriz.
- Falla 2: `HeldCompChemWorker` en `test_close_waits_for_all_descendant_worker_families_and_suppresses_late_ui` no declara `progress_ready`, aunque el controller conecta esa señal obligatoria. Lanza `AttributeError` en `start_job`, antes de los asertos de cierre. La señal de producción se añadió después que el doble; se trata de un doble obsoleto, no de un fallo de producción. Corregir sólo el doble, manteniendo la señal tipada y las verificaciones de cierre/ausencia de UI tardía.
- Falla 3: `test_existing_async_worker_shutdown_contract_in_isolated_pytest_process` vuelve a ejecutar el archivo y propaga el mismo `AttributeError`; es duplicación del fallo 2, no una causa adicional.
- `QOpenGLWidget is not supported on this platform` aparece como aviso bajo `QT_QPA_PLATFORM=offscreen`; no es el origen de esos tres failures. El run no terminó por SIGSEGV ni timeout. Las señales de cierre no quedaron probadas por el test fallido antes de la corrección.

## Decisiones

1. Actualizar el doble de prueba `HeldCompChemWorker` para exponer la señal real `progress_ready`; no añadir fallback dinámico en producción.
2. Desacoplar únicamente la salida del QThread CompChem del event loop GUI mediante `DirectConnection`; dejar intactos el timeout de 3 s, el relay de resultado y las aserciones de `job_finished`/`thread.finished`.
3. El test de regresión debe bloquear deliberadamente el event loop principal entre inicio y cierre, observar `worker.finished` en el hilo emisor y confirmar que QThread termina; luego procesa callbacks GUI y valida que el job se retire.
4. El test de cierre llega a sus aserciones, libera workers controlados y comprueba ausencia de UI tardía. El wrapper subprocess es una segunda verificación del mismo contrato.
5. Reutilizar `tools/ci_pytest_shards.py`, manifiesto de IDs, perfil de duración, `test.yml`, `build-preview.yml`, `release.yml` y validadores existentes. Shards locales seriales, timeout individual máximo 300 s.
6. Conservar JUnit/log/timings/report/plan; no añadir binarios descargados al repo. No editar `KNOWN_BASELINE_EXCEPTIONS` para ocultar resultados.
7. La aceptación/distribución requiere SHA/run/artifact ID/checksum y tipo de evidencia. No hay push, tag, Release, PR ni canal público sin autorización.
8. Un aborto Qt no previsto (aunque sea del mismo grupo histórico `QUndoStack`) es bloqueo independiente: registrar evidencia y solicitar tratamiento/OpenSpec específico, no arreglar oportunísticamente UI ajena a la falla CompChem.

## Riesgos y límites

- El runner GitHub original usó Python 3.11.17/Qt 6.11.2/RDKit 2026.9.1. Se creó una venv temporal Python 3.11.15 con los mismos paquetes PyQt6 6.11.0/Qt 6.11.2/RDKit 2026.9.1; el wrapper de shutdown allí hereda un `PYTHONPATH` a `chem/lib/python3.14/site-packages` hardcodeado por el test y falla por extensión ABI incompatible. El mismo wrapper pasa con `chem/bin/python` Python 3.14.8; se registra como limitación del runner local, no cambio de fuente.
- Una matriz local de 2.598 IDs, ejecutada serialmente en Python 3.14.8, tuvo siete shards verdes y shard 6 abortó a 144.659 s (308/373 timings; 65 sin ejecutar). El mismo shard, repetido una vez en otro directorio, abortó igual en el mismo test. El verifier rechaza la matriz. El monolítico basal de 8 min había expirado sin resumen.
- Sin un nuevo run Actions/preview en el SHA estabilizado no se puede declarar CI remoto verde ni validar cuatro binarios nuevos. Sin Windows/OpenGL interactivo no se afirma aceptación de esos entornos. El aborto local bloquea la decisión final hasta diagnóstico/decisión del owner.
- Los históricos que llaman incorrectamente AppImage al ejecutable se corrigen con nota explícita de documento histórico; no se reescriben resultados como si fueran Type 2.

## Addendum — alcance autorizado de cierre público (2026-10-23)

El propietario autorizó de forma expresa integrar la rama de estabilización en
la rama prep por fast-forward-only, push normal, preparar el cuerpo público de
la beta y generar el preview oficial desde la SHA final. Se conserva el freeze
de funcionalidades. La autorización no elimina gates ni equivale a aceptación
de paquetes: tag, GitHub Release y despliegue Flatpak siguen condicionados a
aceptación P0/P1 exacta y protección/configuración externa de publicación.
La reconciliación entre `AGENTS.md` y la política de pruebas se documenta como
propuesta separada, sin ejecutar la suite masiva ni cambiar sus controles.
