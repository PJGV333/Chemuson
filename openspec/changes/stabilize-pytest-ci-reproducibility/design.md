# Diseño — confiabilidad de pytest y CI

## Baseline y diagnóstico

Ver `baseline.md`. La rama de campaña parte exactamente del commit solicitado; antes de editar, `HEAD` y `origin/fix/chemio-stereo-roundtrip` coinciden en `f5f82a63c84c56c9dd3ee07e7c7f56e7eb700bef`. La colecta local Python 3.11 contiene 2.128 IDs. No se ejecutó pytest monolítico.

Evidencia local Python 3.11, con un proceso nuevo y límite externo por nodo:

- `test_generate_candidates_attempts_rdkit_for_cyclic_graphs` falla porque el candidato `rdkit_isolated` equivalente se elimina al deduplicar; sí aparecen `simple_aromatic_template` y `rdkit_direct`.
- `test_compchem_controller_generates_async_with_fake_backend`, el diálogo molecular, ambos tests del smoke worker y `test_primary_tabs_have_complete_labels_padding_and_separation` pasan aislados.
- Al colectar una prueba GUI junto con los dos nodos del smoke empaquetado, ambos smoke tests fallan porque la colección ya importó 110 módulos `chemuson.gui`; aislados pasan.
- Una prueba acotada que precarga `QApplication` y una preferencia `resolution_method=ai` reproduce que el fixture existente, que sólo cambia `XDG_CONFIG_HOME`, entrega el valor persistido y falla la expectativa predeterminada. La causa es la ruta NativeFormat ya cacheada por Qt.
- El nodo CompChem instrumentado fuera de pytest observa `worker.finished → QThread.finished → job_finished`; queda por conservar esa telemetría en la regresión. El panel no reproduce localmente el `x=-3` de CI; se medirá geometría estable y se conservará el fallo ante clipping persistente, sin tocar `side_panel.py`.

GitHub run #343 confirma el SHA y que falló el paso pytest; `windows-smoke` y `flatpak-smoke` pasaron. La API accesible sólo expone una anotación genérica de exit 1; la descarga del log de Actions respondió 403. El run #342 corresponde a `chemname/iupac-robustness` en `054f6c7`; el usuario reporta los mismos seis nodos más los tres fallos estéreo corregidos después.

## Decisiones

1. El test Clean2D hará spy del backend `clean2d_isolated`, capturará el candidato producido antes de deduplicar y comprobará que su hash de geometría coincide con el template retenido. No exigirá que dos geometrías equivalentes sobrevivan como candidatos distintos ni omitirá la invocación real.
2. El fixture del Asistente aislará y limpiará explícitamente `QSettings.NativeFormat/UserScope` con `QSettings.setPath`, y restaurará la ruta anterior en `finally`; además comprobará los defaults relevantes. Cambiar sólo una variable de entorno después de crear `QApplication` no se considera aislamiento.
3. Los dos contratos de importación del smoke empaquetado se ejecutarán en procesos Python frescos, que reflejan el límite real del ejecutable y no heredan la tabla `sys.modules` de pytest.
4. El test CompChem observará las tres señales y reportará el orden, jobs activos y señales recibidas si no termina. Se mantiene el límite actual; no se declara resuelta una carrera por aumentar el timeout.
5. El test visual esperará hasta observar geometría idéntica en muestras consecutivas con límite finito, registrará offset/rango del scroll y conservará las aserciones actuales de visibilidad, tamaño de texto, padding, gaps y recorte. Si `Inspector` sigue en negativo al estabilizarse, el test falla y el defecto visual permanece abierto; no se relajan cotas ni se modifica producción sin esa evidencia.
6. La colecta completa se compara con un archivo versionado y ordenado de node IDs. Un plan determinista asigna cada ID exactamente a un único shard; cada runner valida su selección y el XML JUnit. El agregador exige todos los reportes exitosos y completos. Ningún nodo se filtra, se salta ni recibe `continue-on-error`.
7. Ocho shards se ejecutan en GitHub Actions; el supervisor termina cada proceso pytest a los 300 s, conserva log/JUnit/reporte y marca explícitamente timeout o señal fatal. El workflow termina en fallo si falta el plan, un shard o un resultado.
8. Se mantienen los jobs Windows y Flatpak existentes. No se introducen dependencias runtime ni imports nuevos bajo `src/`.

## Riesgos

- El manifiesto de IDs requiere una actualización consciente al añadir, retirar o renombrar tests; esta disciplina es necesaria para que una reducción accidental de colección no pase inadvertida.
- Ocho shards paralelos implican instalaciones repetidas; se reutilizará la caché pip de `setup-python`. La carga se debe medir en CI real y cualquier shard que supere cinco minutos deja la campaña no lista.
- El orden de ejecución dentro de una suite fragmentada difiere del monolítico. Los fallos conocidos se mantienen en shards con resultado y logs; los tests de colección/importación tienen regresiones explícitas.
- El `x=-3` reportado por Actions no se reproduce con Python 3.11 local; el nuevo diagnóstico debe determinar si era layout transitorio o clipping estable.

## Validación

Los pasos de verificación y los resultados finales se anotarán en `validation.md`. La aceptación requiere los seis nodos focalizados, grupos de orden pertinentes, plan de 2.128+ IDs con cobertura exacta, shards bajo 300 s, contratos de workflow, arquitectura, compileall, Ruff focal y un run de Actions verde sobre el SHA publicado. Si push/auth o Actions no están disponibles, el estado final será `NOT READY`.
