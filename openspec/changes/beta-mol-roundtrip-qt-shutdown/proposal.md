# Propuesta — Microcampaña P0: interoperabilidad MOL y cierre Qt

## Problema y evidencia

En aceptación manual de ChemUSON v0.3.0-beta.1 se dibujó `(2S)-2-fluoropropan-1-ol` (`C[C@H](F)CO`, `C3H7FO`). La exportación guardó `prueba_mol_1` sin sufijo; RDKit la leyó, pero ChemUSON y Avogadro no la reabrieron correctamente. Las coordenadas de enlace estaban en escala aproximada del lienzo (40 unidades) y no en la escala MOL convencional (~1.5). Al cerrar AppImage se informó `SIGSEGV`; apareció además una observación de etiquetas C superpuestas tras recuperar autosave.

Inspección inicial de `FileController.save_file` confirma que el filtro MOL selecciona la rama MOL pero no añade `.mol`; el filtro de apertura depende de `.mol`, por lo que el archivo sin extensión no aparece en el filtro usual. Una reproducción fuente controlada guardó la ruta exacta sin sufijo; el bloque se reimportó por contenido a un `MolGraph` de 5 átomos/4 enlaces. La serialización inicial emite distancias MOL ~40. No se infiere que esos síntomas compartan causa.

La aplicación ya tenía coordinación de cierre para QThreads descendientes, pero no desconectaba las señales del canvas activo. Una reproducción fuente endurecida (exportar MOL, aceptar cierre y destruir canvas con `DeferredDelete`) terminó en SIGABRT; la traza identifica `QUndoStack::~QUndoStack → clear → indexChanged → ChemusonWindow._on_active_undo_index_changed` y el slot accedía a un `ChemusonCanvas` ya destruido. Es un defecto productivo distinto de la antigua conexión anónima de test. El crash manual AppImage se informó como SIGSEGV, no se obtuvo una traza de ese binario y por tanto su identidad exacta con la reproducción fuente aún requiere retest.

## Alcance

- Añadir `.mol` en Guardar como de forma idempotente y persistir ruta/estado limpio únicamente tras escritura exitosa.
- Normalizar coordenadas MOL mediante una copia del grafo y un único escalado uniforme, compartido por rutas RDKit/fallback, usando longitud de enlace mediana objetivo cercana a 1.5; no mutar el documento ni reflejarlo.
- Reimportar archivos MOL reales, mantener conectividad, átomos/enlaces, cargas, isótopos y estereoquímica especificada; identificar errores de lectura/escritura en el flujo real.
- Tratar SDF explícitamente como no soportado en esta campaña: no ofrecerlo en filtros ni escribir un MOL simple con extensión `.sdf`; una implementación completa multi-registro requiere un cambio separado.
- Reproducir en subprocesos acotados cierre normal, cierre tras exportar, documento modificado, destrucción de canvas y recuperación de autosave; investigar etiquetas C y conexiones de Qt.
- Modificar el ciclo de vida productivo Qt sólo si se reproduce un fallo atribuible y la causa puede demostrarse. De lo contrario, conservar evidencia y riesgo residual, sin declarar el SIGSEGV resuelto.

## Contratos de compatibilidad

- `.cmsn` conserva su formato, ruta de guardado y flujo actual.
- `.mol` es un único registro CTAB válido. El importador conserva la geometría relativa y escala al tamaño interno del lienzo; exportar no altera las coordenadas fuente.
- La representación del estereoisómero se compara químicamente (incluyendo CIP cuando esté definido), no por texto literal `@`/`@@`.
- No se inventa estereoquímica ausente. Las rutas sin capacidad fiel deben fallar explícitamente.
- SDF multi-registro no se declara compatible ni se mezcla su implementación con esta campaña.
- No cambia `.cmsn`, Clean2D, ChemName, reglas IUPAC ni semántica de renderizado.

## Casos reproducibles

1. Guardar MOL con filtro seleccionado como `molecula` y como `molecula.mol`.
2. Exportar/importar `C[C@H](F)CO` y control R `C[C@@H](F)CO`; abrir tanto MOL propio como MOL generado por RDKit mediante el controlador de archivos.
3. Ejecutar las rutas RDKit y fallback admisibles y medir longitudes directamente del CTAB.
4. Forzar fallos de lectura/escritura y comprobar el error mostrado y que no se actualizan ruta ni estado limpio.
5. Cerrar ventana normal, tras exportar, con cambios pendientes y tras recuperar autosave; repetir cada caso en subproceso con timeout externo.
6. Inspeccionar estructura y etiquetas recuperadas para distinguir duplicación de datos de una representación visible normal.

## Criterios de aceptación

- Los nombres MOL adquieren exactamente un sufijo `.mol`; se conserva el nombre/ruta real; el fallo no marca el documento como guardado.
- El archivo usa escala 2D convencional (mediana de enlace próxima a 1.5), con proporciones/ángulos y quiralidad conservados en ambas rutas que soporten el caso.
- RDKit y ChemUSON reabren los archivos reales. R y S permanecen opuestos; elementos, topología, órdenes, cargas, isótopos y CIP definidos no se pierden silenciosamente.
- El guardado `.cmsn` no regresa. SDF no se anuncia como compatible y recibe rechazo explícito.
- Pruebas Qt sujetas a bloqueo tienen timeout externo; no hay nuevos abortos en las reproducciones controladas. Un no-repro no se informa como causa raíz ni como fix demostrado.
- Se ejecutan arquitectura, OpenSpec estricto, `compileall`, Ruff focal y `git diff --check`; nunca se ejecuta la suite monolítica local.

## Limitaciones

La aceptación manual AppImage pertenece al propietario. Un binario anterior no prueba que el fix esté empaquetado; cualquier nuevo AppImage se generará exclusivamente con Build Preview y se registrarán SHA-256/provenance. Sin una reproducción controlada y una traza que vinculen causa y corrección, el SIGSEGV de la prueba manual queda como riesgo abierto. DATA-11 no se declara completado.
