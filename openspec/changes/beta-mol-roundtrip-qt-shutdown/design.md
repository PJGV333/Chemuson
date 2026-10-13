# Diseño — interoperabilidad MOL y diagnóstico de cierre

## Hechos y diagnóstico a confirmar

- `FileController.save_file` distingue el filtro MOL, pero escribe a la ruta elegida literalmente; por eso el filtro no agrega `.mol`. `ChemusonWindow` sólo ofrece MOL/SDF por sufijo en el filtro de apertura.
- `molgraph_to_rdkit_with_map` conserva coordenadas internas del lienzo sin conversión. El fallback MOL escribe esos mismos valores. La importación luego aplica `_scale_to_default(..., 40)`, por lo que una exportación normalizada debe partir de una copia y no alterar el lienzo.
- `molfile_to_molgraph` conserva cuñas CTAB; comprobar si el backend asigna también metadatos CIP en el MolGraph importado. El cambio sólo complementará extracción de metadatos cuando exista estereoquímica MOL explícita.
- `ChemusonWindow.closeEvent` hace confirmaciones antes del cierre asíncrono, detiene autosaves/timers y espera QThreads descendientes. Sin embargo, el cierre aprobado dejaba conectadas señales del canvas activo a la ventana.
- Reproducción source controlada en subproceso: insertar fixture, exportar MOL, cerrar ventana y enviar `QEvent.Type.DeferredDelete` al canvas. Antes del fix, el proceso salió `-6` (SIGABRT). GDB mostró `QUndoStack::~QUndoStack → QUndoStack::clear → QUndoStack::indexChanged → PyQtSlotProxy → ChemusonWindow._on_active_undo_index_changed`; el traceback Python llega a `self.canvas.current_opacity_percent()` y `RuntimeError: wrapped C/C++ object of type ChemusonCanvas has been deleted`. El slot de la ventana seguía conectado mientras se destruía el QUndoStack hijo del canvas.
- Causa raíz confirmada para esta reproducción: el cierre de ventana no llamaba a `_disconnect_canvas_signals` para el canvas activo antes de destruir canvas/QUndoStack. La solución desconecta el canvas activo sólo después de aprobar todas las confirmaciones de cierre; un close cancelado conserva las conexiones.
- El antecedente de CI tenía una conexión anónima de test retenida, y era un defecto diferente. El reporte manual AppImage indica SIGSEGV, mientras la reproducción source causa SIGABRT PyQt ante excepción en slot; la semejanza del camino de teardown no demuestra que sean idénticos. El AppImage nuevo debe retestearse.

## Decisiones

1. **Normalizador MOL compartido.** Calcular longitudes de enlaces estructurales no nulas, usar la mediana robusta, escalar uniformemente al objetivo 1.5 y aplicar el mismo resultado a una copia aislada del MolGraph antes de elegir RDKit, worker o fallback. El factor positivo no refleja la estructura. Si no hay longitudes utilizables, mantener coordenadas. No llamar este normalizador en importación, Clean2D ni SMILES.
2. **Exportación sin mutación.** Copiar el grafo completo (incluidas anotaciones) y cambiar sólo coordenadas de la copia. Mantener writer RDKit y fallback detrás de una entrada común; no duplicar fórmulas de escala.
3. **Save As por contrato de filtro.** Interpretar explícitamente el filtro elegido. Para MOL agregar `.mol` sólo cuando falta; conservar sufijo existente sin duplicarlo. Publicar ruta, limpiar undo stack, cancelar snapshot autosave y añadir recientes sólo después de escritura correcta.
4. **SDF fail-closed.** El diálogo no lista SDF. Las rutas explícitas `.sdf` se rechazan con explicación porque el escritor actual produce un único MOL sin `$$$$` y el lector no implementa contrato multi-registro. No añadir soporte SDF parcial.
5. **Stereo MOL.** Probar el fixture principal `(2S)-2-fluoropropan-1-ol` y el enantiómero R contra RDKit con identidad molecular canónica, fórmula, topología, órdenes, cargas, isótopos y CIP. Abrir archivos propios y RDKit por FileController. Mantener sin asignar los centros no especificados; no degradar en fallback.
6. **Qt investigación primero.** Ejecutar casos de cierre y autosave en procesos separados con timeout externo. Incluir errores C++/Python y código de salida; para SIGSEGV, capturar faulthandler/core/backtrace si se reproduce con gdb disponible. No introducir sleeps arbitrarios, `terminate`, catches globales o cambios cosméticos de Qt. Si no hay fallo reproducible ni causa demostrable, sólo pruebas/documentación y riesgo residual explícito.
7. **Autoguardado mínimo.** Verificar un snapshot con el contrato existente, reabrirlo a través de `RecoveryController` y comparar grafo y número/visibilidad de AtomItems. No cambiar persistencia `.cmsn` ni el esquema autosave.
8. **Arquitectura.** Cambios en módulos existentes M01 ChemIO y M09 GUI/file workflow. No se añaden dependencias ni paquetes; Clean2D y ChemName quedan intactos.

## Validación

Ver `baseline.md` y `validation.md`. Todos los pytest se lanzan con timeout externo inferior a 10 minutos; pruebas Qt con riesgo de bloqueo se encapsulan en subprocess. Ejecutar shards focales, arquitectura y controles de persistencia/autosave, no suite monolítica local. Package acceptance sólo con Build Preview sobre SHA nuevo; propietario realiza el retest manual.
