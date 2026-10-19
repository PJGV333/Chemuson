# Design — corrección de UI 3D y generación compleja

## Diagnóstico de baseline

- `MolecularViewportDock` es `QDockWidget`; el tema global se aplica a la ventana principal, pero el dock no recibe un tema/paleta propia cuando pasa a ventana superior. Los widgets de control no tienen colores explícitos de superficie/texto. Se conectarán `topLevelChanged` y cambios de tema a un QSS generado desde tokens y limitado al objectName del dock.
- El widget declara `edit_mode=False`, mientras que el botón `move` aparece checked. Elegir otra herramienta sí activa edición, pero la discrepancia inicial presenta una acción desactivada como lista. El estado visual debe derivarse conjuntamente del modo y la herramienta; los gestos de navegación tendrán modificadores documentados.
- `_on_compchem_job_finished()` presenta errores en `CompChemDock`, pero el `MolecularViewportDock` conserva una escena vacía y no muestra el motivo ni una acción de reintento. La respuesta del worker pierde parte del diagnóstico estructurado al cruzar `rdkit_safe`/`Conformer3DResult`.
- `_handle_graph_conformer3d_mode()` sanitiza, añade H, ejecuta ETKDGv3 con coordenadas aleatorias y repite con los mismos parámetros/semilla; luego si falla MMFF/UFF aún puede devolver geometría con `embed-only` sin mensaje claro. Se reemplazará por una secuencia limitada, semillas diferenciadas, ruta macro cíclica condicionada por tamaño de anillo, validación y diagnóstico por etapa.
- La ruta `conformer_3d_isolated()` es síncrona y no consume el evento de cancelación del worker de aplicación; una cancelación de QThread no termina el subproceso RDKit. La generación en streaming NDJSON usará el terminador existente, con eventos de etapa reales y timeout monotónico.
- La referencia externa auténtica de tetrandrina reproduce la generación local correctamente en baseline; no reproduce el archivo 2D que el usuario dibujó. El documento original no está en el repositorio.

## Decisiones

### Dock/tema/contexto

- La hoja de estilo específica se genera con tokens existentes y se asigna únicamente a `MolecularViewportDock` y sus controles. No se cambia `QApplication`, el tema KDE ni QSS global.
- `topLevelChanged` reaplica paleta/stylesheet después del cambio de ventana; los cambios claros/oscuros llaman el mismo setter. La escena/cámara/selección/tool state residen en el widget y deben sobrevivir a hide/show/context teardown. La limpieza GPU no debe limpiar el modelo de escena.
- Cuando `edit_mode` es falso se muestra estado de navegación y no se marca un botón de edición como activo. Al elegir una herramienta, ésta activa modo edición; apagarlo cancela solo previews no confirmados y deja la herramienta como preferencia no activa.

### Edición directa

- Picking conserva IDs tipados; la selección múltiple mantiene orden guiado para seleccionar 2/3/4 átomos. Rendering/highlight debe reflejar toda la selección documental y derivada, y los enlaces.
- Construcción, cambio de elemento/enlace y promoción de H pasan por comandos existentes, validación de `MolGraph` antes de mutar y undo/redo. No se permite que un gesto inválido deje topología parcial.
- Los handles geométricos son gestos guiados por selección: 2 átomos para longitud, 3 para ángulo y 4 para dihedro. El lado móvil se obtiene separando el enlace pivote y haciendo búsqueda de conectividad; si permanece conectado por un ciclo, se rechaza sin mutación. La longitud traduce el fragmento en el eje del enlace; ángulo/dihedro rotan el fragmento alrededor del eje geométrico. H derivados ligados al fragmento reciben la misma transformación.
- Escape restaura preview; release válido produce un único commit de coordenadas. Cámara: botón derecho pan, rueda zoom, arrastre izquierdo en modo navegación orbita; Alt+arrastre en herramienta Move controla profundidad; Shift amplía selección y Ctrl selecciona componente.

### Conformeros/diagnóstico

- Instrumentar serialización payload, conteos, sanitización, H, asignación/validación estereoquímica, embedding, minimización, conversión a CoordinateSet/Scene, aceptación de snapshot y upload/render intentado. Usar `time.perf_counter()` por etapa.
- Errores mínimos estables: `invalid_graph`, `embed_failed`, `timeout`, `cancelled`, `unsupported_stereochemistry`, `missing_forcefield_parameters`, `geometry_invalid`, `worker_failed`, `stale_result`, `scene_failed`; incluir etapa, detalle, conteos, retorno RDKit/semillas y tiempos cuando estén disponibles.
- Intentos finitos y reproducibles: ETKDGv3 convencional; parámetros de macro ciclo solo si un anillo los requiere; semillas distintas limitadas; coordenadas aleatorias como último intento cuando proceda. Sanitizar y validar finitud/geometría/estereoquímica antes de aceptar.
- Guardar copia de la conformación embebida validada antes de minimización. Si el campo de fuerza no existe, falla o produce geometría que no supera validación, devolver la conformación no optimizada marcada claramente, no un viewport vacío ni una geometría 2D.
- Los eventos de progreso corresponden a transiciones reales de etapa; el proceso aislado se termina al cancelar o exceder timeout. El controlador y la ventana descartan por snapshot/tab/firma y publican resultado, fallo y progreso solo al documento actual.

## Restricciones

Sin paquetes ni imports nuevos de terceros. La implementación permanece en módulos y capas ya existentes, respeta `architecture/modules.yml` y no altera `persistence.py`, Clean2D o ChemName. Pruebas Qt offscreen validan wiring/estado, no legibilidad real ni recuperación GPU en un compositor KDE. No iniciar GUI; aceptación gráfica queda manual.
