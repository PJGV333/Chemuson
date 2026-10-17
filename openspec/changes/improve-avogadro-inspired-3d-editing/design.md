# Design — edición molecular 3D directa

## Evidencia y estado de partida

El análisis estático existente identifica Avogadro 1.2.0 (tag `1.2.0`, commit `aaa67d012828b35314049a9463b77f70415c29c3`) y funciones fuente/Ghidra en `avogadro-1.2.0-analysis.md`. ChemUSON ya tiene picking CPU, cámara perspectiva/ortográfica, `SceneMolecule3D`, edición de coordenadas, construcción por clic, `QUndoStack`, streaming RDKit y Live MM. El defecto A queda confirmado: `_on_active_undo_index_changed()` cancela trabajos y `_reconcile_geometry3d_for_active_document()` oculta geometría incompatible, pero no agenda una generación nueva. Las escenas 3D incluyen IDs string de H derivados, pero `_GeometryDragState` y los mapas de posición del widget admiten solo IDs `int`; por ello el picking del H no puede iniciar un drag válido.

## Decisiones

### Estado y sincronización 2D/3D

- `MolGraph` continúa siendo la fuente; reutilizar `_geometry3d_graph_signature()` como identidad química, excluyendo coordenadas de dibujo y otros atributos de layout.
- Tras cambios de undo/redo, comparar firma anterior y actual. Mantener coordenadas y trabajos RDKit compatibles ante cambios de layout. Cuando cambia la firma, cancelar/desautorizar el trabajo anterior, invalidar preview, limpiar la escena inmediatamente y mostrar «Actualizando geometría 3D».
- Un `QTimer` single-shot por ventana agrupa ráfagas de edición; al vencer, generar una vez desde la pestaña activa si tiene átomos; en vacío, dejar la escena limpia.
- Cada solicitud lleva `document_token`, firma química y generación monotónica. Solo puede publicar frames/finales si coinciden pestaña, firma y generación actuales. No depender solo de `document_revision`, porque el layout 2D cambia esa revisión sin cambiar la química.

### Identidad/edición de H

- Mantener `SceneAtom3D.id` int/string, pero dar semántica explícita al tipo de referencia: documental `int` frente a derivado `str` con relación padre obtenida de `generated_bonds`. Los strings nunca se envían a `MolGraph` ni se truncan/parsean como ID documental.
- El drag opera en coordenadas de escena, actualiza la metadata derivada en preview y mantiene el `generated_bond` padre. Cancelar no modifica `MolGraph` ni el undo stack.
- Confirmar un H derivado ejecuta promoción validada: obtener un nuevo ID documental, crear átomo H + enlace al padre y transferir su XYZ; quitar el H y enlace derivados de metadata. Agrupar topología y geometría en una macro de undo; deshacer restaura la escena derivada y rehacer no duplica el átomo. Los demás H permanecen derivados.
- Live MM recibe las coordenadas manipuladas; el H promovido pasa a formar parte del payload documental y RDKit no debe añadir un segundo H al mismo sitio. Comprobar conteo e identidad después de conformación/optimización.
- No implementar una falsa eliminación visual que reaparezca en el siguiente conformero. Hasta que exista una operación de deshidrogenación químicamente validada, rechazar eliminar H derivado con estado explícito y sin mutación; permitir promoción y eliminación documental solo cuando la validación química confirme la transacción.

### Interacción y coordinación geométrica

- Reutilizar `MolecularViewportWidget`, `MolecularViewportDock`, `pick_scene`, cámara y malla. Añadir modos de herramienta compactos, con gestos distintos de navegación; Ctrl/Shift conservan selección/componente según atajos documentados.
- Construcción: clic/arrastre inicia en átomo o espacio, previsualiza endpoint y enlace y comitea una operación atómica tras validación de valencia. Modificar un enlace seleccionado cambia el orden mediante comando undoable.
- Manipulación de enlace: distinguir longitud, ángulo y torsión con handles/overlays proyectados desde posiciones 3D. Resolver el lado móvil por puente en `MolGraph`; rechazar rotación independiente en anillos/ciclos cerrados. Una acción genera un comando undoable y H derivados reciben la misma transformación.
- Medidas muestran distancia en Å, ángulo/dihedro en grados mientras el gesto está activo. Indicadores permanecen asociados a geometría de mundo y se proyectan con la cámara.
- Controles por IDs quedan como edición precisa accesible, no como camino principal.

### Render y representaciones

- El renderer existente usa `u_mvp` para posición, pero el shader calcula la luz con normal en espacio de objeto. Cambiar contrato a normal de cámara usando matriz normal/view coherente y luz principal fija a la cámara más relleno/ambiente/especular moderada. No alterar colores CPK para ocultar la causa.
- Añadir preferencias de estilo simples: `wireframe` (átomos y enlaces lineales), `ball_and_stick` y `space_filling` (esferas con radios vdW escalados, enlaces ocultos o secundarios según diseño). Picking debe compartir radios/visibilidad del estilo activo.
- No representar cargas, orbitales, densidad electrónica ni potencial energético mediante color falso. La hoja de ruta registra superficies con contrato, backend y licencia apropiados.

### Live MM y undo

- Mantener el worker de proceso NDJSON como dueño de RDKit; no usar QThread ni cálculo pesado en el hilo GUI.
- Durante drag/construcción, debounce/coalesce trabajos, fijar temporalmente los átomos bajo control y comprobar `document_token`, firma química, coordenadas base/generación antes de aplicar frames.
- Los previews son efímeros. Al soltar/confirmar, la última relajación válida puede producir una sola `SetGeometry3DCommand` que se combina con la operación de construcción/edición lógica; no crear una entrada por frame.
- Medir latencia de señal de preview y tiempo entre cambio de gesto y primer frame/final del worker con reloj monotónico en pruebas headless; distinguir tiempo CPU de validez visual/FPS.

## Restricciones de arquitectura y verificación

Sin dependencias runtime nuevas, paquetes nuevos ni cambios de `.cmsn`. Respetar `architecture/modules.yml`. Aplicar TDD por hito y pruebas Qt offscreen solo para lógica/wiring; no iniciar GUI en esta campaña. La suite completa no se repite por el SIGSEGV de baseline. Aceptación visual real sigue manual y el SIGSEGV Qt/SIP no se considera resuelto.
