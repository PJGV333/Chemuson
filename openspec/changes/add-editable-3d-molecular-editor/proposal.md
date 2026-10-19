# Proposal — editor molecular 3D y mecánica molecular en vivo

## Por qué

Las fases 1–2 ya entregan un `QOpenGLWidget` moderno con cámara, picking y escenas 3D, pero el dock sólo inspecciona una conformación efímera guardada en estado global de ventana. Para editar una molécula en 3D hace falta que las coordenadas pertenezcan al documento, que los gestos sean operaciones reversibles sobre ese mismo documento y que los resultados de mecánica molecular lleguen realmente mientras el cálculo progresa.

El backend RDKit ya se ejecuta aislado y produce frames, pero `rdkit_backend.optimize_iter()` llama a una operación bloqueante que devuelve todos los frames al finalizar el subproceso. Emitir después esos frames en un `QThread` no constituye optimización en vivo. Este cambio añade un protocolo de frames en flujo, control de cancelación y descarte por revisión.

La inspección estática de Avogadro 1.2.0 orienta los contratos de interacción y el manejo de restricciones; no se porta código, algoritmos fuente ni API GPL.

## Alcance

- Incorporar un modo de edición 3D explícito, apagado por defecto, sobre el visor existente.
- Guardar un `CoordinateSet3D` por `ChemusonCanvas`, separado de cámara, selección y coordenadas 2D.
- Editar coordenadas de átomos y fragmentos, distancias, ángulos y dihedros; rechazar rotaciones de enlace que no separan el grafo por ciclos.
- Construir átomos y enlaces compatibles desde 3D usando `MolGraph` y comandos undo/redo existentes.
- Reflejar en el visor cambios estructurales hechos desde 2D, conservando IDs y coordenadas supervivientes.
- Persistir coordenadas 3D explícitas como campo opcional y versionado dentro de `canvas` en `.cmsn`, compatible con documentos anteriores.
- Añadir minimización MM en vivo progresiva con RDKit aislado, MMFF94/MMFF94s y alternativa UFF, restricciones de arrastre, revisión de resultados, debounce y cancelación.
- Probar etano, butano —incluida una torsión seguida de relajación— y benceno; guardar evidencia visual y una secuencia real de coordenadas/energías.

## Fuera de alcance

- Integrar o copiar código, dependencias, gestos o recursos GPL de Avogadro.
- Qt4, Python 2, PyOpenGL o dependencias runtime nuevas.
- Alterar las coordenadas 2D al mover u optimizar una molécula en 3D.
- Hacer persistentes cámara, selección, tareas activas, hidrógenos generados o frames intermedios.
- Presentar minimización como dinámica molecular: esta versión implementa relajación incremental; MD con integración temporal/termostato queda fuera.
- Refactor general del apagado Qt/SIP o del SIGSEGV baseline; documentar y no atribuirlo a esta función sin evidencia diferencial.
- Cambios en `main`, otros worktrees, ramas, push, merge o PR.

## Módulos probables afectados

- M05 `geometry3d`: transformaciones, validación de coordenadas y estado 3D puro.
- M04 `chemio`: flujo NDJSON/proceso supervisado para el worker RDKit existente.
- M08 `gui`: controles del dock, interacción, comandos de geometría, ciclo por documento y sincronización 2D/3D.
- `chemio/persistence.py` y `gui/canvas/canvas_structure.py`: campo opcional de coordenadas en el payload ya propiedad del canvas.
- `architecture/modules.yml`: sólo actualizar ownership/imports si la implementación cambia las dependencias catalogadas; no crear un paquete superior nuevo.

## Criterio de éxito

El usuario puede construir o modificar una molécula desde el dock, mover átomos y fragmentos, editar una torsión, deshacer/rehacer dejando el grafo y la geometría consistentes, guardar/reabrir la conformación y activar la minimización MM para observar coordenadas y energía evolucionar en frames emitidos antes de terminar el cálculo. Una prueba automatizada demuestra que una edición nueva cancela o invalida el trabajo anterior y que una respuesta atrasada no sobrescribe el documento.
