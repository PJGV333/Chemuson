# Tasks — edición molecular 3D directa

## 0. Preparación y evidencia

- [x] 0.1 Confirmar `feature/chemuson-3d`, worktree independiente, HEAD, upstream y preservar cambios previos.
- [x] 0.2 Leer los tres OpenSpecs 3D, sus informes de validación y el análisis Avogadro existente.
- [x] 0.3 Capturar baseline sin repetir pytest monolítico SIGSEGV conocido.
- [x] 0.4 Validar este OpenSpec con `openspec validate improve-avogadro-inspired-3d-editing --strict` antes de editar código.
- [x] 0.5 Completar matriz comparativa de fuente Avogadro → implementación actual → defecto → solución; distinguir código fuente, Ghidra e inferencia.

## Hito A — Sincronización 2D/3D

- [ ] 1.1 Añadir regresiones para butano→vacío→etanol, conectividad/layout, debounce, undo/redo y pestañas.
- [x] 1.2 Añadir invalidación inmediata, cancelación y generación debounced solo ante cambio químico real.
- [x] 1.3 Vincular cada resultado/frame al token de documento, firma química y solicitud vigentes.
- [ ] 1.4 Verificar vacío, ráfagas rápidas, cambio de pestaña y que no se publique ningún frame obsoleto.

## Hito B — Hidrógenos editables

- [x] 2.1 Probar picking/resaltado de H derivado y resolver tipo/parent ID sin colisión documental.
- [x] 2.2 Hacer drag de H derivado en escena, preview, cancelación Escape y continuidad de la identidad del padre.
- [ ] 2.3 Promover H a átomo explícito + enlace + coordenada en una transacción; probar undo/redo, save/load y no duplicación en RDKit. (save/load y undo/redo pasan; la prueba RDKit queda saltada porque no está disponible en este entorno.)
- [x] 2.4 Definir y probar rechazo seguro/diagnóstico al borrar un H derivado sin operación de deshidrogenación válida.
- [x] 2.5 Cubrir Live MM usando la coordenada manipulada y restricciones del H promovido/seleccionado.

## Hito C — Herramientas directas

- [x] 3.1 Añadir select/draw/move/bond/measure/fit en barra compacta con ayudas y atajos.
- [ ] 3.2 Añadir construcción por drag desde átomo, preview, enlace/orden, cambio de elemento y validación atómica de valencia. (drag, preview, creación/cambio de enlace y valencia están implementados; cambiar el elemento de un átomo ya existente queda pendiente.)
- [x] 3.3 Añadir movimiento de átomo/fragmento y selección de componente con comandos undoables.
- [ ] 3.4 Añadir handles directos para longitud/ángulo/torsión, mediciones en pantalla y rechazo de enlaces no puente. (mediciones por picks directos están implementadas; faltan handles y selección visual de operaciones restringidas.)
- [ ] 3.5 Probar picking/render, Escape, undo/redo y selección de enlaces/fragmentos sin IDs escritos. (picking, Escape y undo/redo tienen regresiones focalizadas; la selección guiada para todas las restricciones sigue pendiente.)

## Hito D — Iluminación y estilos

- [ ] 4.1 Añadir tests del contrato de espacio de normales/luz y orientación de cámara. (el shader transforma normales/posición a espacio de vista; falta test de orientación de cámara y aceptación GPU.)
- [x] 4.2 Implementar key/fill/ambient/specular equilibrados manteniendo colores por elemento.
- [x] 4.3 Añadir wireframe/space-filling compatibles con el estilo ball-and-stick actual.
- [x] 4.4 Mantener picking coherente con la geometría visible en cada modo; registrar aceptación visual manual pendiente. (contrato cubierto por tests; inspección gráfica no ejecutada.)

## Hito E — Live MM e integración

- [x] 5.1 Probar preview de construcción/drag con Live MM, restricción temporal y cancelación/reemplazo de jobs.
- [x] 5.2 Asegurar un único commit undoable por acción lógica y continuidad de H/estereoquímica.
- [ ] 5.3 Medir latencias monotónicas acotadas del preview y streaming con moléculas representativas. (no medido en esta campaña.)
- [x] 5.4 Ejecutar regresiones de persistencia, cambio de pestañas, staleness, valencia y contratos de arquitectura.

## Documentación y entrega

- [x] 6.1 Crear hoja de ruta separada de superficies científicas avanzadas.
- [x] 6.2 Ejecutar shards focalizados, compileall, Ruff scoped, arquitectura y OpenSpec strict; documentar no ejecución del full suite.
- [x] 6.3 Documentar resultados/limitaciones y prueba manual CachyOS sin declarar validación gráfica si no se hizo.
- [x] 6.4 Commits locales por hito funcional; nunca añadir `PROMPT_LUNA_TEMP.txt`; verificar worktree limpio y no publicar.
