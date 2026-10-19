## ADDED Requirements

### Requirement: Herramienta activa y navegación se distinguen

El viewport SHALL reflejar una sola fuente de verdad para `edit_mode`/`edit_tool` y controles visibles. Ningún botón SHALL aparecer activo para edición mientras su lógica esté deshabilitada. Elegir una herramienta SHALL activar edición. El usuario SHALL poder diferenciar órbita, pan, zoom, selección, movimiento y edición geométrica con gestos/modificadores accesibles.

#### Scenario: Herramienta elegida activa edición
- **GIVEN** el viewport en modo navegación
- **WHEN** se elige Select, Build, Move, Bond, Measure o Geometry
- **THEN** el estado de edición se activa y la herramienta visible coincide con el widget
- **AND** Escape cancela la preview, sin mutar modelo ni historial

### Requirement: Selección visible y construcción validada

Picking SHALL resaltar átomos documentales/derivados y enlaces; Shift SHALL mantener una selección múltiple ordenada y Ctrl SHALL seleccionar componente. La construcción SHALL permitir crear en vacío, extender desde un átomo y enlazar átomos existentes; cambio de elemento y orden SHALL usar comandos undoables. Antes de commit SHALL validar estructura/valencia y rechazar gestos inválidos atómicamente.

#### Scenario: Seleccionar y editar sin IDs escritos
- **GIVEN** etano, etanol o butano cargado en 3D
- **WHEN** el usuario elige y toca átomos/enlaces/H derivados con herramientas del viewport
- **THEN** la selección visible corresponde al target pickeado y el usuario puede construir, cambiar elemento/enlace o mover sin introducir IDs manualmente
- **AND** undo/redo restaura el grafo y la geometría

### Requirement: Handles de longitud, ángulo y dihedro

La selección guiada de 2, 3 o 4 átomos SHALL habilitar, respectivamente, un gesto de longitud, ángulo o dihedro. El fragmento móvil SHALL derivarse de conectividad tras cortar el enlace pivote; si un ciclo impide separar los lados, la operación SHALL rechazarse con estado y sin mutación. Desplazamiento de profundidad SHALL estar disponible y los H derivados del fragmento SHALL acompañar la transformación.

#### Scenario: Ajustar un dihedro y cancelar
- **GIVEN** una cadena abierta seleccionada en el orden A–B–C–D
- **WHEN** el usuario arrastra el handle de torsión
- **THEN** cambia el diedro con enlace B–C conservado y se presenta la medición en grados
- **WHEN** presiona Escape
- **THEN** escena, `MolGraph` e historial vuelven al estado anterior
- **WHEN** suelta válidamente
- **THEN** se crea una sola entrada undoable y Live MM puede relajarse después

#### Scenario: Un enlace de ciclo no tiene lado móvil único
- **GIVEN** un enlace pivote dentro de un ciclo cerrado
- **WHEN** se intenta girar su lado con Geometry
- **THEN** la operación se rechaza explícitamente sin deformar ni confirmar la molécula
