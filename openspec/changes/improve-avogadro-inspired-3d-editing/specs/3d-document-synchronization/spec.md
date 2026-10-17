## ADDED Requirements

### Requirement: El visor 3D sigue la identidad química vigente del documento

Cuando cambie la conectividad o identidad química de `MolGraph` en la pestaña activa, ChemUSON SHALL invalidar inmediatamente coordenadas y previews incompatibles, cancelar o desautorizar trabajos RDKit antiguos, mostrar un estado de actualización y generar automáticamente un conformero de la estructura vigente tras debounce. SHALL limpiar la vista cuando el documento quede vacío. El resultado solo SHALL publicarse si corresponden el token de documento, firma química y generación de solicitud actuales. Un cambio de layout 2D que conserve la firma química SHALL preservar una geometría 3D válida y no SHALL reiniciar la generación por cada gesto de dibujo.

#### Scenario: Reemplazo de estructura
- **GIVEN** butano 2D con visor 3D abierto y conformero visible
- **WHEN** se borra el butano y se dibuja etanol
- **THEN** la escena anterior se invalida inmediatamente y el estado indica «Actualizando geometría 3D»
- **AND** tras debounce aparece automáticamente etanol con sus H derivados, sin cerrar/reabrir el dock

#### Scenario: Ráfaga de edición y resultados obsoletos
- **GIVEN** una generación 3D en curso para una estructura
- **WHEN** el usuario cambia rápidamente la conectividad varias veces o cambia de pestaña
- **THEN** se cancelan/desautorizan solicitudes antiguas y solo la estructura/firma/documento vigente puede actualizar la escena

#### Scenario: Vacío, layout y undo/redo
- **GIVEN** geometría 3D válida asociada a un `MolGraph`
- **WHEN** cambia solo la posición 2D de átomos o enlaces
- **THEN** se conserva el conformero 3D sin regenerarlo
- **WHEN** se eliminan todos los átomos
- **THEN** la escena queda vacía y ningún worker obsoleto la repuebla
- **AND** undo/redo restauran/regeneran únicamente el estado correspondiente a cada pestaña