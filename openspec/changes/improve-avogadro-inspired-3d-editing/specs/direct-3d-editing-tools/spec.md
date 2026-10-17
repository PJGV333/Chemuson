## ADDED Requirements

### Requirement: Las herramientas 3D convierten gestos directos en operaciones químicas undoables

El viewport SHALL exponer herramientas reconocibles de selección, dibujo, manipulación, edición centrada en enlace y medición, además de centrar/ajustar vista. Los gestos de edición SHALL distinguirse de la navegación de cámara y SHALL incluir preview/ayuda de contexto. La UI SHALL permitir trabajar sin escribir IDs de átomo; los campos numéricos permanecerán disponibles para precisión. Cambios topológicos SHALL validar valencia antes del commit, no dejar mutaciones parciales y registrarse como una operación lógica.

#### Scenario: Crear átomo y enlace por gesto
- **GIVEN** la herramienta de dibujo y un elemento/orden de enlace elegidos
- **WHEN** el usuario arrastra desde un átomo a espacio vacío
- **THEN** aparece un preview de endpoint/enlace y al soltar se crean átomo y enlace validados en una sola acción undoable
- **AND** Escape cancela sin alterar documento ni historial

#### Scenario: Seleccionar y manipular fragmento
- **GIVEN** una escena con componentes conectados
- **WHEN** el usuario selecciona varios átomos o un componente y arrastra
- **THEN** se resalta el conjunto afectado, se mantiene la conectividad y el gesto produce un único cambio geométrico undoable

#### Scenario: Edición centrada en enlace
- **GIVEN** un enlace acíclico seleccionado con átomos vecinos suficientes
- **WHEN** se arrastra un handle para longitud, ángulo o torsión
- **THEN** la parte móvil se transforma respecto del eje correcto y un overlay legible muestra Å/grados durante el gesto
- **AND** la operación aparece como una acción undoable

#### Scenario: Enlace de ciclo no rotatable
- **GIVEN** un enlace cuyo corte no separa la estructura en dos componentes
- **WHEN** se solicita torsión independiente alrededor de ese enlace
- **THEN** la torsión se rechaza con explicación estable y no deja una geometría parcialmente deformada

#### Scenario: Cambiar elemento u orden de enlace
- **GIVEN** un átomo/enlace seleccionado y un elemento/orden objetivo
- **WHEN** se ejecuta la herramienta contextual
- **THEN** se valida la química antes del commit, se sincroniza `MolGraph` y la escena 3D, y undo/redo restaura topología y geometría coherentes