## Why

`CanvasSelectionInputMixin` mezcla consultas de escena, promoción semántica,
selección de anotaciones y la política que combina picking molecular con items
gráficos. Esa política puede probarse sin mover eventos Qt ni el picking
geométrico de átomos/enlaces.

## What Changes

- Crear `gui/canvas/selection_hit_testing.py` como biblioteca consultiva.
- Extraer promoción de padres semánticos, resolución de item bajo un punto,
  prioridad del click y anotación seleccionada superior.
- Mantener los cuatro métodos privados existentes como wrappers compatibles.
- Añadir regresiones funcionales y contratos AST.

## Scope / Non-goals

No se mueven `mousePressEvent`, eventos, `_pick_hover_target`, comandos,
controllers, dialogs, escena ni modelos químicos. No se cambia la prioridad
observable ni la selección del usuario.

## Impact

M09 conserva el canvas interactivo y su MRO. El módulo nuevo sólo consulta la
escena/items y recibe callbacks/clases necesarias del wrapper.
