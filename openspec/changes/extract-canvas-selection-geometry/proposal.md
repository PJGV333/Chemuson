## Why

Después de comenzar la modularización de la GUI, el monolito
`gui/canvas/canvas_selection.py` todavía posee varias transformaciones y
comparaciones deterministas que no dependen de la escena, el modelo, eventos,
comandos ni estado de interacción. Su ubicación mezcla reglas numéricas
reutilizadas por varios mixins con la coordinación sensible del editor 2D.

## What Changes

- Crear `gui/canvas/selection_geometry.py` como implementación interna de M09.
- Trasladar literalmente delta angular, rotación y escalado de puntos,
  comparaciones tolerantes y normalización de escala/grosor.
- Conservar los nombres privados actuales en `CanvasSelectionMixin` mediante
  aliases estáticos o un wrapper mínimo cuando se necesita el grosor heredado.
- Añadir regresiones numéricas directas y contratos AST de ownership.
- Actualizar el catálogo y la ficha M09.

## Capabilities

### New Capabilities

- `canvas-selection-geometry`: reglas geométricas y normalizaciones puras
  usadas durante transformaciones de selección.

### Modified Capabilities

- `module-catalog`: M09 registra el nuevo módulo interno y sus pruebas.

## Impact

El cambio abre la primera costura interna de `canvas_selection.py` sin alterar
selección, drag, undo/redo, transformaciones de items, ramas, 3D, renderizado,
química, persistencia, eventos ni apariencia.
