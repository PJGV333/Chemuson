## Why

La selección del canvas mezcla creación de items, coordinación de escena y
cálculos deterministas de padding, posiciones y hit testing de handles. Extraer
sólo los cálculos abre una costura testeable sin tocar el ciclo de vida visual.

## What Changes

- Crear `gui/canvas/selection_overlay.py` para geometría y consultas de
  overlays/handles.
- Delegar desde `CanvasSelectionMixin` los cálculos puros y conservar allí la
  creación, actualización y visibilidad de `QGraphicsItem`.
- Añadir regresiones numéricas y contratos AST.

## Non-goals

No se crea un mixin, no se mueve dispatch de eventos, no se cambia MRO,
apariencia, escena, drag ni comandos. La prueba Qt manual queda pendiente para
el archivo del OpenSpec.
