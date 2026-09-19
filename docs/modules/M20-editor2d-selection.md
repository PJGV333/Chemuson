# M20 — Selección del editor 2D

## Responsabilidad

M20 contiene las implementaciones canónicas y deterministas que comparten los
mixins del canvas: geometría de selección, bounds, hit testing, overlays y
política de clipboard. No posee escena, eventos Qt, comandos undo/redo ni
mutaciones de items.

## Inventario

- `src/chemuson/gui/editor2d/selection/selection_geometry.py`
- `src/chemuson/gui/editor2d/selection/selection_bounds.py`
- `src/chemuson/gui/editor2d/selection/selection_hit_testing.py`
- `src/chemuson/gui/editor2d/selection/selection_overlay.py`
- `src/chemuson/gui/editor2d/selection/selection_clipboard.py`

M20 no representa todo `gui.editor2d`: el directorio padre es un namespace
para futuros módulos hermanos como drawing, text, rendering e interactions.
M20 posee únicamente la selección del editor 2D.

## Compatibilidad

Los imports históricos desde `chemuson.gui.canvas.selection_*` siguen siendo
válidos mediante shims que reexportan desde
`chemuson.gui.editor2d.selection`. Las nuevas dependencias internas del canvas
importan siempre desde `chemuson.gui.editor2d.selection`.

## Verificación

Los tests funcionales y AST de las cinco políticas están registrados en M20 en
`architecture/modules.yml`. El namespace padre `gui/editor2d/__init__.py` no
pertenece a M20. La validación manual de la GUI debe completarse antes de
archivar el OpenSpec de migración.
