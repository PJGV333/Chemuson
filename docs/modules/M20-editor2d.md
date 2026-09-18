# M20 — Políticas del editor 2D

## Responsabilidad

M20 contiene las implementaciones canónicas y deterministas que comparten los
mixins del canvas: geometría de selección, bounds, hit testing, overlays y
política de clipboard. No posee escena, eventos Qt, comandos undo/redo ni
mutaciones de items.

## Inventario

- `src/chemuson/gui/editor2d/selection_geometry.py`
- `src/chemuson/gui/editor2d/selection_bounds.py`
- `src/chemuson/gui/editor2d/selection_hit_testing.py`
- `src/chemuson/gui/editor2d/selection_overlay.py`
- `src/chemuson/gui/editor2d/selection_clipboard.py`

El paquete no depende de otros módulos ChemUSON. PyQt6 sólo aparece donde los
contratos geométricos necesitan `QPointF`/`QRectF`.

## Compatibilidad

Los imports históricos desde `chemuson.gui.canvas.selection_*` siguen siendo
válidos mediante shims que reexportan M20. Las nuevas dependencias internas del
canvas importan siempre desde `chemuson.gui.editor2d` para que la propiedad
canónica permanezca explícita.

## Verificación

Los tests funcionales y AST de las cinco políticas están registrados en M20 en
`architecture/modules.yml`. La validación manual de la GUI debe completarse
antes de archivar el OpenSpec de migración.
