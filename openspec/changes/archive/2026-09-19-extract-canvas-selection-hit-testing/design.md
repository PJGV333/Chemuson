## Boundaries

`selection_hit_testing.py` contiene cuatro funciones de consulta:

- `semantic_diagram_parent(item, composite_type)` recorre `parentItem()` hasta
  encontrar la raíz semántica.
- `get_item_at(...)` conserva el orden de `scene.items(scene_pos)`, promoción
  semántica, descarte de átomos huérfanos, tipos seleccionables y conversión de
  texto hijo a átomo.
- `resolve_click_item(...)` conserva la prioridad de anotaciones, luego el pick
  geométrico de átomo, item de escena y finalmente enlace.
- `selected_annotation_item_at(...)` filtra por escena, visibilidad,
  `sceneBoundingRect`, `contains(mapFromScene(...))` y el mayor `zValue`.

Las clases y el callback de descarte son argumentos explícitos para evitar que
la biblioteca consulte el canvas o importe `canvas_selection_input`.

## Compatibility

`CanvasSelectionInputMixin` mantiene `_semantic_diagram_parent`, `_get_item_at`,
`_resolve_click_item` y `_selected_annotation_item_at`; sus cuerpos delegan a
las funciones nuevas y conservan exactamente sus firmas públicas/privadas.
`_pick_hover_target` permanece en `canvas_bond_hit_testing.py`.

## Risks

El riesgo principal es alterar el orden de prioridad o capturar errores Qt de
forma distinta. No se agregan catches: la única tolerancia existente en la
consulta de anotaciones sigue siendo `RuntimeError`.

## Tests

Las pruebas funcionales cubren promoción, texto hijo, huérfanos, prioridad,
fallbacks y filtros de anotaciones. Las pruebas AST prohíben mutaciones,
handlers y dependencias inversas, y verifican wrappers y consumidor de
`mousePressEvent`.
