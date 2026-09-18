## Classification

- `_ensure_selection_overlay`: scene mutation; permanece en el mixin.
- `_apply_selection_overlay_bbox`: scene mutation con cálculo de layout;
  permanece como orquestador y delega padding/posiciones.
- `_update_selection_overlay` y `_update_drag_selection_overlay`: event/view
  coordination; permanecen.
- `_hit_selection_handle`, `_hit_selection_move_handle` y
  `_hit_selection_scale_handle`: coordinación de estado; permanecen.
- `_hit_handle_item`: query wrapper; delega distancia/radio.
- `_handle_item_distance_sq`, `_handle_item_hit_radius` y
  `_selection_handle_hit_kind`: pure/query; se extraen.

## Extracted contracts

`selection_overlay.py` posee `padded_selection_bbox`, `offset_scene_point`,
`selection_handle_scene_positions`, `handle_item_distance_sq`,
`handle_item_hit_radius` y `selection_handle_hit_kind`. Las fórmulas, offsets,
radii, orden scale/rotate/move, desempate estable y manejo de excepciones se
trasladan literalmente.

## Tests

Se cubren padding, tres posiciones, drag translation, distancia en view,
radio mínimo, prioridad/desempate, handles ausentes/invisibles y RuntimeError.
El módulo no importa widgets/items ni llama mutaciones.
