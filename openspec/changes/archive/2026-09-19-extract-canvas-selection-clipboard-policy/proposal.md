## Why

La selección del canvas contiene contratos de portapapeles mezclados con
renderizado, creación de items, comandos y undo. Se puede extraer la política
MIME, el codec JSON y decisiones deterministas sin mover pegado ni cambiar
formatos.

## What Changes

- Crear `gui/canvas/selection_clipboard.py`.
- Centralizar MIME strings, detección de formatos pegables, codec de payload,
  política de selección grande y deduplicación/prioridad de enlaces.
- Mantener `copy_to_clipboard`, `paste_from_clipboard`,
  `has_copyable_selection`, `can_paste_from_clipboard` y
  `_paste_selection_payload` en el coordinador.

## Non-goals

No se mueve `_build_selection_graph`, que tiene consumidores en controllers y
canvas_structure. No se mueve `_paste_selection_payload`, no se cambian JSON,
IDs, offsets, imágenes, Molfile/SMILES, comandos, undo, escena ni selección
post-paste.
