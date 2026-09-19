## Ownership

`selection_clipboard.py` sólo contiene valores y política determinista:

- MIME constants con los strings existentes.
- `mime_has_pasteable_format(mime)` con el orden/formatos actuales.
- `encode_selection_payload` y `decode_selection_payload` para el JSON UTF-8.
- `is_large_clipboard_structure(atom_count, bond_count, ...)` sobre counts.
- `bond_copy_priority` y `unique_bonds_for_copy` para la política de copia.

El coordinador conserva toda interacción con `QApplication.clipboard()`,
QMimeData, render/export, comandos, undo stack, escena y payload paste.

## Compatibility

Los wrappers privados permanecen y las firmas públicas no cambian. El codec
rechaza datos inválidos devolviendo `None`; el wrapper conserva el fallback
existente a los otros formatos. No se normalizan ni renombrarán keys.

## Tests

Se cubren todos los MIME strings, codec roundtrip/malformed data, límites de
selección grande, deduplicación y prioridad. AST prohíbe commands, controllers,
dialogs, undo/redo y mutación de escena en el módulo extraído.
