# Baseline: extract-canvas-selection-bounds

## Estado inicial

- Rama: `architecture/phase11-extract-canvas-selection-bounds`
- HEAD: `6aa3a5d3` Normalize archived OpenSpec formatting
- Contiene fase 10: `402569a7` Extract canvas selection geometry
- Árbol: limpio

## canvas_selection.py

- Líneas: 5128
- Métodos objetivo: `_selected_items_bbox` (línea 1976), `_targets_bbox` (línea 4549), `_selected_atom_ids_for_transform` (línea 2060)

## Baseline de compilación

```
python -m compileall -q src tests tools packaging
```

Resultado: limpio (sin errores)

## Baseline de Ruff (selección preexistente)

Los errores E402 en `packaging/release/set_version.py` y `src/chemuson/chemio/cml_io.py` son preexistentes y ajenos a este cambio.
