# Tasks: Extract canvas selection bounds

## Tarea 1: Crear módulo selection_bounds.py

- [x] Crear `src/chemuson/gui/canvas/selection_bounds.py`
- [x] Implementar `resolve_selected_atom_ids()`
- [x] Implementar `selection_bounds()`
- [x] Validar compilación: `python -m compileall -q src/chemuson/gui/canvas/selection_bounds.py`

## Tarea 2: Refactorizar canvas_selection.py

- [x] Reemplazar `_selected_atom_ids_for_transform()` con delegación
- [x] Reemplazar `_selected_items_bbox()` con delegación
- [x] Reemplazar `_targets_bbox()` con delegación
- [x] Eliminar duplicación de `extend` y `extend_atom_bounds`
- [x] Validar: `python -m compileall -q src/chemuson/gui/canvas/canvas_selection.py`

## Tarea 3: Actualizar catálogo arquitectónico

- [x] Añadir `selection_bounds` a `internal_api` de M09
- [x] Añadir `tests/test_canvas_selection_bounds.py` a `tests` de M09
- [x] Añadir `tests/architecture/test_canvas_selection_bounds.py` a `tests` de M09
- [x] Actualizar conteo de líneas de canvas_selection.py en notas

## Tarea 4: Pruebas funcionales

- [x] Crear `tests/test_canvas_selection_bounds.py`
- [x] Cubrir `resolve_selected_atom_ids`: átomos, enlaces, enlace inexistente, no mutación, vacío
- [x] Cubrir `selection_bounds`: vacío, QRectF inválidos, unión numérica, átomo visible, átomo oculto, label, charge_label, H implícito, bond item, graphic item, scene removal, combinación

## Tarea 5: Pruebas arquitectónicas AST

- [x] Crear `tests/architecture/test_canvas_selection_bounds.py`
- [x] Verificar existencia del módulo
- [x] Verificar funciones propietarias
- [x] Verificar prohibición de imports
- [x] Verificar ausencia de llamadas de mutación
- [x] Verificar wrappers conservados
- [x] Verificar catálogo actualizado

## Tarea 6: Validación final

- [x] `python -m compileall -q src tests tools packaging`
- [x] `pytest -q tests/test_canvas_selection_bounds.py tests/architecture/test_canvas_selection_bounds.py`
- [x] `ruff check src/chemuson/gui/canvas/selection_bounds.py src/chemuson/gui/canvas/canvas_selection.py tests/test_canvas_selection_bounds.py tests/architecture/test_canvas_selection_bounds.py`
- [x] `openspec validate --all --strict`
