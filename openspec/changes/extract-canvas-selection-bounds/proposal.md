# extract-canvas-selection-bounds

## Problem

`CanvasSelectionMixin._selected_items_bbox()` y `CanvasSelectionMixin._targets_bbox()` contienen una lógica duplicada de consulta de límites visuales de selección. Ambas funciones implementan:

- `extend()`: unión de `QRectF` válidos.
- `extend_atom_bounds()`: recopilación de límites de átomos, etiquetas, cargas e hidrógenos implícitos.

Además, `_selected_atom_ids_for_transform()` resuelve IDs de átomos desde enlaces seleccionados de forma independiente.

Esta duplicación dificulta la prueba unitaria y la evolución del comportamiento visual.

## Proposal

Extraer dos funciones puras de consulta hacia un nuevo módulo interno:

1. `resolve_selected_atom_ids()` — resuelve el conjunto de IDs de átomos desde selección de átomos y enlaces.
2. `selection_bounds()` — calcula el bounding box visual a partir de atom IDs, bond IDs y grupos de graphic items.

Los tres wrappers privados de `CanvasSelectionMixin` se convertirán en delegadores mínimos:

- `_selected_atom_ids_for_transform()` → `resolve_selected_atom_ids()`
- `_selected_items_bbox()` → `selection_bounds()` con items seleccionados de la escena
- `_targets_bbox()` → `selection_bounds()` con los grupos explícitos recibidos

## Design

### Módulo `selection_bounds.py`

- Ubicado en `src/chemuson/gui/canvas/`
- Solo importa `collections.abc`, `PyQt6.QtCore.QRectF`, `PyQt6.QtCore.Qt`
- Trabaja por duck typing (protocols estructurales), sin imports de QtGui, QtWidgets ni items concretos
- Sin efectos secundarios

### Contratos

**`resolve_selected_atom_ids()`**

```python
def resolve_selected_atom_ids(
    *,
    selected_atom_ids: Iterable[int],
    selected_bond_ids: Iterable[int],
    model: object,
) -> set[int]:
```

- Copia IDs de átomos seleccionados
- Añade extremos de enlaces existentes en el modelo
- Ignora enlaces inexistentes

**`selection_bounds()`**

```python
def selection_bounds(
    *,
    scene: object,
    atom_items: Mapping[int, object],
    bond_items: Mapping[int, object],
    implicit_h_overlays: Mapping[int, Iterable[tuple[object, object]]],
    atom_ids: Iterable[int] = (),
    bond_ids: Iterable[int] = (),
    graphic_items: Iterable[object] = (),
) -> QRectF | None:
```

- Para cada átomo: verifica item existe y pertenece a la escena; incluye cuerpo (si pen != NoPen o brush != NoBrush), label visible, charge_label visible, overlays visibles de H implícitos
- Para enlaces: incluye `sceneBoundingRect()` del bond item si existe
- Para graphic_items: incluye solo los que pertenecen a la escena

### Pruebas

- Funcionales: cobertura completa de ambos contratos
- Arquitectónicas AST: existencia de funciones, prohibición de imports externos, integridad de wrappers

## Tasks

1. Crear `src/chemuson/gui/canvas/selection_bounds.py`
2. Refactorizar `canvas_selection.py` para delegar en el nuevo módulo
3. Actualizar `architecture/modules.yml`
4. Crear `tests/test_canvas_selection_bounds.py`
5. Crear `tests/architecture/test_canvas_selection_bounds.py`
6. Validar con `openspec validate --all --strict`

## Impact

- `canvas_selection.py` reduce su complejidad ciclomática
- El nuevo módulo es 100% testeable sin Qt widgets
- Sin cambios observables en la GUI
