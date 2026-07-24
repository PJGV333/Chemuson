# Design: Extract canvas selection bounds

## Architecture

### Módulo `selection_bounds.py`

El módulo es una biblioteca de funciones puras de consulta. No importa:

- `PyQt6.QtGui` — ningún tipo de pintado
- `PyQt6.QtWidgets` — ningún item de widget
- `canvas_selection`, `canvas_view` — no acoplamientos al mixin
- `comandos`, `diálogos`, `controladores` — no lógica de aplicación
- Modelos químicos concretos — solo protocolos estructurales
- Renderizado, undo/redo, exportación

### Protocolos estructurales

Ambas funciones operan sobre objetos que cumplen contratos estructurales:

**`model`** para `resolve_selected_atom_ids()`:
- `model.bonds` — mapeable/iterable de `{bond_id: bond}`
- `bond.a1_id`, `bond.a2_id` — atributos de extremo

**`scene`** para `selection_bounds()`:
- `scene.sceneBoundingRect()` — callable que retorna QRectF

**`atom_items`** — `Mapping[int, object]`:
- `.scene()` retorna `scene` o `None`
- `.pen().style()` retorna Qt.PenStyle
- `.brush().style()` retorna Qt.BrushStyle
- `.sceneBoundingRect()` retorna QRectF
- `.label` con `.isVisible()`, `.sceneBoundingRect()`
- `.charge_label` con `.isVisible()`, `.sceneBoundingRect()`

**`bond_items`** — `Mapping[int, object]`:
- `.sceneBoundingRect()` retorna QRectF

**`implicit_h_overlays`** — `Mapping[int, Iterable[tuple[object, object]]]`:
- Cada `(line_item, text_item)` tiene `.scene()` y `.isVisible()` y `.sceneBoundingRect()`

**`graphic_items`** — `Iterable[object]`:
- `.scene()` retorna `scene` o `None`
- `.sceneBoundingRect()` retorna QRectF

### Decoración con `typing.Protocol`

No se definen Protocols formales de typing para mantener compatibilidad con versiones anteriores de Python y evitar imports innecesarios. Se usa duck typing directo.

### Envoltorios

Los tres métodos envoltorio en `CanvasSelectionMixin`:

```python
def _selected_atom_ids_for_transform(self) -> set[int]:
    return resolve_selected_atom_ids(
        selected_atom_ids=self.state.selected_atoms,
        selected_bond_ids=self.state.selected_bonds,
        model=self.model,
    )

def _selected_items_bbox(self) -> Optional[QRectF]:
    from PyQt6.QtCore import Qt
    from chemuson.gui.composite_diagram_item import CompositeDiagramItem
    from chemuson.gui.items import (
        ArrowItem, BracketItem, TextAnnotationItem,
        EnergyDiagramItem, OrbitalAnnotationItem,
        ImageAnnotationItem, WavyAnchorItem,
    )
    from chemuson.gui.plate_items import TLCPlateItem, GelElectrophoresisItem

    atom_ids = self._selected_atom_ids_for_transform()
    graphic_items = [
        item for item in self.scene.selectedItems()
        if isinstance(item, (
            ArrowItem, BracketItem, TextAnnotationItem,
            EnergyDiagramItem, CompositeDiagramItem,
            OrbitalAnnotationItem, ImageAnnotationItem,
            WavyAnchorItem, TLCPlateItem, GelElectrophoresisItem,
        ))
    ]
    return selection_bounds(
        scene=self.scene,
        atom_items=self.atom_items,
        bond_items=self.bond_items,
        implicit_h_overlays=self._implicit_h_overlays,
        atom_ids=atom_ids,
        bond_ids=self.state.selected_bonds,
        graphic_items=graphic_items,
    )

def _targets_bbox(self, **kwargs) -> Optional[QRectF]:
    atom_ids = list(kwargs.get("atom_ids", ()))
    bond_ids = list(kwargs.get("bond_ids", ()))
    graphic_items = []
    for key in ("text_items", "arrow_items", "bracket_items",
                "energy_diagram_items", "semantic_diagram_items",
                "orbital_items", "image_items"):
        graphic_items.extend(kwargs.get(key, ()))
    return selection_bounds(
        scene=self.scene,
        atom_items=self.atom_items,
        bond_items=self.bond_items,
        implicit_h_overlays=self._implicit_h_overlays,
        atom_ids=atom_ids,
        bond_ids=bond_ids,
        graphic_items=graphic_items,
    )
```

Los imports específicos de tipo se mantienen dentro del método para evitar circularidades.

### Conservación del comportamiento

- `extend()` idéntico: ignora `QRectF` inválidos/nulos, usa `united()`
- `extend_atom_bounds()`: misma lógica de pen/brush, label, charge_label, overlays
- El orden de unión no afecta el resultado final de `QRectF.united()`
