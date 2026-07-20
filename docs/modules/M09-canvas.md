# M09 - Canvas de Edición Molecular

## Descripción General

El `ChemusonCanvas` es el componente central de la interfaz de usuario. Es un lienzo interactivo basado en `QGraphicsView` que permite la visualización, edición y manipulación directa de estructuras moleculares 2D. Su diseño se basa en una arquitectura de mixins para separar responsabilidades de entrada, selección, renderizado y estructura.

## Composición de Mixins (MRO)

La clase `ChemusonCanvas` hereda de los siguientes componentes, en orden de base directa:

```python
class ChemusonCanvas(
    CanvasInputMixin,
    CanvasSelectionMixin,
    CanvasTextMixin,
    CanvasRenderMixin,
    CanvasStructureMixin,
    QGraphicsView,
):
```

**Bases directas:**
1. **`CanvasInputMixin`**: Gestiona la entrada de usuario (mouse/teclado).
   - `CanvasSelectionInputMixin`
   - `CanvasKeyboardMixin`
   - `CanvasContextMenuMixin`
   - `CanvasToolsAnnotationsMixin`
   - `CanvasToolsBondingMixin`
   - `CanvasToolsRingsChainsMixin`
2. **`CanvasSelectionMixin`**: Controla la lógica de selección de átomos, enlaces y otros elementos.
3. **`CanvasTextMixin`**: Gestiona la edición y visualización de texto y anotaciones.
4. **`CanvasRenderMixin`**: Controla el proceso de renderizado visual y exportación.
5. **`CanvasStructureMixin`**: Sincroniza el estado visual con el modelo químico (`MolGraph`).
6. **`QGraphicsView`** (Base de Qt).

**`CanvasToolsBondingMixin`** agrega 5 sub-mixins: `CanvasBondStateMixin`, `CanvasBondModelOpsMixin`, `CanvasBondDragMixin`, `CanvasBondHitTestingMixin`, `CanvasBondGeometryMixin`.

## Inventario de Archivos Relacionados

El canvas y su lógica de soporte se distribuyen en 20 archivos:

- `src/chemuson/gui/canvas/canvas_view.py`: Clase principal `ChemusonCanvas` (1035 líneas).
- `src/chemuson/gui/canvas/canvas_constants.py`: Constantes de interacción y roles.
- `src/chemuson/gui/canvas/canvas_selection.py`: Lógica de selección y hit testing (5154 líneas).
- `src/chemuson/gui/canvas/canvas_text.py`: Gestión de etiquetas de texto.
- `src/chemuson/gui/canvas/canvas_render.py`: Lógica de dibujo y exportación.
- `src/chemuson/gui/canvas/canvas_structure.py`: Sincronización modelo-escena (3875 líneas).
- `src/chemuson/gui/canvas/canvas_input.py`: Manejo de eventos de entrada.
- `src/chemuson/gui/canvas/canvas_commands.py`: Implementación de comandos de edición.
- `src/chemuson/gui/canvas/canvas_items.py`: Implementación de los elementos gráficos.
- `src/chemuson/gui/canvas/canvas_annotations.py`: Manejo de anotaciones y formas.
- `src/chemuson/gui/canvas/canvas_tools_bonding.py`: Sub-mixin para enlaces.
- `src/chemuson/gui/canvas/canvas_bond_drag.py`: Sub-mixin de drag para enlaces.
- `src/chemuson/gui/canvas/canvas_bond_geometry.py`: Geometría de enlaces.
- `src/chemuson/gui/canvas/canvas_bond_hit_testing.py`: Hit testing de enlaces.
- `src/chemuson/gui/canvas/canvas_bond_model_ops.py`: Operaciones de modelo en bonding.
- `src/chemuson/gui/canvas/canvas_bond_state.py`: Estado de bonding.
- `src/chemuson/gui/canvas/canvas_chem_data.py`: Datos químicos del canvas.
- `src/chemuson/gui/canvas/canvas_context_menu.py`: Sub-mixin de menú contextual.
- `src/chemuson/gui/canvas/canvas_keyboard.py`: Sub-mixin de teclado.
- `src/chemuson/gui/canvas/canvas_tools_annotations.py`: Sub-mixin para anotaciones.
- `src/chemuson/gui/canvas/canvas_tools_rings_chains.py`: Sub-mixin para anillos/cadenas.
- `src/chemuson/gui/canvas/canvas_selection_input.py`: Sub-mixin de entrada para selección.
- `src/chemuson/gui/canvas/canvas_selection_item.py`: Lógica de interacción con items seleccionados.

## Estado Compartido y Flujo de Datos

El canvas mantiene un estado interno que coordina la interacción entre la vista y el modelo:
- `model`: Referencia al `MolGraph` actual.
- `scene`: La `QGraphicsScene` que contiene los items.
- `selection`: Conjunto de elementos seleccionados.
- `drawing_style`: Estilos actuales de dibujo (colores, grosores).
- `state`: Variables de control para procesos activos (arrastre, selección, edición de texto).

**Flujo de Eventos**: Los eventos de mouse/teclado son interceptados por `CanvasInputMixin`, que los despacha a los mixins correspondientes (ej. `CanvasSelectionInputMixin` para clicks de selección, `CanvasTextMixin` para edición de etiquetas).

**Sincronización**: `CanvasStructureMixin` asegura que cualquier cambio en el `MolGraph` se refleje en la escena y viceversa, manteniendo la integridad de la representación visual.

## Capacidades de Interacción

- **Selección e Hit Testing**: Soporta selección por átomo, enlace, grupo o área. El hit testing considera la geometría de los átomos y la topología de los enlaces.
- **Undo/Redo**: Todos los cambios en la estructura o el estilo son encapsulados en comandos de `QUndoCommand` y gestionados por el `QUndoStack` que cada `ChemusonCanvas` crea y posee como `self.undo_stack`.
- **Clipboard**: Permite copiar y pegar estructuras (vía SMILES/Molfile) y elementos gráficos (anotaciones, imágenes).
- **Texto y Edición**: Soporte para etiquetas de texto con edición rica y control de formato.
- **Render y Export**: Capacidad de exportar la vista actual en formatos PNG, SVG y PDF.
- **Serialización**: El estado completo del canvas se guarda en archivos `.cmsn`.

## Dependencias Internas

El canvas importa los siguientes módulos de fuera de `canvas/`:

- **Desde `core` (M00)**: `Bond`, `BondStyle`, `BondStereo`, `MolGraph`, `ChemState`, `normalize_opacity`, `normalize_optional_opacity`, `SIMPLE_HYDROGEN_GROUP_LABELS`, `bond_is_structural`, `bond_affects_valence`.
- **Desde `chemio` (M01)**: `molgraph_to_molfile`, `molgraph_to_smiles`, `rdkit_io` (funciones de importación/exportación), `SUPPORTED_ABBREVIATION_LABELS`.
- **Desde `chemname` (M04)**: `MolView`, `iupac_name`, `NameOptions`, `find_rings_simple`, `ring_bonds`.
- **Desde `geometry3d` (M05)**: `Rotation3D`, `conformer_3d_for_graph`, `project_conformer_to_2d`.
- **Desde `gui.items` (M13)**: Elementos gráficos `QGraphicsItem` (referenciados indirectamente a través de canvas_items).
- **Desde `gui.commands` (M11)**: Comandos undo/redo (referenciados indirectamente).
- **Desde `gui.dialogs` (M12)**: Diálogos (referenciados indirectamente).
- **Desde PyQt6**: `QGraphicsView`, `QGraphicsScene`, `QGraphicsItem`, y todas las clases de widgets y eventos.

## Zonas de Riesgo

- **Complejidad de Mixins**: El orden del MRO es crítico; un cambio en la jerarquía de herencia puede romper la propagación de eventos.
- **Interacción con el Modelo**: La sincronización entre la vista gráfica y el modelo químico es una zona de alta sensibilidad a regresiones.
- **Rendimiento de la Escena**: La gestión de miles de elementos gráficos y su actualización puede afectar la fluidez de la UI.
- **Archivos grandes**: `canvas_selection.py` (5154 líneas), `canvas_structure.py` (3875 líneas), `items.py` (5704 líneas).

## Extracciones Futuras

- **Seguras**: Helpers puros de cálculo de hit-testing, lógica de formato de texto y utilidades de exportación de imagen.
- **Peligrosas**: Lógica de despacho de eventos Qt, gestión del ciclo de vida de la `QGraphicsScene` y sincronización profunda del modelo químico.
