# Refactorizacion Phase 5

## Hotspots antes de esta fase

- `src/chemuson/gui/main_window.py`: 4104 lineas. Seguía concentrando wiring, workflows de archivo, updater y rich text.
- `src/chemuson/gui/canvas/canvas_input.py`: 6242 lineas. Mezclaba herramientas de enlace, anillos/cadenas, anotaciones, seleccion, teclado y menu contextual.
- `src/chemuson/gui/template_library.py`: 265 lineas. Seguía mezclando fachada publica, persistencia e integracion quimica.
- `src/chemuson/gui/controllers/document_controller.py`: ya usaba contextos explicitos para recientes/descartes, pero conservaba rutas legacy para gestion de pestañas/rutas.
- `src/chemuson/gui/canvas/*` y `src/chemuson/gui/commands/*`: seguian existiendo `from ._shared import *`, con trazabilidad pobre y dependencias ocultas.

## Que se extrajo

### Canvas input

`canvas_input.py` quedo como fachada de compatibilidad y el comportamiento se redistribuyo por responsabilidad:

- `src/chemuson/gui/canvas/canvas_selection_input.py`
- `src/chemuson/gui/canvas/canvas_keyboard.py`
- `src/chemuson/gui/canvas/canvas_context_menu.py`
- `src/chemuson/gui/canvas/canvas_tools_annotations.py`
- `src/chemuson/gui/canvas/canvas_tools_bonding.py`
- `src/chemuson/gui/canvas/canvas_tools_rings_chains.py`

### Main window

Se sacaron tres workflows fuera de la shell principal:

- `src/chemuson/gui/controllers/file_controller.py`
  Responsabilidad: abrir/cargar/guardar/copy-as con `FileWorkflowContext`.
- `src/chemuson/gui/controllers/update_controller.py`
  Responsabilidad: preferencias, chequeo y aplicacion de updates con `UpdateControllerContext`.
- `src/chemuson/gui/rich_text_dialog_service.py`
  Responsabilidad: valor de editor y flujo del dialogo rich text.

`main_window.py` quedo como shell de composicion y wiring, delegando esos flujos a servicios/controladores.

### DocumentController

Se endurecio el contrato explicito:

- `RecentFilesContext` sigue siendo el contrato para recientes.
- `DocumentDiscardContext` sigue siendo el contrato para descarte.
- `DocumentTabsContext` se agrego para la parte de pestañas/rutas.
- Se eliminaron los fallbacks legacy basados en diccionarios internos para `set_canvas_file_path` y `update_tab_title`.

Compatibilidad temporal que queda:

- `set_canvas_file_path(...)` y `update_tab_title(...)` aun aceptan un target legacy solo si expone `_tab_manager`.
- `main_window.py` ya no depende de fallbacks genericos del controlador para recientes/descartes.

### Template library

`TemplateLibrary` paso a ser una fachada mas delgada:

- `src/chemuson/gui/template_repository.py`
  Responsabilidad: load/save/import/export, categorias y CRUD de plantillas.
- `src/chemuson/gui/template_chem_adapter.py`
  Responsabilidad: `graph_from_template` y `add_template_from_graph`.
- `src/chemuson/gui/template_library.py`
  Responsabilidad: API publica delgada y compatibilidad.

## Tamano aproximado antes/despues

| Archivo / area | Antes | Despues | Comentario |
| --- | ---: | ---: | --- |
| `src/chemuson/gui/main_window.py` | 4104 | 3178 | Baja real de responsabilidad; workflows de archivo/update/rich text salen del hub principal |
| `src/chemuson/gui/canvas/canvas_input.py` | 6242 | 19 | Queda solo como fachada de compatibilidad |
| Equivalente repartido de `canvas_input` | 6242 | 6449 | El total se redistribuye en 6 modulos; el mayor queda en 2089 lineas (`canvas_tools_bonding.py`) |
| `src/chemuson/gui/template_library.py` | 265 | 227 | La logica de persistencia y adaptacion quimica sale a componentes dedicados |

Detalle del reparto de `canvas_input`:

- `canvas_selection_input.py`: 1518
- `canvas_keyboard.py`: 548
- `canvas_context_menu.py`: 1166
- `canvas_tools_annotations.py`: 651
- `canvas_tools_bonding.py`: 2089
- `canvas_tools_rings_chains.py`: 477

## `import *` eliminados

Verificacion actual: no quedan coincidencias de `from ._shared import *` en `src/chemuson/gui/canvas` ni en `src/chemuson/gui/commands`.

Archivos limpiados en canvas:

- `src/chemuson/gui/canvas/canvas_render.py`
- `src/chemuson/gui/canvas/canvas_selection.py`
- `src/chemuson/gui/canvas/canvas_structure.py`
- `src/chemuson/gui/canvas/canvas_text.py`
- `src/chemuson/gui/canvas/canvas_view.py`
- `src/chemuson/gui/canvas/__init__.py`

Archivos limpiados en commands:

- `src/chemuson/gui/commands/annotation_commands.py`
- `src/chemuson/gui/commands/atom_commands.py`
- `src/chemuson/gui/commands/bond_commands.py`
- `src/chemuson/gui/commands/diagram_commands.py`
- `src/chemuson/gui/commands/plate_commands.py`
- `src/chemuson/gui/commands/text_commands.py`
- `src/chemuson/gui/commands/transform_commands.py`

Adicionalmente, los comandos dejaron de importar tipos/Qt basicos desde `_shared`; ahora esos imports viven en sus modulos reales.

## Tabla final

| Archivo original | Nuevos archivos creados | Responsabilidad nueva |
| --- | --- | --- |
| `src/chemuson/gui/canvas/canvas_input.py` | `canvas_selection_input.py`, `canvas_keyboard.py`, `canvas_context_menu.py`, `canvas_tools_annotations.py`, `canvas_tools_bonding.py`, `canvas_tools_rings_chains.py` | Separacion de seleccion, teclado, menu contextual, anotaciones, bonding y anillos/cadenas |
| `src/chemuson/gui/main_window.py` | `src/chemuson/gui/controllers/file_controller.py` | Workflow de archivo desacoplado con contexto explicito |
| `src/chemuson/gui/main_window.py` | `src/chemuson/gui/controllers/update_controller.py` | Workflow de update desacoplado con preferencias, chequeo y aplicacion |
| `src/chemuson/gui/main_window.py` | `src/chemuson/gui/rich_text_dialog_service.py` | Flujo del dialogo rich text y helpers asociados |
| `src/chemuson/gui/template_library.py` | `src/chemuson/gui/template_repository.py`, `src/chemuson/gui/template_chem_adapter.py` | Fachada publica ligera sobre persistencia y conversion quimica |

## Tests ejecutados

### No-GUI

Comando:

```bash
PYTHONPATH=src pytest -q tests/test_document_controller.py tests/test_template_library.py tests/test_update_ui_text.py -x
```

Resultado:

- `19 passed in 0.52s`

### GUI relevantes

Comando:

```bash
PYTHONPATH=src pytest -q tests/test_main_window_tabs.py tests/test_text_tool_interaction.py tests/test_template_insertion_mode.py tests/test_double_bond_orientation_toggle.py tests/test_chain_tool_zigzag.py -x
```

Resultado:

- `30 passed in 1.44s`

Chequeos adicionales:

- `python -m compileall src/chemuson/gui/main_window.py src/chemuson/gui/controllers src/chemuson/gui/canvas src/chemuson/gui/template_library.py src/chemuson/gui/template_repository.py src/chemuson/gui/template_chem_adapter.py tests/test_document_controller.py`
- `ruff check --select F401,F811,F821,E741 ...` sobre los archivos tocados por la fase

## Compatibilidad preservada

- `src/chemuson/gui/canvas/canvas_input.py` sigue existiendo como fachada compatible.
- `src/chemuson/gui/canvas/__init__.py` sigue reexportando el API publico del paquete canvas sin volver a usar `import *`.
- `src/chemuson/gui/main_window.py` reexporta `FLATPAK_APP_ID`, `format_no_update_message` y `format_update_disabled_message` para no romper imports existentes.
- `TemplateLibrary` conserva la API publica previa.

## Deuda tecnica restante

- `src/chemuson/gui/canvas/canvas_tools_bonding.py` quedo como el mayor hotspot nuevo. Ya no mezcla todo `canvas_input`, pero sigue siendo candidato claro para una fase 6.
- `src/chemuson/gui/main_window.py` sigue siendo grande porque aun concentra construccion de menus/toolbars/docks y varios handlers de UI de composicion.
- `src/chemuson/gui/canvas/_shared.py` sigue siendo un modulo de compatibilidad grande. Ya no hay `import *`, pero todavia concentra muchas dependencias comunes.
- `src/chemuson/gui/controllers/recovery_controller.py` aun opera contra una interfaz tipo-window. No era el foco de esta fase y quedo pendiente para un desacople posterior.

## Pendiente y por que

- No se partio `canvas_tools_bonding.py` en esta fase para no convertir la deuda en microarchivos artificiales sin una frontera clara todavia.
- No se intento una reduccion cosmetica adicional de `main_window.py`; solo se movieron workflows con independencia real.
- No se eliminaron todas las rutas legacy del ecosistema de tabs/documentos: se mantuvo compatibilidad minima en `DocumentController` para targets que aun exponen `_tab_manager`.
