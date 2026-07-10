## Why

El flujo de plantillas de usuario ya existe, pero combina persistencia JSON, normalización, migración de archivos MOL legados, vistas de menú/dock, acciones CRUD e inserción en canvas sin un contrato funcional explícito. Esta propuesta estabiliza ese contrato para reducir regresiones en guardado, importación/exportación, recuperación de bibliotecas existentes e inserción desde la UI.

## What Changes

- Definir el comportamiento esperado de la biblioteca de plantillas de usuario: creación inicial, ruta de almacenamiento, normalización, categorías, plantillas y sincronización de plantillas integradas.
- Especificar el flujo de guardado desde selección o documento completo, incluyendo validaciones de contenido vacío, nombre/categoría y persistencia inmediata.
- Especificar gestión de categorías y plantillas desde dock/menú: crear, renombrar, eliminar y refrescar vistas sin dejar referencias obsoletas.
- Especificar importación/exportación JSON en modos combinar y reemplazar, incluyendo desduplicación, resolución de IDs y conservación de categorías.
- Especificar migración tolerante de plantillas MOL legadas y fallback químico `molblock` -> `smiles` al convertir una plantilla a `MolGraph`.
- Especificar inserción estable desde menú o galería, incluyendo modo de colocación por clic y conexión al átomo cuando aplique.
- Fuera de alcance: rediseño visual del dock, cambios en el formato público de documento Chemuson, nuevo editor gráfico de plantillas, sincronización en nube o refactors globales de `main_window.py`.

## Capabilities

### New Capabilities
- `user-template-flow`: Flujo estable de biblioteca de plantillas de usuario, desde persistencia/importación/exportación hasta gestión UI e inserción en canvas.

### Modified Capabilities
- Ninguna. No hay specs base existentes en `openspec/specs/` para modificar.

## Impact

- Código probablemente afectado: `src/chemuson/gui/template_library.py`, `template_repository.py`, `template_service.py`, `template_store.py`, `template_conversion.py`, `template_chem_adapter.py`, `template_browser_service.py`, `controllers/template_controller.py`, `docks.py`, `main_window.py` y rutas de inserción en canvas.
- Pruebas probablemente afectadas o ampliadas: `tests/test_template_library.py`, `tests/test_template_insertion_mode.py` y nuevas pruebas unitarias/GUI ligeras para controlador, browser service y dock cuando sea viable.
- APIs internas a estabilizar: `TemplateLibrary`, `TemplateRepository`, `TemplateController`, `TemplateBrowserService`, payloads del dock y formato JSON de biblioteca `version/categories/templates`.
- Dependencias: no se esperan dependencias nuevas; se mantiene PySide6/Qt y ChemIO/RDKit como integraciones existentes.
