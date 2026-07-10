## 1. Baseline y Cobertura de Persistencia

- [x] 1.1 Ejecutar pruebas actuales de plantillas para establecer baseline: `PYTHONPATH=src pytest -q tests/test_template_library.py tests/test_template_insertion_mode.py`.
- [x] 1.2 Añadir pruebas de normalización para bibliotecas parciales: categorías duplicadas/blancas, categorías faltantes en plantillas, entradas sin contenido químico y categoría de usuario garantizada.
- [x] 1.3 Añadir pruebas de persistencia inmediata para crear categoría, renombrar categoría, borrar categoría, añadir plantilla, renombrar plantilla y borrar plantilla.
- [x] 1.4 Ajustar `template_store.py`, `template_service.py` o `template_repository.py` sólo si las pruebas muestran inconsistencias con el contrato.

## 2. Importación, Exportación y Migración

- [x] 2.1 Añadir pruebas de importación merge que cubran categorías nuevas, duplicados por firma e IDs duplicados resueltos sin sobrescribir plantillas existentes.
- [x] 2.2 Añadir pruebas de importación replace que verifiquen reemplazo completo con datos normalizados y biblioteca actual sin mezcla residual.
- [x] 2.3 Añadir prueba de error de importación para asegurar que un archivo inválido no muta la biblioteca actual.
- [x] 2.4 Añadir pruebas de migración MOL legada en `TemplateBrowserService`: importar archivo nuevo, omitir duplicado, ignorar archivo ilegible/vacío y no propagar errores de discovery.
- [x] 2.5 Ajustar import/export o migración legacy sólo donde sea necesario para cumplir los escenarios especificados.

## 3. Recuperación Química e Inserción

- [x] 3.1 Añadir pruebas de conversión `molblock` primero, fallback a `smiles` y error controlado cuando no existe contenido químico usable.
- [x] 3.2 Reforzar pruebas de inserción por ID en `TemplateController` con contextos fake: iniciar modo por clic, colocación inmediata con posición conocida/fallback y error de carga sin entrar en modo.
- [x] 3.3 Mantener o ampliar prueba de canvas para conexión al átomo y verificar que se crea exactamente un enlace de conexión.
- [x] 3.4 Ajustar `template_conversion.py`, `template_chem_adapter.py`, `controllers/template_controller.py` o canvas sólo si la cobertura descubre divergencias.

## 4. Menú, Dock y Refresh de Vistas

- [x] 4.1 Añadir pruebas de `TemplateBrowserService.refresh_templates_menu` para categorías con plantillas, categorías vacías y acciones fijas de guardar/categoría/import/export.
- [x] 4.2 Añadir pruebas de `TemplateBrowserService.refresh_template_views` para payloads con IDs estables, limpieza de cache y tolerancia a fallo de preview.
- [x] 4.3 Añadir pruebas ligeras de `PlantillasDock.set_templates` para placeholders, payload de categoría, payload de plantilla y señales de selección/contexto cuando sea viable.
- [x] 4.4 Ajustar browser service o dock sólo con cambios mínimos necesarios para que menú y dock reflejen la biblioteca tras mutaciones.

## 5. Validación Final

- [x] 5.1 Ejecutar el conjunto específico de plantillas: `PYTHONPATH=src pytest -q tests/test_template_library.py tests/test_template_insertion_mode.py` y las pruebas nuevas añadidas.
- [x] 5.2 Ejecutar pruebas relacionadas con UI/controladores si se tocaron rutas Qt: `PYTHONPATH=src pytest -q tests/test_main_window_tabs.py tests/test_update_ui_text.py` cuando apliquen.
- [x] 5.3 Ejecutar `openspec validate estabilizar-flujo-plantillas-usuario --strict` y corregir cualquier problema de especificación.
- [x] 5.4 Revisar diff final para confirmar que no hubo refactor global ni cambios fuera del flujo de plantillas de usuario.
