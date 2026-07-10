## Context

El flujo actual de plantillas está distribuido en módulos con responsabilidades ya separadas: `TemplateLibrary` como fachada, `TemplateRepository` para persistencia, `template_service` para CRUD puro, `template_store` para almacenamiento/normalización, `TemplateChemAdapter` y `template_conversion` para conversión química, `TemplateBrowserService` para menú/dock/migración MOL, `TemplateController` para acciones de usuario y `PlantillasDock` para interacción visual.

El código ya cubre casos importantes: inicialización de biblioteca por defecto, ruta por plataforma, importación/exportación JSON, fallback de `molblock` a `smiles`, categorías, CRUD de plantillas y colocación de plantillas en canvas. La estabilización debe convertir ese comportamiento en contrato verificable y cerrar huecos de regresión sin reestructurar globalmente la UI.

## Goals / Non-Goals

**Goals:**

- Formalizar requisitos verificables para guardar, gestionar, importar/exportar, migrar, recuperar e insertar plantillas de usuario.
- Mantener compatibilidad del formato JSON actual: `version`, `categories`, `templates`, con entradas que contienen `id`, `name`, `category`, `molblock` y `smiles`.
- Consolidar el comportamiento tolerante ante bibliotecas incompletas, MOL legados, IDs duplicados en importación y fallos de conversión de MOL con fallback a SMILES.
- Asegurar que menú y dock reflejen el estado persistido tras operaciones que mutan la biblioteca.
- Añadir pruebas pequeñas y enfocadas en los módulos existentes antes de tocar arquitectura.

**Non-Goals:**

- No cambiar el formato público de biblioteca salvo reparaciones normalizadas compatibles.
- No mover el flujo a una base de datos ni introducir sincronización remota.
- No rediseñar visualmente `PlantillasDock` ni cambiar el modelo de interacción principal.
- No hacer un refactor amplio de `main_window.py` más allá de ajustes mínimos necesarios para estabilizar el flujo.
- No implementar un editor de plantillas separado.

## Decisions

1. Mantener `TemplateLibrary` como API estable del dominio de plantillas.

   Rationale: ya concentra la fachada usada por la UI y preserva compatibilidad con tests existentes. Alternativa considerada: exponer `TemplateRepository` directamente al controlador; se descarta porque filtraría detalles de persistencia y haría más frágiles los consumidores.

2. Tratar `template_service` y `template_store` como núcleo testeable sin Qt.

   Rationale: las reglas de normalización, agrupación, desduplicación y CRUD pueden verificarse sin UI. Alternativa considerada: probar sólo vía `TemplateController`; se descarta porque mezclaría diálogos Qt con reglas de datos y dificultaría cubrir errores.

3. Mantener persistencia inmediata tras cada mutación de biblioteca.

   Rationale: el comportamiento actual reduce pérdida de datos si la app se cierra después de guardar, renombrar, borrar o importar. Alternativa considerada: persistencia diferida por sesión; se descarta por ser un cambio de experiencia y riesgo de pérdida.

4. Definir importación JSON como operación tolerante y normalizada.

   Rationale: bibliotecas exportadas o editadas manualmente pueden contener categorías faltantes, IDs duplicados o MOL con encabezado reparable. En modo merge se deben agregar sólo entradas nuevas por firma y resolver colisiones de ID; en reemplazo se debe instalar una biblioteca normalizada completa. Alternativa considerada: rechazar archivos parcialmente inválidos; se descarta porque rompe recuperación de datos de usuario.

5. Mantener migración MOL legada como operación best-effort y no bloqueante.

   Rationale: `TemplateBrowserService.migrate_legacy_templates` ya ignora archivos ilegibles o vacíos; estabilizar ese contrato evita que plantillas antiguas impidan abrir la aplicación. Alternativa considerada: mostrar errores por archivo; se descarta por ruido y porque la migración ocurre al arranque.

6. Mantener inserción como contrato entre `TemplateController` y canvas.

   Rationale: la UI selecciona por ID, carga `MolGraph` desde la biblioteca y delega al modo de inserción del canvas. El canvas conserva la responsabilidad de colocar, seleccionar y conectar al átomo. Alternativa considerada: que el dock inserte directamente; se descarta porque acoplaría presentación con modelo químico.

## Risks / Trade-offs

- [Risk] Reforzar normalización podría alterar bibliotecas editadas manualmente. → Mitigation: cubrir con pruebas de datos incompletos, categorías faltantes, entradas sin contenido químico y MOL con header reparable.
- [Risk] El fallback `molblock` -> `smiles` puede ocultar MOL corruptos. → Mitigation: mantener el fallback sólo cuando existe SMILES; si no hay respaldo, propagar error y mostrar mensaje claro en UI.
- [Risk] Refrescar menú/dock tras cada mutación puede recalcular iconos y afectar rendimiento en bibliotecas grandes. → Mitigation: conservar cache de previews y limpiar sólo en refresh; no introducir renderizado asíncrono en este cambio.
- [Risk] Tests GUI pueden ser frágiles por diálogos Qt. → Mitigation: preferir pruebas unitarias de servicios y pruebas de controlador con dependencias fake/monkeypatch; limitar pruebas end-to-end Qt a interacciones críticas.
- [Risk] Migración MOL legada usa rutas históricas bajo `src/chemuson/templates`. → Mitigation: estabilizar como compatibilidad best-effort, no como fuente primaria de plantillas.

## Migration Plan

1. Añadir o reforzar pruebas alrededor del contrato especificado sin cambiar comportamiento público.
2. Implementar ajustes mínimos donde las pruebas revelen inconsistencias entre contrato y código actual.
3. Verificar que bibliotecas existentes cargan, se normalizan y se vuelven a guardar sin perder plantillas válidas.
4. Mantener rollback simple: al no cambiar formato de archivo ni dependencias, revertir los ajustes de código vuelve al comportamiento previo y las bibliotecas JSON siguen siendo legibles.

## Open Questions

- ¿Debe la UI informar cuántos archivos MOL legados fueron omitidos por estar vacíos o ser ilegibles, o se mantiene el silencio actual para no bloquear el arranque?
- ¿Debe la eliminación de una categoría vacía requerir confirmación igual que una categoría con plantillas, o basta con el diálogo actual uniforme?
