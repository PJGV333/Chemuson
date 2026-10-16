# Tareas — visibilidad de la fila de menús

## 1. Inspección y baseline
- [x] 1.1 Leer los OpenSpec activos y trazar QActions, QMenuBar, QToolBar, hamburguesa y QSettings antes del cambio.
- [x] 1.2 Registrar el estado y comandos baseline en `baseline.md`, con timeout máximo local de 10 minutos.

## 2. Corrección mínima
- [x] 2.1 Reusar la QAction persistida de visibilidad de QMenuBar con el texto «Mostrar barra de herramientas» y retirar del menú Ver el toggle de QToolBar.
- [x] 2.2 Mantener la QToolBar antigua oculta aun con preferencias obsoletas; conservar visible por defecto la AppBar compacta y el menú hamburguesa.
- [x] 2.3 Añadir regresiones Qt para default/estado checkable, solo-menú toggle, restauración desde hamburguesa, reinicio/QSettings, QToolBar antigua, temas/tamaños y acciones existentes.

## 3. Verificación y entrega
- [x] 3.1 Tests GUI focalizados, preferencias, guardado/cierre, `tests/architecture`, Ruff focal, OpenSpec estricto y `git diff --check`; cada pytest <10 min.
- [x] 3.2 Validar el set completo (2235 IDs) con los ocho shards CI bounded; reconciliar node IDs y actualizar el manifiesto sólo por los cambios de tests revisados.
- [x] 3.3 Commit pequeño y push normal sólo a `fix/beta-final-ui-clean2d`; Actions verde, SHA remoto verificado y árbol limpio.
- [x] 3.4 Dejar aceptación manual del propietario pendiente; no integrar en prep ni publicar tags/releases.

## 4. Seguimiento — popup compacto y accesible
- [x] 4.1 Reorganizar las opciones de Ver en submenús reutilizando las QActions/conexiones existentes; la acción de recuperación queda primera y marcada.
- [x] 4.2 Retirar la opción y wiring obsoletos de Copiar/Zoom hacia la QToolBar; preservar las acciones funcionales de Copiar/Zoom y sus rutas/atajos.
- [x] 4.3 Añadir pruebas de mouse QTest, límites de `QScreen.availableGeometry()` para menú y submenús, temas claro/oscuro, tamaños solicitados, restauración por hamburguesa y persistencia tras reinicio; guardar evidencia visual.
- [ ] 4.4 Ejecutar validaciones focalizadas acotadas, arquitectura, Ruff, compileall, OpenSpec estricto y CI por shards; publicar en esta misma rama y reportar SHA/resultados.
- [ ] 4.5 Esperar la revisión manual del propietario; no integrar en beta prep, no tag/release y no declarar beta lista.
