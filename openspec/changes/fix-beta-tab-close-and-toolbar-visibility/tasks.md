# Tareas

- [x] 1.1 Confirmar rama/SHA/árbol limpio; leer las especificaciones UI activas y ejecutar baseline antes de editar producción.
- [x] 1.2 Reproducir en Qt el botón inicial y pestaña nueva: icono inicial válido, botón recién creado sin icono, icono presente después de refrescar tema; documentar QSS/padding y límites baseline.
- [x] 2.1 Inicializar en toda pestaña el icono X theme-aware del proveedor existente; mejorar tamaño/accesibilidad/hitbox/QSS sin mezclar el indicador dirty.
- [x] 2.2 Cubrir icono visible antes de hover en ambas themes y en tabs nuevas/existentes/reordenadas; verificar clic cierra sólo el documento objetivo con las confirmaciones actuales.
- [x] 3.1 Añadir `UiPreferences.show_main_toolbar` por defecto true y la clave estable `ui/show_main_toolbar`; actualizar cada escritura UI para conservar las dos visibilidades de barra.
- [x] 3.2 Añadir Mostrar barra de herramientas checkable al inicio de Ver, controlar el `main_toolbar` original y persistir; mantener aparte mostrar menús y acciones auxiliares.
- [x] 3.3 Verificar arranque nuevo, alternancia/reinicio, independencia de menú/panel lateral/rail, QActions/controladores y tema claro/oscuro en tamaños compactos.
- [x] 4.1 Añadir pruebas de atajos Guardar y Undo/Redo, sin modificar el flujo Ctrl+S ya aceptado.
- [x] 4.2 Actualizar sólo los node IDs nuevos, ejecutar regresiones focalizadas, arquitectura, suite completa por shards acotados, Ruff focalizado, OpenSpec estricto y `git diff --check`; investigar y documentar cualquier fallo baseline.
- [ ] 5.1 Crear commit(s) correctivos en `fix/beta-save-shortcut-menu-access`, publicar sólo esa rama por push normal y verificar SHA remoto/CI verde; no integrar a prep/main ni crear tags/releases.
