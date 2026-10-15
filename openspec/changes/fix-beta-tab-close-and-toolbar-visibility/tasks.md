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
- [x] 5.1 Crear commit(s) correctivos en `fix/beta-save-shortcut-menu-access`, publicar sólo esa rama por push normal y verificar SHA remoto/CI verde; no integrar a prep/main ni crear tags/releases.

## Ronda 2 — shell compacto solicitado tras la segunda aceptación visual

- [x] 6.1 Baseline de la ronda: reproducir y medir en Qt la altura real de `QMenuBar` (43 px) y `main_toolbar` (65 px, nueve botones) y la geometría de la pestaña limpia vs. dirty (X desplazada de `x=2` a `x=14` y desbordando el rect de 28 px).
- [x] 6.2 A3: bloque lateral de tamaño fijo con el espacio del punto dirty reservado permanentemente; pintura condicional por propiedad dinámica; QSS sin `min-width/min-height` en el botón de cierre. Verificar geometría idéntica (rect de pestaña, tamaño de bloque, tamaño y posición de la X, holgura al borde derecho) en limpio/sucio/limpio, con temas y tamaños de ventana.
- [x] 6.3 A1: barra de menús tradicional visible por defecto y compacta reutilizando los seis `QMenu`/`QAction` originales; opción `Ver → Mostrar barra de menús` y acceso hamburguesa/Alt intactos cuando está oculta. Medir altura ≤ 36 px sin solapamiento con la app bar.
- [x] 6.4 A1: default nuevo migrable una sola vez (`ui/shell_defaults_version`) sobre INIs previos con `ui/show_menu_bar=false` persistido, respetando después la elección del usuario.
- [x] 6.5 A2: `main_toolbar` reducido a Nuevo/Abrir/Guardar con las mismas QActions; retirar undo/redo, rotaciones, volteos, Clean2D y SMILES de esa barra sin eliminarlos del programa; compactar con QSS dedicada (`objectName("mainToolBar")`, icono 18 px) y botones uniformes.
- [x] 6.6 Actualizar los tests que codificaban el modo compacto por defecto y la barra de nueve botones; aislar las preferencias del shell en tests (`isolated_shell_settings`) para que los defaults sean deterministas.
- [x] 6.7 Añadir pruebas focalizadas de las tres correcciones (`tests/test_beta_compact_shell_and_tabs.py`): menús por defecto/compactos/migración/hamburguesa, barra de tres botones uniforme y theme-stable, geometría estable de pestañas y cierre correcto.
- [x] 6.8 Volver a ejecutar baseline comparado, arquitectura, Ruff focalizado, OpenSpec estricto, `git diff --check` y la suite completa por shards acotados; documentar el fallo intermitente CompChem preexistente sin relajar contratos.
- [x] 6.9 Publicar sólo la rama `fix/beta-save-shortcut-menu-access` y registrar SHA remoto y resultado de CI; sin merge a `release/v0.3.0-beta.1-prep`/`main`, sin tags ni releases.
