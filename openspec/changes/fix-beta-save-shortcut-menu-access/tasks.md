# Tareas

- [x] 1.1 Registrar SHAs/estado y crear `fix/beta-save-shortcut-menu-access` desde `release/v0.3.0-beta.1-prep`.
- [x] 1.2 Leer el OpenSpec beta activo, contratos de shell/pestañas/preferencias/persistencia y los siete adaptadores solicitados; capturar baseline antes de modificar producto.
- [x] 1.3 Reproducir Ctrl+S, acción Guardar y comprobar la hipótesis de asociación al QMenuBar oculto con el contenido real de `.cmsn`.
- [x] 2.1 Hacer que el QAction existente de Guardar y sus atajos funcionen con menú compacto; conservar `Ctrl+G` y el controlador común.
- [x] 2.2 Añadir regresiones de guardado Qt: Ctrl+S/menú, documento nuevo/existente, activo/inactivo, contenido guardado, dirty tras cancelar/fallar, sin confundir exportación MOL.
- [x] 3.1 Añadir la preferencia Mostrar barra de menús en Ver y persistir/aplicar el modo compacto/clásico mediante settings existentes.
- [x] 3.2 Añadir regresiones de interacción para menú visible/compacto, misma QAction/QMenu, hamburguesa, tema claro/oscuro y tamaños distintos.
- [x] 4.1 Verificar botón X visible/accionable y Save/Discard/Cancel, pestaña limpia y pestaña dirty no activa; comprobar que cancelación/error no cierra ni pierde datos.
- [x] 4.2 Ejecutar pruebas focalizadas (timeout externo ≤600 s), ocho shards completos y acotados, arquitectura, OpenSpec estricto, Ruff focalizado y `git diff --check`; registrar resultados y actualizar el manifiesto CI de node IDs sólo con los IDs revisados.
- [ ] 4.3 Crear commits pequeños, publicar sólo la rama correctiva por push normal si hay autenticación y comprobar SHA/CI. No integrar a beta ni crear tag/release.