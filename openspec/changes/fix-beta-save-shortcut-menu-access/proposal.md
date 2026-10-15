# Propuesta — microcorrección de guardado y acceso a menús

## Problema

En la aceptación de `v0.3.0-beta.1`, la misma acción `Archivo → Guardar` funciona al activarse, pero Ctrl+S no alcanza la `QAction` mientras el `QMenuBar` está oculto en el modo compacto. Ctrl+G está asociado a esa misma acción y debe conservarse. Se solicita además el modo opcional de menús clásicos siempre visibles. Los cierres de pestañas dirty requieren regresión de interacción para Save/Discard/Cancel.

La reproducción y causa raíz están registradas en `baseline.md`; no se presume el origen antes de observarlo.

## Alcance

- Asociar el `QAction` existente de guardar directamente con la ventana para que sus atajos se procesen también con menús compactos ocultos; mantener el mismo `_on_file_save` y `FileController`.
- Cubrir `.cmsn` nuevo/existente, sólo documento activo, contenido escrito, estado dirty tras cancel/error y cierre explícito de pestañas dirty/limpias/no activas.
- Añadir en Ver la acción verificable **Mostrar barra de menús**, reutilizando las seis `QMenu`/`QAction` y persistiendo el modo por el sistema de preferencias UI existente.
- Mantener el botón hamburguesa en modo compacto; no cambiar la barra de herramientas de dibujo, la Command Palette ni la disposición de documentos.

## Fuera de alcance

No modernizar otras áreas de interfaz ni alterar química, ChemName, Clean2D, ChemIO, estereoquímica, el formato/persistencia `.cmsn`, canales beta/stable, tags, releases o workflows de publicación. No merge a beta.