# Propuesta — opción de visibilidad de la barra superior

## Problema

La ventana ofrece dos QAction distintas en `Ver`: `Mostrar barra de menús`, conectada a la visibilidad de `QMenuBar`, y `Mostrar barra de herramientas`, el `toggleViewAction()` de la `QToolBar` antigua. La hamburguesa ya abre un agregador que contiene los mismos seis `QMenu` originales. La segunda acción muestra Nuevo/Abrir/Guardar en una fila adicional, mientras la barra de menús que el propietario quiere controlar usa otra etiqueta.

## Alcance

- Reutilizar una sola acción de Ver con el texto solicitado «Mostrar barra de herramientas», pero conectada exclusivamente a la visibilidad de la `QMenuBar`.
- Conservar la preferencia existente `ui/show_menu_bar`, el default visible y el callback de persistencia; una configuración nueva no necesita migración adicional.
- Mantener la `QToolBar` histórica oculta, incluso si `ui/show_main_toolbar` quedó en `true` en una configuración anterior, y no ofrecer su toggle como opción de usuario.
- Usar el mismo `QMenu` de Ver y la misma QAction desde la hamburguesa para recuperar la fila cuando esté oculta.
- Añadir pruebas Qt acotadas de alternancia, persistencia, toolbar oculta y preservación de AppBar/acciones/temas/tamaños.

## Fuera de alcance

No rediseñar la cabecera, no cambiar acciones de archivo, pestañas, atajos, panel lateral o comportamiento químico. No cambiar Clean2D, ChemIO, ChemName, otras preferencias, `release/v0.3.0-beta.1-prep`, `main`, tags, releases ni paquetes.
