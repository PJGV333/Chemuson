# Propuesta — cierre visible y barra de herramientas para la beta

## Problema

La aceptación visual de la rama `fix/beta-save-shortcut-menu-access` encontró que algunas pestañas se crean con el botón Cerrar sin icono hasta el siguiente cambio de tema. Además, `main_toolbar` continúa ocultándose incondicionalmente aunque existe una opción distinta para incluir acciones auxiliares en esa barra.

## Alcance

- Inicializar el icono de cierre en cada pestaña al construirla; garantizar contraste, legibilidad, área de interacción cómoda, tooltip y accesibilidad, manteniendo separado el punto dirty y reutilizando el flujo de cierre existente.
- Exponer en Ver una única acción **Mostrar barra de herramientas** que gobierna la visibilidad del `QToolBar main_toolbar`, visible por defecto y persistida en `ui/show_main_toolbar`.
- Agrupar al inicio de Ver las acciones existentes de barra de menús y la nueva de barra de herramientas; preservar el acceso hamburguesa y la opción independiente copiar/pegar/zoom.
- Cubrir temas, documentos nuevos y reordenados, guardado/Undo/Redo, persistencia, panel lateral y tamaños compactos con tests Qt acotados.

## Fuera de alcance

No cambiar Ctrl+S, handlers/controladores, persistencia CMSN, `CanvasTabManager`, flujos de cierre, Undo/Redo, rail de dibujo, panel lateral, iconos de rail, química, Clean2D, ChemName, empaquetado, canales, workflow de release, beta, `main` o tags. No añadir dependencias ni duplicar QActions de comandos o toolbars.
