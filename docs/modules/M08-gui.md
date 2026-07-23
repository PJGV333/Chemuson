# M08 — GUI y application shell

## Responsabilidad

M08 posee la interfaz PyQt6, incluida la ventana pública
`chemuson.gui.main_window.ChemusonWindow`, el editor 2D, acciones, docks,
controllers y composición visual de alto nivel.

## Application shell

`src/chemuson/gui/shell/` ensambla una instancia ya creada de
`ChemusonWindow`. Su responsabilidad se limita a:

- estado y colaboradores de la ventana;
- pestañas y canvas inicial;
- docks y toolbars;
- menús y conexiones globales;
- barra de estado, tema inicial y chequeo diferido de actualizaciones.

El shell no crea `QApplication`, no controla el proceso y no sustituye al
composition root de M19. Tampoco posee handlers de química, persistencia,
documentos, templates, validación o edición del canvas.

## Límites

- M19 crea el proceso gráfico y la instancia de `ChemusonWindow`.
- M08 ensambla y coordina la ventana.
- `ChemusonWindow` conserva handlers y lifecycle propio de QMainWindow.
- M09 conserva la coordinación específica del canvas/editor 2D.

La extracción es estructural: constructores, parentage Qt, conexiones,
visibilidad inicial y orden de arranque visual se conservan.
