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

## Workers de segundo plano

`src/chemuson/gui/background_workers.py` contiene únicamente los workers Qt
privados para descriptores RDKit y Name→Structure. Sus imports de backend son
diferidos y el módulo no importa la ventana.

`ChemusonWindow` conserva el lifecycle de `QThread`, los registros de jobs,
cancelación, diálogos de progreso y callbacks que actualizan la interfaz. Los
workers no forman parte de la API pública de M08.

## Geometría de integración Clean2D

`src/chemuson/gui/clean2d_geometry.py` posee las transformaciones puras usadas
para centrar, reescalar, alinear y completar coordenadas al integrar resultados
Clean2D o proyecciones 3D. El módulo no importa Qt, canvas, ventana ni backends
de dominio.

`ChemusonWindow` conserva aliases estáticos privados para no cambiar todavía
los puntos de llamada de controllers y pruebas. Estas funciones no forman
parte de la API pública de M08.
