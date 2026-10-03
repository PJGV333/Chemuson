# Proposal: Panel lateral moderno y refinamiento de la barra de estado (Fase 5)

## Why

Tras las fases 3 y 4, la región derecha sigue formada por siete `QDockWidget`
clásicos independientes. La fase 5 debe integrarlos en una única superficie
lateral de 340 px (ajustable dentro de 336–340 px y del rango permitido
300–340 px), coherente con el spike PyQt6 aprobado, sin reescribir ni duplicar
el contenido funcional que ya implementan esos docks.
También debe sustituir las acciones históricas de visibilidad del menú **Ver**
por navegación explícita a páginas del mismo panel, persistir su estado y pulir
la composición de la barra de estado moderna ya existente.

## What Changes

- Nuevo `src/chemuson/gui/side_panel.py`: fila de cinco tabs principales,
  control de overflow para Espectroscopía y CompChem y un `QStackedWidget` que
  contiene las siete instancias históricas de docks como páginas.
- `shell/assembly.py`: montar el panel junto al lienzo en la región central,
  sin registrar los docks en el área de docking de `QMainWindow`; conservar
  instancias, widgets, señales, handlers y conexiones existentes.
- `main_window_ui_builder.py` y `main_window.py`: reemplazar las acciones
  `toggleViewAction()` de los siete docks por acciones **Ver** que delegan en
  `side_panel.show_page(key)`, y dirigir a Validación las rutas existentes que
  abrían directamente su dock.
- `platform/settings.py`: persistir `ui/side_panel/active_tab` y
  `ui/side_panel/visible`, normalizando valores inválidos a `inspector` y
  `true` respectivamente.
- Tokens/QSS: métricas y estilo de tabs, overflow, páginas embebidas y labels
  semánticas de la barra de estado, usando los tokens light/dark existentes.
  Se conservan `QStatusBar`, `showMessage()`, el indicador IUPAC y la carga.
- Polish visual de `SideTabRow`: `sideW=340`, padding lateral de 3 px y gap de
  3 px entre tabs, conservando fuente de 10 px y underline/acento activo; no
  cambia API, persistencia, docks ni comportamiento funcional.
- Tests de contrato de Fase 5, registro en `architecture/modules.yml` y
  capturas reproducibles bajo `docs/ui-modernization/side-panel-phase-shots/`.

## No Changes

- No se copian ni rediseñan los contenidos funcionales de Inspector,
  Validación, Propiedades, Plantillas, Apariencia, Espectroscopía o CompChem.
- No se altera la lógica química, Clean2D, canvas, geometría, selección,
  persistencia CMSN, orbitales, diagramas, ToolRail, Flyouts ni shortcuts.
- No se implementa empty state ni se inicia la Fase 6.
- No se añaden fórmula, posición del cursor, autosave ni cálculos/polling a la
  barra de estado: no se ampliarán contratos existentes para mostrarlos.
- La acción actual de preferencias del AppBar conserva su handler
  `PreferencesDialog`; la navegación a Apariencia desde la superficie del
  AppBar se hace mediante su hamburguesa → **Ver → Apariencia**. No se
  rediseña AppBar ni se crea una QAction de preferencias paralela.
- Las capturas offscreen son evidencia reproducible, no aprobación visual.
  La aprobación final queda pendiente de la comprobación manual del usuario en
  KDE/Wayland real.

## Compatibility / Rollback

- Los siete atributos históricos de la ventana (`templates_dock`,
  `inspector_dock`, `validation_dock`, `chemical_properties_dock`,
  `spectroscopy_dock`, `compchem_dock`, `appearance_dock`) continúan siendo
  las mismas instancias con sus APIs, widgets y señales originales.
- Cada acción **Ver → página** muestra el panel, activa la página solicitada y
  no muestra ni flota un dock clásico.
- La reversión del cambio restaura el montaje y las acciones previas; no cambia
  modelos, química ni el formato de persistencia de documentos.
