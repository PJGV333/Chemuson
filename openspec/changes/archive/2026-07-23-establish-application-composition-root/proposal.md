## Why

El entry point `chemuson.__main__:main` conserva correctamente el parsing CLI
y el import diferido de la GUI, pero la creación de `QApplication`, la
instalación del crash reporter, la construcción de `ChemusonWindow` y el ciclo
de vida completo todavía viven al final de `gui/main_window.py`. Esto mezcla la
composición de la aplicación con la definición de una ventana de más de 2600
líneas y deja a M19 como un único archivo que no representa su responsabilidad
catalogada.

La fase de límites arquitectónicos ya terminó con cero excepciones y cero
ciclos. El siguiente incremento debe establecer un composition root real y
pequeño antes de extraer el application shell o dividir el editor 2D.

## What Changes

- Crear `src/chemuson/app/` como paquete de composición perteneciente a M19.
- Trasladar allí, sin alterar su orden ni efectos observables, el arranque de
  Qt, la configuración de metadata, la creación de la ventana, la recuperación
  de autosaves, el manejo de fallos y la terminación del proceso.
- Mantener `src/chemuson/__main__.py` como parser y dispatcher CLI mínimo con
  import diferido del composition root.
- Retirar `run_app()` de `gui/main_window.py`; ese módulo conservará la ventana,
  no el ciclo de vida de la aplicación.
- Actualizar M19, las dependencias reales del catálogo y su documentación.
- Añadir pruebas de importación/dispatch que demuestren que `--version` y la
  carga del entry point no importan PyQt6 ni la GUI.
- Corregir la descripción obsoleta de deuda temporal en
  `docs/architecture.md`.

## Capabilities

### New Capabilities

- `application-composition-root`: Composición y ciclo de vida explícitos de la
  aplicación de escritorio, separados de la definición de la ventana.

### Modified Capabilities

- `module-catalog`: M19 pasa a poseer `src/chemuson/app/` y declara sus
  dependencias directas reales.
- `architecture-boundaries`: El bootstrap permanece como capa superior; ningún
  módulo M00-M18 puede depender de M19.

## Impact

- Afecta M19 (`bootstrap`) y reduce responsabilidad de M08
  (`gui.main_window`) sin cambiar la jerarquía de widgets ni controllers.
- Las dependencias objetivo de M19 serán M08 (ventana), M15 (crash reporter) y
  M18 (metadata de versión).
- Conserva los entry points públicos `chemuson` y `python -m chemuson`, además
  de `chemuson.__main__:main`.
- No reorganiza menús, toolbars, docks, canvas, `items.py`, controllers,
  serialización o lógica química.
- No incorpora dependencias externas ni cambios visuales.
