## Why

El composition root ya separa el ciclo de vida de la aplicación, pero
`gui/main_window.py` conserva 2603 líneas y `ChemusonWindow.__init__` todavía
ensambla en un solo bloque controllers, pestañas, docks, toolbars, barra de
estado y conexiones globales. Esta concentración dificulta revisar el orden de
construcción, probar el shell sin recorrer lógica química y continuar
reduciendo la ventana de manera segura.

La aplicación inicia correctamente en un entorno Qt real después del cambio
anterior. El siguiente incremento debe extraer únicamente la composición
visual de alto nivel, sin modificar apariencia, acciones ni comportamiento.

## What Changes

- Crear `src/chemuson/gui/shell/` dentro de M08.
- Extraer la construcción y el wiring de alto nivel de
  `ChemusonWindow.__init__` a un ensamblador explícito del application shell.
- Mantener `ChemusonWindow` como propietario de handlers, lifecycle de la
  ventana y API pública.
- Reutilizar `MainWindowUiBuilder` y las factorías de acciones existentes; no
  duplicarlas ni sustituirlas.
- Caracterizar por AST la delegación única al shell, su ownership y la ausencia
  de lógica química/canvas movida.
- Actualizar catálogo y documentación de M08.

## Capabilities

### New Capabilities

- `application-shell`: Ensamblaje explícito de la ventana principal, sus
  regiones visuales y conexiones globales dentro de M08.

### Modified Capabilities

- `module-catalog`: M08 registra el paquete `gui/shell/` como implementación
  interna y sus pruebas de arquitectura.

## Impact

- Reduce responsabilidad y tamaño de `gui/main_window.py` sin cambiar la ruta
  pública `chemuson.gui.main_window.ChemusonWindow`.
- Conserva exactamente menús, toolbars, docks, pestañas, status bar, señales,
  orden de construcción y chequeo diferido de actualizaciones.
- No toca canvas, items, controllers, química, persistencia, formatos,
  renderizado ni estilos visuales.
- No añade dependencias externas, excepciones ni ciclos.
