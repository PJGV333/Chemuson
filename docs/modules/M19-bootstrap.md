# M19 — Arranque y composición de la aplicación

## Responsabilidad

M19 es la capa superior que transforma la invocación del CLI en una aplicación
de escritorio. Posee el parser ligero, el composition root de Qt y el ciclo de
vida del proceso; no posee la ventana ni comportamiento de dominio.

## Rutas

- `src/chemuson/__main__.py`: entry point público, parsing de `--version` y
  dispatch diferido.
- `src/chemuson/app/__init__.py`: paquete liviano, sin imports gráficos.
- `src/chemuson/app/bootstrap.py`: composición interna de `QApplication`,
  crash reporter y `ChemusonWindow`.

## Dependencias

M19 depende directamente de:

- M08 para `ChemusonWindow`.
- M15 para `crash_reporter`.
- M18 para `get_app_version`.

La dirección es terminal: M00–M18 no pueden importar `chemuson.app` ni
`chemuson.__main__`.

## API y límites

La única API pública es `chemuson.__main__:main`, usada también por el script
instalado `chemuson`. `app.bootstrap.run_app()` es una operación interna.

Importar el entry point, solicitar ayuda o consultar `--version` no debe cargar
PyQt6, `chemuson.gui` ni el bootstrap gráfico. M19 no contiene menús, toolbars,
docks, controllers, canvas, serialización ni lógica química. La extracción
futura del application shell pertenece a otro cambio.

## Pruebas

`tests/test_application_startup_non_qt.py` verifica aislamiento de imports,
delegación única, ownership, secuencia de lifecycle, códigos de salida y
comportamiento de fallos sin necesitar un display Qt.
