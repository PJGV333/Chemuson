# Startup Recognition

## Consumidores y entry points

- `src/chemuson/__main__.py` es el único consumidor de `run_app`.
- `src/chemuson/gui/main_window.py` es la única definición de `run_app`.
- `pyproject.toml` registra `chemuson = "chemuson.__main__:main"`.
- No existe consumidor que requiera conservar un alias de compatibilidad en
  `gui.main_window`.

## Dependencias directas observadas

La función histórica de arranque usa directamente:

- M08: `QApplication`, `QMessageBox` y `ChemusonWindow`;
- M15: `crash_reporter.install()` y `write_crash_log()`;
- M18: `get_app_version()`.

Esto confirma las dependencias objetivo M19 → M08/M15/M18 descritas en el
diseño. El reconocimiento no encontró una dependencia inversa que obligue a
introducir una excepción o ciclo.

## Caracterización no Qt

`tests/test_application_startup_non_qt.py` ejecuta las rutas CLI en procesos
Python frescos y caracteriza la función histórica mediante su AST y
colaboradores falsos. De ese modo verifica aislamiento de imports, delegación,
orden del lifecycle, recuperación de autosaves, event loop, códigos de salida
y los dos caminos de reporte de fallos sin importar PyQt6 ni construir widgets.
