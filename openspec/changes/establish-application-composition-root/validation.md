# Validación

Fecha: 2026-07-23

## Resultado

- `openspec validate establish-application-composition-root --strict`: válido.
- `openspec validate --all --strict`: 19/19 válidos.
- `pytest -q tests/test_application_startup_non_qt.py tests/architecture`:
  206 passed.
- `compileall` focalizado y global: correcto.
- `git diff --check`: correcto.
- Ruff global con `F401,F811,F821,E722,E741`: conserva únicamente el F401
  histórico de `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py`;
  el diff de este cambio no introduce diagnósticos.
- Suite completa: la colección se detiene con 56 errores porque la imagen no
  contiene `libEGL.so.1`. Todos son errores de importación de PyQt6 y reproducen
  la limitación ambiental registrada en `baseline.md`.

## Evidencia de preservación

Las pruebas no Qt ejecutan el cuerpo trasladado con colaboradores falsos y
confirman el orden histórico completo, autosaves, event loop, exit code y las
dos rutas de error. También comprueban que ayuda, versión e importación del
entry point no cargan PyQt6, GUI ni `chemuson.app.bootstrap`.

El cambio no modifica menús, canvas, items, controllers, serialización,
renderizado ni lógica química. La aplicación gráfica real no puede abrirse en
esta imagen hasta que esté disponible la biblioteca de sistema `libEGL.so.1`.
