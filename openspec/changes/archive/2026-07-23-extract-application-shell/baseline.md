# Baseline Report

- **Fecha**: 2026-07-23
- **Rama base**: `main`
- **HEAD base**: `dfaca3e`
- **Rama de trabajo**: `architecture/phase7-extract-application-shell`
- **Entorno**: Python 3.12; imagen Work sin `libEGL.so.1`

## Comandos y resultados

### `git status --short`

Árbol limpio antes de archivar el cambio anterior y crear esta rama.

### `python -m compileall -q src tests tools packaging`

Exit code 0.

### `pytest --collect-only -q`

Con `/workspace/.venvs/chemuson/bin/python`: exit code 2. Pytest detectó 1033
pruebas y se detuvo con 56 errores de importación por ausencia de
`libEGL.so.1`, la misma limitación ambiental ya documentada.

### Ruff focalizado

Con `/workspace/.venvs/chemuson/bin/ruff`: conserva únicamente el F401
histórico de
`tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.

### OpenSpec

El cambio anterior fue archivado después de validación individual. La
validación global posterior pasó 19/19 y no quedaron cambios activos antes de
crear `extract-application-shell`.

## Reconocimiento

- `src/chemuson/gui/main_window.py`: 2603 líneas.
- `ChemusonWindow.__init__`: líneas 146–361 antes de la extracción.
- El constructor mezcla estado/controllers, tabs, docks, toolbars, status bar
  y conexiones.
- `MainWindowUiBuilder` y `gui/actions/` ya poseen construcción detallada y
  deben reutilizarse.
- La prueba manual del usuario confirma que `main` inicia normalmente antes de
  esta fase.
