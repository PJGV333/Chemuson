# Baseline — unificación progresiva de releases y actualizaciones

## Identidad y estado

- Fuente: rama `release/v0.3.0-beta.3-prep`, tag/HEAD `01652e1c9e79945a050f3733adc307b47866b7bb` (beta.3).
- Estado inicial de esa rama: `## release/v0.3.0-beta.3-prep...origin/release/v0.3.0-beta.3-prep`, sin cambios locales.
- Rama de trabajo creada desde ese SHA: `refactor/unified-release-updates`; estado inicial: `## refactor/unified-release-updates`, sin cambios locales.
- `origin/main` apunta a `0f8fd6d7847756bad81da6ee62cb81332b1f4cd9`; merge-base con beta.3: `f0fde603371255902bf0e630ca1ce032b8f16ad3`. Las dos pruebas de ancestralidad (`main` ancestro de beta.3 / beta.3 ancestro de `main`) devolvieron código 1. `origin/main` tiene un único commit divergente `docs: portada bilingüe y manual actualizado para beta.3`, limitado a `README.en.md`, `README.md`, `docs/MANUAL_USUARIO.md` y `docs/manual_usuario.md`; se verificó el diff antes de incorporarlo a esta rama.

## Comprobaciones previas a cambios

- `python -m compileall -q src tests tools packaging`: exit 0, sin salida.
- Ruff focal sobre `src/chemuson/update`, tests de update/release y `packaging/release`, `packaging/linux`, selector `F401,F811,F821,E722,E741`: exit 0, `All checks passed!`.
- Primer intento focal: 116 tests; 114 pasaron y 2 fallaron por entorno, no por cambios de esta campaña:
  - `tests/test_update_portable.py::test_detect_portable_update_context_for_windows_portable`: el entorno exporta `APPIMAGE=/home/ccachyavgp/Aplicaciones/Strata-0.2.3-linux-amd-rx9070.AppImage`, que contamina la prueba Windows.
  - `tests/test_appimage_validation.py::test_graphics_smoke_disables_background_geometry_generation`: el Python activo es `/usr/bin/python` y no tiene `rdkit` (`ModuleNotFoundError`).
- Reejecución focal con `env -u APPIMAGE`, omitiendo sólo `tests/test_appimage_validation.py` que necesita RDKit: `103 passed in 1.88s` (exit 0).
- No se ejecutó `pytest --collect-only` ni la suite completa; la instrucción de esta campaña limita la verificación a tests focales y prohíbe pytest masivo. El caso de AppImage que importa RDKit queda como limitación ambiental y no se alterará empaquetado para ocultarlo.
- Exploración pre-cambios de `git diff --check origin/main...HEAD` detectó whitespace final únicamente en snapshots OpenSpec históricos traídos por la línea beta; no son cambios de esta campaña. La verificación final se hará sobre los commits/diff propios, sin refactorizar evidencia ajena.

## Salidas resumidas exactas

```text
python -m compileall -q src tests tools packaging
(exit 0; stdout vacío)

ruff check ... --select F401,F811,F821,E722,E741
All checks passed!

env -u APPIMAGE python -m pytest -q [14 archivos update/release/Linux, sin test_appimage_validation.py]
103 passed in 1.88s
```
