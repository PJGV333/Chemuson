# Baseline antes de cambios de campaña

## Identidad y protección

- Fecha UTC: `2026-10-26T17:19:23Z`.
- Repo: `PJGV333/Chemuson`; remote `origin=https://github.com/PJGV333/Chemuson.git`.
- Estado inicial: `main`, `git status --short` vacío, SHA `9c171a8267d632e9bd914392680c316c28fe5e87`.
- `git fetch --prune origin main`: correcto; `HEAD == origin/main`, `git rev-list --left-right --count HEAD...origin/main` = `0 0`.
- Rama creada desde ese commit: `release/v0.3.0-beta.4-prep`.
- Tags existentes/anotados preservados: `v0.3.0-beta.1` → `85fa4ed2bd0afe1d31953440228fb73a69b54e7e`; `v0.3.0-beta.2` → `6aecba5a5585d9560b693186578b76bea9e623b8`; `v0.3.0-beta.3` → `01652e1c9e79945a050f3733adc307b47866b7bb`. Remote y local beta tags coinciden. Ningún tag se modificó.
- `gh auth status`: sin sesión autenticada. La API pública de GitHub sí respondió; esto no habilita push, PR ni Actions dispatch.

## PNG oficial

- `/home/unison-pjgv/Descargas/logo_ChemUSON.png`: existe; `file` lo identifica como PNG RGB, 1254×1254, no interlazado.
- Pillow `Image.verify()`: PASS; un frame; SHA-256 `aab2600d2d46d87dba9dcea4266ae89ee4bbceafef0b05892b10713c310e29e6`.
- Inspección de imagen: emblema circular aprobado con búho, matraz, montañas/paisaje sonorense, libro y texto ChemUSON; fondo blanco. Sin reinterpretación.

## Verificaciones locales previas

Entorno: `.venv/bin/python` Python 3.14.8, pytest 9.1.1, Ruff 0.16.10.

| Comando | Resultado |
|---|---|
| `.venv/bin/python -m compileall src tests tools packaging` | Exit 0 |
| `.venv/bin/python -m pytest -q tests/test_branding.py tests/test_ui_icon_provider_hidpi.py tests/test_ui_svg_icons.py tests/test_packaged_icon_smoke.py tests/test_appimage_validation.py tests/test_linux_distribution_scripts.py tests/test_preview_build.py tests/test_preview_workflow_contract.py tests/test_release_version_script.py tests/test_version_metadata.py tests/test_release_artifacts.py tests/test_release_workflow_contract.py tests/test_update_core.py tests/test_update_policy.py tests/test_update_provider.py tests/test_update_semver.py tests/test_update_windows.py` | 179 passed in 6.19s |
| `.venv/bin/python -m pytest -q tests/architecture` | 283 passed in 16.72s |
| `.venv/bin/ruff check src/chemuson/gui/main_window.py src/chemuson/gui/shell/assembly.py src/chemuson/app/bootstrap.py tests/test_branding.py packaging/release/set_version.py --select F401,F811,F821,E722,E741` | All checks passed |
| `git diff --check` | Baseline sin cambios |
| `pytest -q` | No ejecutado, por instrucción del propietario |

## Diagnóstico de artefactos antes de limpiar

Comandos ejecutados: `git status --short --ignored -- <rutas>`, `git ls-files <rutas>`, `du -sB1 <rutas>`, `find ... -type f`, `file ...`, `git clean -ndX`.

- `dist/`: executable ELF PyInstaller antiguo, 356,786,176 bytes asignados.
- `dist-appimage/`: AppImage beta.1, updateinfo, update.json y zsync antiguos, 355,786,752 bytes.
- `dist-appimage-final/`: AppImage beta.1 final y sidecars antiguos, 355,790,848 bytes.
- `dist-preview/`: dos AppImages preview beta.1, 710,340,608 bytes.
- `dist-preview-final/`: un AppImage preview beta.1, 355,172,352 bytes.
- `build/`: TOCs, PyInstaller PKG/PYZ, caché y temporales, 366,108,672 bytes.
- Suma exacta de las seis salidas generadas (sin diagnostic): `2,499,985,408` bytes asignados (`2.328 GiB`, `2.500 GB`).
- `dist-appimage-diagnostic/`: 8,192 bytes, contiene `checksums.sha256` y `diagnostic-provenance.json`, ambos registrados por Git; preservar.
- Artefacto externo a esas carpetas: `Chemuson-v0.3.0-beta.1-linux-x86_64.AppImage.zsync`, ignorado, 607,305 bytes; sidecar generado antiguo sin referencias documentales.
- No hay ejecutables/paquetes trackeados en esas rutas; sólo los dos manifiestos anteriores. `git clean -ndX` listó 57 candidatos, incluyendo `.venv`, caches, pycache y el zsync; no se usará limpieza global.

## Estado beta.3 verificado remotamente (API pública, 2026-10-26)

- API `/releases/tags/v0.3.0-beta.3`: tag beta.3, `prerelease=false`, `draft=false`.
- API `/releases/latest`: también devuelve beta.3; actualmente tiene el distintivo GitHub **Latest** por la excepción manual del propietario.
- El asset publicado `.AppImage.updateinfo` conserva `gh-releases-zsync|PJGV333|Chemuson|prerelease|Chemuson-v0.3.0-beta.3-linux-x86_64.AppImage.zsync`.
- La herramienta `appimageupdatetool`/`appimageupdate` no está instalada localmente; no se afirma ejecución dinámica. Revisión estática del upstream `AppImageCommunity/AppImageUpdate` en commit `a211784dfc746fdb6d8d32d6bb39add451c1ddeb`, archivo `src/updateinformation/GithubReleasesZsyncUpdateInformation.cpp`: reconoce `latest`, `latest-pre` y `latest-all`; cualquier otro valor se concatena como tag en `/releases/tags/<valor>`. Por lo tanto, el `prerelease` literal actual no es un selector válido. `latest-all` recorre releases de cualquier clasificación; `latest-pre` filtra `prerelease=true`. Para beta, `latest-all` es coherente con la política existente de permitir graduación a final y con beta.3 etiquetada manualmente `prerelease=false`.
