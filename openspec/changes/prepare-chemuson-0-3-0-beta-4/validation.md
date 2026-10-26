# Validación beta.4

## Identidad y baseline

- Baseline pre-cambio, árbol limpio y SHA inicial `9c171a8267d632e9bd914392680c316c28fe5e87`: [`baseline.md`](baseline.md).
- Commit de implementación validado: `4d91879d8aa7e154e5258de274d156302d60158c` (`feat(release): prepare beta.4 branding and update tracks`), rama `release/v0.3.0-beta.4-prep`.
- Entorno: `.venv`, Python 3.14.8, pytest 9.1.1, Ruff 0.16.10. Estado limpio en el SHA validado.
- No se ejecutó el suite completo `pytest -q`, conforme a la política del repositorio.

## Verificaciones locales en el SHA indicado

| Comando | Resultado |
|---|---|
| `.venv/bin/python -m compileall src tests tools packaging` | PASS |
| `.venv/bin/python -m pytest -q tests/test_branding.py tests/test_ui_icon_provider_hidpi.py tests/test_ui_svg_icons.py tests/test_packaged_icon_smoke.py tests/test_appimage_validation.py tests/test_linux_distribution_scripts.py tests/test_preview_build.py tests/test_preview_workflow_contract.py tests/test_release_version_script.py tests/test_version_metadata.py tests/test_release_artifacts.py tests/test_release_workflow_contract.py tests/test_update_core.py tests/test_update_policy.py tests/test_update_provider.py tests/test_update_semver.py tests/test_update_windows.py tests/test_iupac_ui.py` | PASS — 184 tests |
| `.venv/bin/python -m pytest -q tests/architecture` | PASS — 283 tests |
| `.venv/bin/ruff check src/chemuson/app/bootstrap.py src/chemuson/gui/theme/icon_smoke.py packaging/release/validate_appimage.py packaging/release/validate_packaged_icons.py packaging/release/validate_release_artifacts.py tests/test_branding.py tests/test_iupac_ui.py tests/test_linux_distribution_scripts.py tests/test_packaged_icon_smoke.py tests/test_appimage_validation.py tests/test_release_artifacts.py tests/test_update_provider.py --select F401,F811,F821,E722,E741` | PASS |
| `bash -n packaging/linux/build_appimage.sh packaging/linux/build_flatpak.sh packaging/linux/build_linux_all.sh` | PASS |
| `appstreamcli validate packaging/flatpak/io.github.PJGV333.Chemuson.metainfo.xml` | PASS (one non-fatal redundant-metadata notice) |
| `openspec validate prepare-chemuson-0-3-0-beta-4 --strict` | PASS |
| Local Markdown link scan, 14 changed/OpenSpec Markdown files | PASS — no missing local targets |
| `git diff --check` | PASS |

Los tests de branding ejercitan renderizado en temas claro/oscuro y DPR 2 desde el código fuente. No equivalen a smoke de un ejecutable congelado ni a aceptación manual de paquetes instalables.

## Logo, versión y mantenimiento

- PNG maestro: 1254×1254 RGB, SHA-256 `aab2600d2d46d87dba9dcea4266ae89ee4bbceafef0b05892b10713c310e29e6`; las variantes y el SVG incrustan el arte raster original, sin vectorizarlo.
- `_version.py` y AppStream: `0.3.0-beta.4`, fecha `2026-10-26`.
- README español/inglés y notas bilingües documentan la identidad como emblema de ChemUSON, no como escudo institucional, y agregan a Jesús David Albino-Espíndola sin asterisco.
- La lista pública de beta.3 conserva sus enlaces, porque beta.4 no se ha publicado. README, política y manual aclaran que `prerelease=false`/**Latest** fue una excepción manual de beta.3, no una política futura.
- Los contratos beta.3/preferencias del updater integrado se preservan. AppImageUpdate upstream fue inspeccionado estáticamente en `AppImageCommunity/AppImageUpdate` commit `a211784dfc746fdb6d8d32d6bb39add451c1ddeb`: `prerelease` es tratado como tag literal; las nuevas builds beta usan `latest-all` y stable conserva `latest`. La herramienta local `appimageupdatetool` no está instalada; no se afirma una actualización dinámica.

## Limpieza local

Se eliminaron selectivamente `dist/`, `dist-appimage/`, `dist-appimage-final/`, `dist-preview/`, `dist-preview-final/`, `build/` y el sidecar root obsoleto `Chemuson-v0.3.0-beta.1-linux-x86_64.AppImage.zsync`. La medición previa registró `2,499,985,408` bytes asignados para los seis directorios y `607,305` bytes para el sidecar (total `2,500,592,713` bytes asignados). No se usó `git clean`; caches y `.venv` permanecen.

Se conservaron los dos manifests trackeados de `dist-appimage-diagnostic/`, byte por byte idénticos a `main`:

- `checksums.sha256`: `09b1596d737272d0bb0509aed43502e2bd92b363b4f2058922e24fcbcc45b100`
- `diagnostic-provenance.json`: `dca186f9e76c8dc30e767fcfcb607fea7c690773447b9f52df421b229ccaa953`

No hay ejecutables o paquetes trackeados en las rutas limpiadas; los tags/releases beta.1–beta.3 se dejaron intactos.

## Preview y aceptación

- `gh auth status`: sin sesión autenticada. No se hizo push, PR ni dispatch de GitHub Actions.
- PyInstaller, `flatpak-builder`, `appimagetool` y `appimageupdatetool` no están disponibles localmente. **Los cuatro grupos oficiales de paquetes preview — Windows portable, Windows installer, AppImage Type 2 y Flatpak — están NOT BUILT.** Las comprobaciones estáticas/locales no se reportan como compilaciones.
- No se modificaron permisos, tags, Releases, canales ni `gh-pages`; no se creó `v0.3.0-beta.4` ni Release.
- Aceptación manual visual/funcional de los paquetes exactos: pendiente del propietario.

El siguiente paso autorizado es publicar la rama conforme a las protecciones del repositorio, ejecutar Build Preview oficial y revisar conjuntamente provenance, checksums y los cuatro artifacts. La creación del tag permanece bloqueada hasta aceptación explícita del propietario.
