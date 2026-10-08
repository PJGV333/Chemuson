# Validación — preparación v0.3.0-beta.1 y builds preview

Fecha: 2026-10-08. Rama local: `release/v0.3.0-beta.1-prep`. Base: `8640bed18f0961ef9582f8376e98e9a33dbc3c01`.

## Código y contratos

- `openspec validate prepare-chemuson-0-3-0-beta-1 --strict`: PASS.
- `timeout 5m python -m pytest -q tests/architecture`: **280 passed**.
- `timeout 8m python -m pytest -q` con suites de versión, updater, paquete, release, preview, workflows y scripts Linux: **104 passed in 1.76s**.
- El test preview ejecuta un stub portable para comprobar nombres y ausencia de sidecars; **no** equivale a PyInstaller/Flatpak/Inno real.
- `python -m compileall -q src tests tools packaging`: PASS.
- Ruff scoped sobre scripts/tests nuevos y editados (`F401,F811,F821,E722,E741`): PASS.
- Ruff global (`src tests tools packaging`, misma selección): exit 1 sólo por el F401 histórico `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no se modificó.
- `timeout 5m python -m pytest --collect-only -q`: **2072 tests collected**.
- `git diff --check`: PASS. `bash -n` para los tres builders Linux: PASS.
- AppStream XML, Flatpak manifest YAML y los dos workflows parsean con los parsers locales. No hay `actionlint`; por tanto esto es verificación sintáctica, no validación de esquema GitHub Actions.

## Pruebas no ejecutadas / límite de baseline

- No se ejecutó la suite monolítica completa, ni se repitió el crash Qt/SIGSEGV histórico. `baseline.md` conserva la evidencia anterior; no se atribuye causa ni se declara resuelto.
- No se hicieron compilaciones reales: `pyinstaller`, `flatpak-builder`, `appimagetool` e Inno Setup (`iscc`) no están instalados localmente. El `.AppImage`-named sigue siendo un ejecutable portable PyInstaller, no Type 2.
- `gh auth status`: no hay sesión autenticada. El workflow preview no se ha subido a GitHub ni se inició mediante CLI/API. No se improvisó una vía con permisos mayores.
- Manual acceptance: 70 casos preparados, todos `NOT TESTED`; no se ejecutaron contra instaladores.

## Aislamiento, estado Git y alcance

- El workflow preview usa `contents: read`, rama `release/**-prep`, versión canónica y checkout de un SHA completo compartido. Sus contratos estáticos cubren artefactos, existencia, checksums/procedencia, no-publicación y separación del workflow oficial.
- El workflow preview no está disponible en GitHub desde este árbol local. No existen artifacts reales de esta campaña.
- No se creó tag ni GitHub Release; no hubo comando de publicación, cambio de `gh-pages`, remoto Flatpak ni manifest público beta/stable. El contenido actual no llegó al updater.
- Diff de `src/chemuson/clean2d/`, dependencias (`pyproject.toml`, `requirements*.txt`) y persistencia `.cmsn`: vacío. No se cambiaron algoritmos químicos ni Molecular Assistant.
- No hubo commit ni push: el árbol queda para revisión del propietario. `gh auth status` confirmó ausencia de credenciales. No se hizo merge, rebase ni force-push.

## Dictamen

- `PREVIEW BUILD INFRASTRUCTURE: NOT READY` para ejecución remota hasta que la rama/workflow esté disponible en GitHub y se revise el aislamiento allí. Contratos estáticos locales: PASS.
- `PREVIEW ARTIFACTS: NOT BUILT` (Windows portable, Windows setup, Linux portable y Flatpak).
- `MERGE READINESS: NOT READY` — pendiente de revisión/decisión del propietario, publicación normal de la rama si se autoriza, ejecución real de previews y aceptación manual; no se solicita merge aquí.
- `RELEASE READINESS: NOT READY TO PUBLISH` — no se creó `v0.3.0-beta.1`; los builds reales y la matriz manual siguen pendientes, además de la deuda Qt documentada.
