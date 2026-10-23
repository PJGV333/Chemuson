# Validación — simplificar pruebas y organización de builds

Entorno local: CachyOS, Python 3.14.8, rdkit 2026.03.6, PyQt6 6.11.0 (pruebas locales, no atribuidas a CI). SHA comprobado: `7d935669300053ee222fc664444bc00a37e558c9` (stabilize) sobre `201a81cb…` para los contratos.

## Verificaciones locales

- `tests/test_linux_distribution_scripts.py`, `test_test_workflow_contract.py`, `test_preview_workflow_contract.py`, `test_release_workflow_contract.py`, `test_appimage_validation.py`, `test_ci_pytest_shards.py`: **49 passed** con `test.yml` solo-dispatch y scripts con `CHEMUSON_OUTPUT_DIR`.
- `bash -n` en los tres scripts modificados; `compileall src tools packaging` OK; Ruff focal: solo el `F401 math` basal preexistente en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py`; `git diff --check` limpio; `openspec validate simplify-beta1-build-test-pipeline --strict` y del change de reparación: válidos.
- Plan local de shards (`tools/ci_pytest_shards.py plan --shards 8`): colección **2610 node IDs == manifiesto actualizado** (los 3 IDs de regresión del worker congelado se registraron; el plan fail-closed había detectado correctamente 2607≠2610 en los runs automáticos antiguos).

## Comprobación de la política en GitHub

- Pushes a `stabilize/v0.3.0-beta.1` (`20b8a191…`, `201a81cb…`, `7d935669…`), `chore/…` y refs RC: **ningún run nuevo del workflow `test`** tras quitar los triggers; los últimos runs de `test` corresponden a pushes anteriores (`20b8a191…` y previos).
- Protección de ramas consultada antes de tocar triggers: **ninguna rama protegida**, ningún status check obligatorio referenciaba los jobs de `test.yml`.
- `build-preview.yml` y `release.yml` sin `needs:` hacia `test.yml`; empaquetan solo producción; los validadores ligeros post-build incluyen ahora el smoke esencial 2D→3D (`geometry_flow_smoke`) — confirmado en el log del run #39881978237: `Frozen AppImage normal 2D->3D conformer flow smoke passed`.

## Directorios de salida

- Prueba local: `CHEMUSON_OUTPUT_DIR=/tmp/chemuson-out-check bash packaging/linux/build_appimage.sh … diagnostic stabilize/v0.3.0-beta.1` generó `Chemuson-v0.3.0-beta.1-diagnostic-201a81cb-linux-x86_64.AppImage` fuera del repositorio; `git status --short` quedó limpio. Con `BUILD_TYPE=diagnostic` sobre una rama `release/**-prep`, el validador fail-closed rechazó correctamente el build (control de canal intacto).
- `.gitignore` ampliado (`dist-appimage-*/`, `dist-preview*/`, `dist-flatpak-preview/`, `pyinstaller-work/`); los manifiestos de procedencia/checksums trackeados del build diagnóstico permanecen trackeados y sin reescribir.

## Resultado de la política aplicada

- Build Preview #39881978237 (beta integrada `7d935669…`): success en los cuatro formatos sin ejecutar la suite completa; report `complete=true`, `publication=false`. Artefactos descargados y verificados localmente (SHA-256 vs checksums del run; AppImage Type 2 `AI\x02`): AppImage `5611d269…`, Windows portable `9871633f…`, instalador `1771cc9b…`, Flatpak `458d6e0f…`; provenance de los cuatro con `git_sha=7d935669…`.
