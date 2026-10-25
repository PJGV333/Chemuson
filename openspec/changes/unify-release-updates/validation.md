# Validación — unificación progresiva de releases y actualizaciones

## Rama y preparación de integración

- Rama de trabajo: `refactor/unified-release-updates`, creada desde el tag/HEAD beta.3 `01652e1c9e79945a050f3733adc307b47866b7bb`.
- `origin/main` apuntaba a `0f8fd6d7847756bad81da6ee62cb81332b1f4cd9`; merge-base original `f0fde603371255902bf0e630ca1ce032b8f16ad3`. Se inspeccionó antes de integrar: la divergencia de `main` era un commit documental (`README.en.md`, `README.md`, `docs/MANUAL_USUARIO.md`, `docs/manual_usuario.md`).
- La integración documental fue un merge explícito **hacia la rama de trabajo**, no hacia `main`. Commit `5914d0332a9718c5ee13a3c5fa351f46f0cea311`, padres: OpenSpec inicial `f826688ef4070f54d928eae6291a6b8928601115` y `origin/main` `0f8fd6d7847756bad81da6ee62cb81332b1f4cd9`. README/manuales de beta.3 conservados; sólo se normalizaron saltos de línea Markdown que producían trailing whitespace.
- Tag `v0.3.0-beta.3` permanece apuntando a `01652e1c9e79945a050f3733adc307b47866b7bb`. No se creó ni movió tag, no se editó Release/paquetes y no se integró el trabajo en `main`.

## Resultado técnico

- El proveedor acepta como identidad sólo un tag oficial `v<SemVer>` válido; no rescata una etiqueta inválida desde el título. Ordena por precedencia SemVer y filtra canales desde el SemVer, ignorando el booleano GitHub `prerelease` como autoridad. Fechas/flags contradictorios, `beta.10`/`beta.3` y títulos válidos con tags inválidos tienen regresiones.
- SemVer rechaza identificadores prerelease/build vacíos, ceros iniciales en identificadores prerelease numéricos y dígitos Unicode en la versión base; metadatos build numéricos con ceros iniciales siguen siendo válidos.
- Las preferencias beta/estable, downgrade guard, SHA-256, selección de asset Windows, tracks AppImageUpdate y remotos Flatpak se conservaron. El test simula `beta.2 → beta.3 → 0.3.0`; verifica no downgrade, asset ausente, fallo de red y checksum discrepante.
- Un único `release.yml` sigue publicando los mismos paquetes/gates/permisos; las notas ahora usan `docs/release/<version>-public.md`. Preview y diagnóstico siguen siendo no publicadores. La validación de diagnóstico acepta `main`, `feature/*` (mantiene `feat/*`) y `fix/*`, y sigue rechazando `release/*` para usar el Preview correspondiente.
- No hubo cambios a la clasificación de beta.3, builder AppImage, manifest/remotos Flatpak ni selector installer/portable de Windows.

## Verificaciones focales

- `env -u APPIMAGE python -m pytest -q tests/test_update_semver.py tests/test_update_policy.py tests/test_update_provider.py tests/test_update_core.py tests/test_update_security.py tests/test_update_portable.py tests/test_update_windows.py tests/test_update_telemetry.py tests/test_preview_build.py`: **58 passed in 0.95s**.
- `env -u APPIMAGE python -m pytest -q tests/test_release_tag_validation.py tests/test_release_artifacts.py tests/test_release_workflow_contract.py tests/test_preview_workflow_contract.py tests/test_linux_distribution_scripts.py tests/test_appimage_validation.py -k 'not graphics_smoke_disables_background_geometry_generation'`: **70 passed, 1 deselected in 1.68s**. El caso excluido importa RDKit, no instalado en `/usr/bin/python` de este entorno; no se modificó ni se ocultó el test.
- `python -m compileall -q src tests tools packaging`: exit 0.
- Ruff focal (`F401,F811,F821,E722,E741`) sobre updater, tests de update/release y `packaging/release`, `packaging/linux`: **All checks passed!**.
- `openspec validate unify-release-updates --strict --no-interactive`: **Change is valid**.
- `git diff --check`: exit 0, sin salida. La revisión histórica `git diff --check origin/main...HEAD` había señalado whitespace en snapshots OpenSpec preexistentes de la línea beta; no se alteraron esos archivos ajenos.
- No se ejecutó `pytest --collect-only` ni la suite completa de aproximadamente 2.600 tests, según alcance explícito. Baseline detallado y fallos ambientales están en `baseline.md`.

## Publicación de rama y cierre

- Se publicó por push normal `refactor/unified-release-updates`; GitHub devolvió `https://github.com/PJGV333/Chemuson/pull/new/refactor/unified-release-updates`. No se creó PR.
- Commit de implementación verificado: `ad9c3ed196bb63df78c3af9de7fdd513eb210b97`. Después del push, `HEAD` y `origin/refactor/unified-release-updates` coincidían en ese SHA y `git status -sb` mostró únicamente el tracking limpio.
- El commit documental de `main` quedó como segundo padre del merge de preparación; `origin/main` continúa en `0f8fd6d7847756bad81da6ee62cb81332b1f4cd9` y no recibió cambios de esta campaña.
- `v0.3.0-beta.3` continúa apuntando en GitHub al commit `01652e1c9e79945a050f3733adc307b47866b7bb` (tag objeto `69da3ab3e0c122da1dc8dc6e706d95eb7f9f8f9b`). No se crearon tags ni se ejecutó workflow de publicación.
- Se detuvo antes de merge a `main`, sin activar `release.yml`. Cualquier commit documental de cierre posterior es sólo el registro de esta verificación; el SHA exacto final de la rama se informa en la entrega.
