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

## Corrección del gate de inicio no Qt — revisión adicional

- Baseline de esta revisión: rama `release/v0.3.0-beta.4-prep`, SHA `32ea9f05045a771874ec5e525e301bf855e8e9f6`, árbol limpio y sincronizado con `origin/release/v0.3.0-beta.4-prep`.
- El propietario reportó Build Preview #40371397035 sobre ese SHA: `Release gate bounded checks rehearsal`, 220 passed y 1 failed en `test_startup_success_order`. La API pública del run respondió HTTP 504 al intentar consultar sus detalles; el estado del run citado se atribuye al reporte del propietario.
- Reproducción local pre-cambio: `pytest -q tests/test_application_startup_non_qt.py::test_startup_success_order` falló con el `NameError: _application_icon_path` reportado y el `AttributeError` secundario por `write_crash_log` ausente. La colección del archivo tenía 11 node IDs; los contratos de branding/iconos previos pasaron (10 tests).
- Corrección sólo en tests: commit candidato `2a629aac1355d96f8885c6ccdc05c7e1136d00d5`. No se modificó código de producción. El harness ahora inyecta el path existente, `QIcon` no nulo, `QApplication.setWindowIcon`, `window.setWindowIcon`, `QMessageBox` y `crash_reporter.write_crash_log`; la aserción sigue comprobando el orden completo y que el mismo icon object llegue a ambos.
- Se añadió un caso de archivo ausente: verifica `FileNotFoundError`, reporte en crash log y mensaje crítico existente, sin construir icono ni ventana. Los tests previos de fallo de arranque genérico permanecen intactos.
- En `2a629aac...`: `.venv/bin/python -m compileall src tests tools packaging` PASS; `.venv/bin/python -m pytest -q tests/test_application_startup_non_qt.py tests/test_branding.py tests/test_packaged_icon_smoke.py` PASS — 22 tests; Ruff focal en esos tres archivos PASS; OpenSpec strict y `git diff --check` PASS.
- La colección completa (sin ejecutar tests) encontró 2,653 node IDs. El manifiesto previo contenía 2,630. Se agregó únicamente el nuevo node ID: el manifiesto final contiene 2,631; la comparación detecta 22 node IDs adicionales de archivos preexistentes ajenos a este cambio. No se añadieron para ocultar o rebaselinear deuda fuera de alcance, y no se ejecutó el suite masivo. Esa discrepancia preexistente queda documentada; el preview gate ejecuta su lista bounded explícita, no el plan de shards de `test.yml`.

## Preview y aceptación

- `gh auth status`: sin sesión autenticada. `git push --dry-run origin release/v0.3.0-beta.4-prep` falló antes de cualquier escritura con `could not read Username for 'https://github.com'`. La rama remota continúa en el SHA base `32ea9f05045a771874ec5e525e301bf855e8e9f6`; el commit corregido `2a629aac1355d96f8885c6ccdc05c7e1136d00d5` permanece local. No se abrió PR ni se hizo merge.
- El intento de consultar por API el run #40371397035 recibió HTTP 504; el estado del nuevo código no está verificado en GitHub Actions. No hay un Build Preview nuevo ni URL de run que reportar.
- PyInstaller, `flatpak-builder`, `appimagetool` y `appimageupdatetool` no están disponibles localmente. **Los cuatro grupos de paquetes corregidos — Windows portable, Windows installer, AppImage Type 2 y Flatpak — están NOT BUILT.** Las comprobaciones estáticas/locales no se reportan como compilaciones.
- No se modificaron permisos, tags históricos, Releases, canales ni `gh-pages`; no se creó `v0.3.0-beta.4` ni Release. Beta.4 permanece identificada como beta y beta.3 conserva su clasificación excepcional.
- Aceptación manual del propietario sobre los paquetes exactos del SHA corregido: pendiente.

Próximo paso: habilitar una sesión GitHub autenticada para push normal, abrir/revisar PR y ejecutar Build Preview oficial en un SHA único. No integrar si un gate falla o falta un artifact/provenance/checksum. La creación del tag permanece bloqueada hasta revisar los cuatro paquetes y confirmar su aceptación explícita.

## Reparación del timeout Flatpak — preparación del entorno (2026-10-26)

### Baseline antes del cambio

- Rama `release/v0.3.0-beta.4-prep`, HEAD `670a7670235ec1f8833cb03c7d2475bbf5a2589f`, `git status --short` vacío y sincronizado con `origin/release/v0.3.0-beta.4-prep` (mismo SHA en la rama remota).
- Evidencia pública de Actions (API, sin sesión): run #40375931525 `cancelled`; job `Linux Flatpak bundle` pasos `Install Flatpak build tools` 18:26:17Z→19:10:16Z (43 min 59 s), `Verify exact source identity` correcto, `Build local Flatpak bundle only` 19:10:16Z→19:12:40Z (2 min 24 s) `cancelled`, upload `skipped`; `Preview build report` `failure` (fail-closed).
- Comparación del mismo job en runs previos: #40371397035 preparación 37 min 24 s y job completo 41 min 28 s; #40163508452 0,52/3,70 min; #40162819727 0,43/4,18; #40156368359 0,67/4,00; #40144590463 0,48/4,52; #39900808044 0,32/4,68; #39881978237 0,35/4,80; #39679607716 0,67/5,02. En el run cancelado, el paso apt del job AppImage tardó 22 s y el del gate 28 s: no hubo caída general de apt.
- Contexto externo: los runners fijan `mirror+file:/etc/apt/apt-mirrors.txt` con `azure.archive.ubuntu.com` prioridad 1; los cierres grandes se estrangulan en ese nodo sin error (actions/runner-images #14823, #14594, #12949). El propio repo ya registró 123 paquetes/99 MB y ~38 min para el job Flatpak de `test.yml`. Los logs del run no son descargables sin permisos de administración (`Must have admin rights`), por lo que el detalle de la URL lenta no se afirma.

### Cambio implementado

Commits locales `5274c2e6e332eaab7ce00463c9dfe86a19fdebad` (preparación, workflows, builder y contratos) y `f9c72514…` (pasar las opciones `-o` de apt como argumentos separados, tras verificar en el arnés que el formalismo `-oKey=Value` quedaba como un solo token). HEAD local de la rama tras el registro: `f9c72514…`.

Archivos: `.github/workflows/build-preview.yml`, `.github/workflows/release.yml`, `packaging/linux/build_flatpak.sh`, `packaging/flatpak/toolchain.lock` (nuevo), `packaging/linux/install_flatpak_tools.sh` (nuevo), `tests/test_flatpak_toolchain_preparation.py` (nuevo), `tests/ci/expected_pytest_nodeids.txt`, `docs/release/PREVIEW_BUILDS.md` y los cuatro artefactos OpenSpec del cambio.

No se tocaron: código de la aplicación, Clean2D, visor 3D, RDKit, ChemName, asistente, persistencia, logo, política del actualizador, manifiesto Flatpak (runtime `org.kde.Platform` 6.10, cinco ruedas con `sha256`, smoke empaquetado, permisos), ni los controles de release (verificación de SHA, checksums, provenance, metadata OSTree/GPG, Pages, gates del reporte). Ningún `timeout-minutes` se aumentó.

### Validaciones locales ejecutadas en `5274c2e6`

- `.venv/bin/python` 3.14.8, pytest 9.1.1, OpenSpec CLI instalado localmente.
- YAML parseado con PyYAML (`BaseLoader`): build-preview 7 jobs, release 6 jobs, test 5 jobs. No hay `actionlint` en el equipo.
- `bash -n packaging/linux/build_flatpak.sh` y `bash -n packaging/linux/install_flatpak_tools.sh`: PASS.
- `pytest -q` focal (contracts preview/release/test-workflow, linux distribution scripts, appimage validation, preview build, release artifacts, tag validation, branding, packaged icon smoke, flatpak toolchain): **100 passed**.
- `pytest -q tests/architecture`: **283 passed**.
- `python -m compileall src tests tools packaging`: PASS. Ruff focal (`F401,F811,F821,E722,E741`) sobre el test nuevo: PASS.
- `openspec validate prepare-chemuson-0-3-0-beta-4 --strict`: valid.
- `git diff --check 670a7670..5274c2e6`: sin avisos; `git status --short` limpio tras el commit.
- Colección completa sin ejecutar tests: 2.664 node IDs; el manifiesto queda en 2.642 (2.631 + los 11 nuevos). La diferencia de 22 node IDs preexistentes se conserva documentada y no se rebaselinea.

### Arnés del preparador y del empaquetador

- Arnés local con `apt-get`, `curl` y `dpkg-query` simulados (fuera del repositorio, en `/tmp/fp-harness`) sobre el script real: ocho escenarios PASS — espejos sanos con instalación offline, nodo azure estrangulado excluido de la selección, rebote de `apt-get update` al segundo origen, rebote de la descarga acotada, caché cálida reutilizada, descarga fallida con repliegue acotado, paquetes de publicación (`ostree gnupg`) y fallo cerrado con `exit 3` cuando ningún origen oficial responde. En el modo `URIs:` directas la reescritura cambia sólo el stanza principal y conserva `Signed-By` y `noble-security`.
- Sonda real contra los cuatro orígenes: `azure.archive.ubuntu.com` responde sólo por HTTP (200), los otros tres por HTTPS; la sonda del script cubre ambos esquemas.
- Mock-run de `build_flatpak.sh` en modo preview con `flatpak`/`flatpak-builder` simulados: bundle `Chemuson-v0.3.0-beta.4-preview-670a7670-linux-x86_64.flatpak`, resumen de fases escrito y salida 0. Con `CHEMUSON_FLATPAK_BUILDER_TIMEOUT_SECONDS=1` devuelve 124 con el diagnóstico específico; con `CHEMUSON_FLATPAK_REPO_URL` definido sigue fallando cerrado (`exit 2`) porque un preview no puede usar remoto público.

### Pendiente de evidencia en GitHub

- `git push origin release/v0.3.0-beta.4-prep` vuelve a fallar sin credenciales (`could not read Username for 'https://github.com'`, `GIT_TERMINAL_PROMPT=0`); `gh auth status` sigue sin sesión. El commit `5274c2e6` permanece local y la rama remota continúa en `670a7670`.
- Por tanto el nuevo Build Preview todavía no existe: no se declara ninguna mejora de tiempos medida en Actions hasta comparar los pasos `Install Flatpak build tools` y `Build local Flatpak bundle only` del run nuevo contra los 43 min 59 s y 2 min 24 s del run cancelado.
- Los cuatro paquetes del candidato corregido siguen **NOT BUILT**: no hay artifact, checksum ni provenance que revisar, y la aceptación manual del logotipo en paquetes congelados continúa pendiente.

### Iteración 1 en GitHub: workflow rechazado por el parser (2026-10-26 19:58Z)

- La rama se publicó en `f7abf93cede468c06b2738c71a7797fac407936a`. El push creó dos runs, `Build Preview` #40 (`40391503539`) y `release.yml` #14 (`40391505542`), ambos `completed/failure` con **0 jobs**, `name` igual a la ruta del archivo y `created_at == updated_at == run_started_at`. Es la firma de un archivo rechazado antes de crear los jobs: no se ejecutó ningún paso, por lo que no hay tiempos de preparación que comparar.
- El detalle sólo es visible en la página pública del run: `Line: 469, Col: 39): Unrecognized named-value: 'runner'. Located at position 1 within expression: runner.temp`. Causa: `CHEMUSON_FLATPAK_APT_CACHE_DIR: ${{ runner.temp }}` estaba en el `env` del job, y el contexto `runner` no está disponible a nivel de job (sí en `env`/`with` de un paso). Los tres jobs afectados eran `build_linux_flatpak`, `build_flatpak` y `publish_flatpak_remote`.
- Corrección en `f6955f12…`: la lista de paquetes y el directorio de caché pasan al `env` del paso de instalación Flatpak; el paso de caché conserva `${{ runner.temp }}` en su `with.path`, que sí es un contexto válido ahí. Se añadieron dos guardias de contrato: ningún `env` de job puede referenciar contextos exclusivos de paso (`runner.`, `steps.`, `job.`, `env.`, `inputs.`) y ningún `uses:` de paso puede contener una expresión. El manifiesto de node IDs queda en 2.643 (12 pruebas del contrato Flatpak).
- La API sin sesión no expone logs ni anotaciones (`check-runs` vacío), así que el mensaje se obtuvo de la página pública del run. Los runs rechazados no generan artefactos ni consumen empaquetado.

### Refuerzo del control del origen apt (escenario I del arnés)

- Riesgo encontrado al revisar el script: si `/etc/apt/apt-mirrors.txt` del runner no fuera escribible, el script declaraba el origen aplicado sin comprobar que la reescritura hubiera tenido efecto, y con `set -e` el `tee` fallido abortaba el paso con salida 1. En ese caso el job seguiría dependiente del nodo asignado, que es justamente la causa del bloqueo.
- Corrección: `apply_origin` verifica la primera línea del `mirror+file` tras escribirlo; si no coincide, cambia a `sources-direct` y fija la línea `URIs:` del stanza principal al origen oficial sondeado. `rewrite_sources_uri` falla cerrado con `exit 3` si la lista de fuentes no puede leerse o escribirse, y la regla directa sólo sustituye una línea `URIs: mirror+file:...`, conservando `Signed-By` y el stanza `noble-security`.
- Arnés ampliado a nueve escenarios: el nuevo escenario I (`chmod 0444` del mirror list) termina en `exit 0`, instala offline y deja `URIs: https://archive.ubuntu.com/ubuntu/` en el stanza principal con el stanza de seguridad y las firmas intactas. Los otros ocho escenarios conservan su resultado esperado.
- Contrato estático nuevo `test_installer_verifies_the_origin_files_it_changes` (13 node IDs del archivo, manifiesto en 2.644) y verificador local de contextos por campo sobre los tres workflows: 0 avisos en la versión corregida, 3 avisos exactos en la versión rechazada (`jobs.build_linux_flatpak.env`, `jobs.build_flatpak.env`, `jobs.publish_flatpak_remote.env`).

### Plazo global de preparación (presupuesto predecible)

- Límite por fase sin presupuesto común: `update` 300 s × 2 rebotes + descarga 600 s × 2 + instalación 300 s × 2 = 3.000 s, por encima del `timeout-minutes: 45` del job. Por eso el preparador añade `prep_deadline_seconds=900` en `packaging/flatpak/toolchain.lock`: `phase_limit` usa el mínimo entre el límite de la fase y lo que queda de plazo, `deadline_exhausted` falla cerrado nombrando la fase, y `total` se mide desde el inicio del paso (antes del sondeo).
- Arnés con dos escenarios nuevos sobre el script real: J con `CHEMUSON_FLATPAK_PREP_DEADLINE_SECONDS=0` aborta al instante con `the 0s preparation deadline was exhausted before apt-get update` y `exit 4`; K con plazo de 2 s y un `apt-get update` colgado (simulado con `sleep 30`) se corta a los 2 s, hace el rebote sin presupuesto restante y termina en `exit 3` a los ~2 s en vez de los 600 s que permitiría el límite de fase.
- Contrato estático ampliado: se exige que ninguna fase use su límite bruto (`timeout "${UPDATE_TIMEOUT}"`, `…DOWNLOAD_TIMEOUT…`, `…INSTALL_TIMEOUT…` no aparecen), que existan `phase_limit`/`budget_left`/`deadline_exhausted` y que `prep_deadline_seconds` sea numérico y menor que la mitad de los 45 min del job. Once escenarios del arnés en su resultado esperado y 13 contratos del archivo.
