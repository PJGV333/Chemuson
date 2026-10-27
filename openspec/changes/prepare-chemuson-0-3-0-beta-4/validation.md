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

### Revisión del diff: tres defectos propios corregidos

- `stage_cached_debs` y la publicación de `.deb` incrementaban su contador aunque `cp` fallara (`|| true`), de modo que el resumen del paso habría afirmado una caché publicada inexistente. Ahora cada contador comprueba el archivo de destino y, si ningún `.deb` llega a la caché, se registra `WARNING: N downloaded archives could not be published …; the next run will start with a cold cache.` con el recuento `0/N`. Escenario L del arnés (caché en `0555`): `published for reuse: 0/3`, aviso presente e instalación offline igualmente correcta (`exit 0`).
- Tras dos intentos de descarga sin presupuesto restante el script continuaba hacia la instalación offline, que sólo podía fallar de forma confusa. Ahora `download_ok=false` con `budget_left=0` invoca `deadline_exhausted "the package download"`, y el aviso de rebote informa el límite efectivo (`LAST_DOWNLOAD_LIMIT`) en lugar del límite teórico de 600 s.
- `APT_OPTS` fijaba reintentos y timeouts pero no `DPkg::Lock::Timeout`, el fallo típico cuando un runner está actualizando paquetes en segundo plano; se añadió `apt_lock_timeout_seconds=60` al lockfile, con validación numérica y contrato estático.
- Auditoría del manifiesto de node IDs tras los cambios: colección 2.666, manifiesto 2.644; los 22 node IDs fuera del manifiesto siguen siendo preexistentes y pertenecen a `tests/test_preview_workflow_contract.py` (3), `tests/test_release_tag_validation.py` (6) y los cuatro archivos de actualizaciones (13); ninguno proviene de archivos creados o modificados en esta campaña. Cero node IDs del manifiesto sin colección.

### Revisión del builder: fases cortadas con diagnóstico

- `build-update-repo` (300 s) y `build-bundle` (300 s) tenían límite de reloj pero no traducían la salida 124: un origen lento habría terminado el job con un código de salida sin causa, justo el tipo de evidencia que faltó en el run cancelado. Se añadió `run_phase`, que mide la fase, captura el código y reporta `exceeded its Ns wall-clock limit after Xs` o `failed (exit N)`, conservando el código original.
- Comprobado con el mock (`FAKE_REPO_SLEEP`/`FAKE_BUNDLE_SLEEP` y límites de 1 s): `ERROR: build-update-repo exceeded its 1s wall-clock limit after 1s.` y `ERROR: build-bundle exceeded its 1s wall-clock limit after 1s.`, ambos con `exit 124` y cortados a 1 s. La ruta feliz sigue generando `Chemuson-v0.3.0-beta.4-preview-<sha8>-linux-x86_64.flatpak` con `exit 0` y su tabla de fases en `$GITHUB_STEP_SUMMARY`.

### Ensayo integrador del job `Linux Flatpak bundle` (sin GitHub)

- Script local (`/tmp/rehearse_flatpak_job.sh`) que extrae del `build-preview.yml` real los cuerpos `run:` de los pasos 3-5 y los ejecuta en orden con `apt-get`/`curl`/`dpkg-query` y `flatpak`/`flatpak-builder`/`ostree` simulados, respetando el `env` declarado en el propio workflow (`CHEMUSON_FLATPAK_APT_PACKAGES`, `CHEMUSON_FLATPAK_APT_CACHE_DIR` bajo `runner.temp`) y la ruta de la caché del paso 2.
- Resultado: paso 3 `exit 0` (origen oficial sano, `published for reuse: 3/3`, instalación offline, versiones verificadas), paso 4 `exit 0` (`Verified preview source 08f90ee2…`), paso 5 `exit 0` con bundle `Chemuson-v0.3.0-beta.4-preview-08f90ee2-linux-x86_64.flatpak`, rama OSTree `preview-08f90ee2` y `Validated and checksummed 1 preview package(s)`. El `preview-provenance.json` generado registra `version: 0.3.0-beta.4`, `git_sha` único, `build_type: preview`, `publication: false`, y su `sha256` coincide con el `sha256sum` calculado aparte sobre el bundle.
- El primer intento del ensayo falló con `Flatpak repo metadata was not generated correctly` porque un stub de versión del arnés tapaba al generador de metadatos del repo: no es un defecto del producto, pero confirma que la guarda de metadatos OSTree falla cerrado cuando el repositorio no quedó bien generado. También se comprobó que la verificación de versiones aborta con `could not read the version of flatpak` si la CLI no responde (`exit 4`).

## Evidencia en GitHub Actions: Build Preview #41 en verde (2026-10-26 20:39–20:45Z)

- Run `40397595176` (Build Preview #41) sobre `08f90ee22c78cc451d7e7b363dc891f673c10f18`, disparado por push a `release/v0.3.0-beta.4-prep`. **7/7 jobs en `success`**: `Validate preview source`, `Linux Flatpak bundle`, `Windows portable and installer`, `Linux AppImage Type 2`, `Release preflight rehearsal (non-publishing)`, `Release gate bounded checks rehearsal` y `Preview build report`.
- Job `Linux Flatpak bundle` completo en **5,78 min** dentro de `timeout-minutes: 45` (sin cambiar):
  | paso | duracion |
  | --- | --- |
  | Restore verified Flatpak toolchain package cache | 0,00 min (cachefria) |
  | Install Flatpak build tools | **0,55 min** |
  | Verify exact source identity | 0,00 min |
  | Build local Flatpak bundle only | 5,03 min |
  | Upload Linux Flatpak preview | 0,07 min |
- Comparativa con el run cancelado #39 (`670a7670`): `Install Flatpak build tools` 43 min 59 s y `Build local Flatpak bundle only` cancelado a los 2 min 24 s. La preparacion bajo de 43 min 59 s a 0,55 min (48x) y el empaquetado vuelve a su rango historico (3,7-5,0 min). El run #38 habia quedado en 37 min 24 s de preparacion.
- Artefactos del run (todos del mismo SHA `08f90ee2`), con tamano y caducidad 2026-11-09: `chemuson-preview-windows-portable` 126.174.095 B, `chemuson-preview-windows-installer` 128.360.242 B, `chemuson-preview-linux-appimage` 171.550.885 B, `chemuson-preview-linux-flatpak` 106.949.974 B y `chemuson-preview-build-report` 461 B.
- La caché de `.deb` se restaura por prefijo (`restore-keys`), asi que este primer run fue frio (0,55 min de preparacion); el siguiente push medira la caché caliente.
- Verificacion de contenido de los paquetes aun pendiente de credenciales: la descarga sin sesion falla (`401` en la API, `404` en la URL web). Herramienta local lista (`objects` OSTree indexados por SHA-256 de contenido): el SVG empaquetado `packaging/flatpak/io.github.PJGV333.Chemuson.svg` en beta.4 vale `c7d017240cd89050…` y en beta.3 valia `23cad7927c43088f…`, y `assets/branding/logo_ChemUSON.svg` comparte el hash de beta.4, de modo que la presencia de ese objeto en el bundle prueba el logotipo nuevo.
- Consecuencia inmediata del fix del parser: el run anterior (#40) tuvo 0 jobs porque GitHub rechazo el archivo; este run si creo los 7 jobs, lo que confirma la correccion de los contextos de paso.

## Reproduccion local de los controles de identidad y del contrato de paquetes

Ejecutado en un worktree detached exacto en `08f90ee2` para que las guardas de
SHA coincidan con el candidato publicado (`git worktree add --detach`).

- `python packaging/release/preview_build.py prepare` con el entorno del job
  `Validate preview source` (`GITHUB_REF_TYPE=branch`, `GITHUB_REF_NAME=release/v0.3.0-beta.4-prep`,
  `GITHUB_SHA=08f90ee2…`): salida `Prepared preview 0.3.0-beta.4 …` y outputs
  `branch`, `sha=08f90ee2…`, `short_sha=08f90ee2`, `version=0.3.0-beta.4`.
  La guarda falla cerrado si el SHA esperado no es el del checkout: con el HEAD
  local (dos commits de documentacion por delante) aborta con
  `Checked-out preview source does not match the expected Git SHA`.
- Contrato de los cuatro paquetes del preview, con los nombres exactos que generan
  los jobs y contenido sintetico: `manifest --output-dir …` produce
  `checksums.sha256` que pasa `sha256sum -c` (4/4) y `preview-provenance.json`
  con `version: 0.3.0-beta.4`, `git_sha: 08f90ee2…`, `build_type: preview`,
  `publication: false`. Nombres:
  `Chemuson-v0.3.0-beta.4-preview-08f90ee2-windows-x86_64-portable.exe`,
  `…-windows-x86_64-setup.exe`, `…-linux-x86_64.AppImage`, `…-linux-x86_64.flatpak`.
- El manifiesto tambien falla cerrado si la salida queda fuera del repositorio
  (`Preview artifacts must be staged inside the checked-out repository`).
- Estas pruebas verifican el contrato de nombres/checksums/procedencia, no el
  contenido de los paquetes que genero Actions: esa verificacion requiere
  descargar los artefactos del run `40397595176`, y la descarga anonymous devuelve
  `401` en la API y `404` en la URL web; los resumenes de paso tampoco son
  accesibles sin sesion (la pagina publica del job solo expone los nombres).
- Referencia para verificar el logotipo dentro del bundle: el SVG empaquetado en el
  candidato (`packaging/flatpak/io.github.PJGV333.Chemuson.svg`, 2.326.332 B) vale
  `c7d017240cd890502157c7d3…` y es idéntico a `assets/branding/logo_ChemUSON.svg`;
  en `v0.3.0-beta.3` el mismo archivo valia `23cad7927c43088f…`.

## Logs reales del job Flatpak del run #41 (prueba de la reparacion)

Extraidos con `gh run view 40397595176 --log --job 121834716716`:

- Sonda de origenes oficiales: `azure.archive.ubuntu.com OK in 6084 ms (http)`,
  `archive.ubuntu.com OK in 689 ms`, `security.ubuntu.com OK in 691 ms`,
  `us.archive.ubuntu.com OK in 202 ms` -> `apt origin selected:
  https://us.archive.ubuntu.com/ubuntu/ (fastest of 4 healthy official mirrors)`.
  Es el caso exacto que provoco la cancelacion: el nodo prioridad 1 estaba lento y
  el script eligio un espejo oficial sano.
- Tiempos de preparacion declarados por el propio script:
  `probe=7 update=7 download=4 install=12 verify=0 total=32 deadline=900`
  (32 s frente a los 43 min 59 s del run cancelado).
- `cached .deb archives published for reuse: 123/123`, y el guardado de la cache
  ocup6 `98.956.656` B con clave
  `flatpak-apt-v1-preview-1e325240f9c956e5a22d2fa261ed9a30269171fc3770fd400b47624df77a9601-40397595176`
  (99 MB y 123 paquetes: coincide con el cierre historico que atasco al nodo throttled).
- Versiones verificadas contra el lockfile: `flatpak=1.14.6-1ubuntu0.1`,
  `flatpak-builder=1.4.2-1build2`, `ostree=2024.5-1build2`.
- Empaquetado real: `builder=139 repo_update=0 bundle=162 packaging_total=301`
  (5,02 min) y `Built dist-flatpak-preview/Chemuson-v0.3.0-beta.4-preview-08f90ee2-linux-x86_64.flatpak`.
- El primer run con la cache nueva la encontro ausente (`Cache not found for input keys: …`),
  como se esperaba; el guardado posterior deja el prefijo disponible para los siguientes runs.

## Verificacion de los cuatro artefactos del run #41 (descarga con sesion)

Descarga con `gh run download 40397595176` en `/tmp/paquetes41`. Los cuatro
paquetes comparten el SHA `08f90ee22c78cc451d7e7b363dc891f673c10f18` y la
`version: 0.3.0-beta.4` en su procedencia, con `build_type: preview` y
`publication: false`.

| paquete | SHA-256 verificado con `sha256sum -c` | tamano | runner |
| --- | --- | --- | --- |
| `Chemuson-v0.3.0-beta.4-preview-08f90ee2-linux-x86_64.flatpak` | `c8a9272edb4af9c2…` coincide | 106.918.304 B | ubuntu24 `20261023.334.3` |
| `Chemuson-v0.3.0-beta.4-preview-08f90ee2-linux-x86_64.AppImage` | coincide con `checksums.sha256` | 171.550.885 B | ubuntu24 |
| `Chemuson-v0.3.0-beta.4-preview-08f90ee2-windows-x86_64-portable.exe` | `259dc6cc856ea6c1…` coincide | 127.457.809 B | win25-vs2026 `20261021.278.1` |
| `Chemuson-v0.3.0-beta.4-preview-08f90ee2-windows-x86_64-setup.exe` | `d0fdb240cecbce48…` coincide | 128.826.659 B | win25-vs2026 |

En cada artefacto el `preview-provenance.json` tambien esta listado en
`checksums.sha256` y su suma coincide, de modo que el par paquete+procedencia
queda atado por SHA-256.

## Contenido real del bundle Flatpak (prueba byte a byte)

El bundle se abrio con la herramienta oficial en un entorno aislado
(`XDG_DATA_HOME=/tmp/fp-verify/data`, despues borrado; la instalacion Flatpak del
usuario no se modifico: `flatpak list --user` no muestra Chemuson y el repo del
usuario tiene 0 referencias de Chemuson):

- Ref OSTree del bundle: `app/io.github.PJGV333.Chemuson/x86_64/preview-08f90ee2`,
  commit `7888014cb3e3b560c010c04fd4b14daceab08e30f3060bbea56378897c77b85b`,
  ContentChecksum `da9b53a82726ec8cebec51e4a744e00305b89046684f5a9ea11d30c7dca34bb2`,
  fecha 2026-10-26 20:42:37 UTC.
- Logotipo nuevo: `/files/share/icons/hicolor/scalable/apps/io.github.PJGV333.Chemuson.svg`
  vale `c7d017240cd890502157c7d32cfbb1ea673f7ebf65b3651fef27fa7dc46e2d2f`, identico
  byte a byte (`cmp`) al SVG del candidato; el SVG de `v0.3.0-beta.3` valia
  `23cad7927c43088f7a01ae50bdf0a26dc5f15238cc157aa18228f5316afddb9e` y no es el que
  esta dentro. El PNG nuevo `/files/.../chemuson/assets/branding/chemuson-icon-512.png`
  coincide con el del repositorio (`d131782bfc04a9dc…`).
- Version declarada dentro del paquete: `_version.py` empaquetado contiene
  `__version__ = "0.3.0-beta.4"` y existe `chemuson-0.3.0b4.dist-info`; el metainfo
  empaquetado declara `<release version="0.3.0-beta.4" date="2026-10-26" />`.
- Componentes cientificos intactos en `/app/lib/python3.13/site-packages`:
  `rdkit` 2025.9.5 (+`rdkit.libs`, `rdkit-stubs`), `numpy` 2.3.2 (+`numpy.libs`),
  `PyQt6` 6.10.0 con `pyqt6_qt6` 6.10.0 y `pyqt6_sip` 13.10.2, y los modulos
  `chemuson/clean2d` (16 entradas), `chemuson/chemname` (21), `chemuson/gui` (58),
  `chemuson/chemio` (9).
- El log del job muestra el humo de ChemName ejecutado contra el arbol empaquetado:
  `"ok": true` con ethanol/benzene/acetamide/ethane/cyclohexane/unsupported_element
  en `pass` y las plantillas `.mol` verificadas por SHA-256 en `/app`.

## Run #42: la cache de paquetes funciona en CI (cachefria y calientefria)

Run `40399880083` (Build Preview #42) sobre `8d0fff77`, 7/7 jobs en verde.

- `Cache restored from key: flatpak-apt-v1-preview-1e325240…-40397595176`: acierto
  por prefijo `restore-keys` sobre la cache guardada por el run anterior (99 MB).
- Preparacion: `probe=7 update=6 download=1 install=10 verify=0 total=28`
  (28 s; la fase de descarga baja de 4 s a 1 s gracias a los `.deb` restaurados).
- `cached .deb archives published for reuse: 123/123` y guardado nuevo con clave
  `flatpak-apt-v1-preview-1e325240…-40399880083`.
- Empaquetado: `builder=134 repo_update=0 bundle=140 packaging_total=274` (4,57 min).
- Origen elegido de nuevo `https://us.archive.ubuntu.com/ubuntu/`.

## Correcciones P1/P2 — validación local previa al commit

- Base exacta `a95cf44969c938fabe10eca78450b515b2bcfa74`, con el diff de seguridad de esta revisión en el árbol; no son resultados de Actions. Python 3.14.8, pytest 9.1.1, Ruff 0.17.1, CachyOS; apt no está disponible localmente.
- P1: `apt-cache dumpavail` proporciona registros de paquetes de índices autenticados; tras `apt-get update --error-on=any` se comparan SHA-256 **y nombre** del archivo copiado a un directorio privado de root. El directorio de archivos global no se reutiliza. Se descartan archivos no coincidentes y apt descarga los restantes; si falla la adquisición no se instala.
- P2: al cambiar de origen se actualizan los índices, se comprueban sus destinos, se vuelve a validar la caché y sólo entonces se reintenta. Los límites usan el presupuesto global. El caso sin origen válido falla antes de instalar.
- `python -m pytest tests/test_flatpak_toolchain_preparation.py tests/test_preview_workflow_contract.py tests/test_release_workflow_contract.py tests/test_linux_distribution_scripts.py tests/test_appimage_validation.py tests/architecture -q`: 345 aprobadas; incluye siete escenarios del script real con apt simulado: archivo íntegro, alteración de igual tamaño, hash auténtico de otro paquete, rebote, ausencia de mirror, origen incorrecto y hashes ausentes.
- `python -m compileall -q src tests tools packaging`, Ruff focal `--select F401,F811,F821,E722,E741`, `bash -n` de ambos scripts Flatpak, parseo YAML de preview/release/test, `openspec validate prepare-chemuson-0-3-0-beta-4 --strict` y `git diff --check`: correctos.
- `pytest --collect-only -q`: 2673 node IDs; manifiesto 2651, siete nuevos incluidos, sin node IDs eliminados. Las 22 discrepancias previas no se modificaron. Sin suite masiva.
- Pendiente de comprobar con apt real de Ubuntu en el siguiente Build Preview y de nueva revisión de Codex; no se afirma ni publicación ni aprobación.
