# Baseline — preparación v0.3.0-beta.1

Fecha de captura: 2026-10-08 UTC. Capturada desde la rama nueva `release/v0.3.0-beta.1-prep`, creada limpia desde `origin/main`.

## Git y publicaciones observadas

- `git fetch origin --prune`: completado.
- `git status --short --branch` antes de crear la rama: `## main...origin/main` (limpio).
- `origin/main` y HEAD inicial: `8640bed18f0961ef9582f8376e98e9a33dbc3c01`; coincide con el HEAD de referencia recibido.
- `origin/gh-pages`: `c856fc1445a11576a5a2718cd40380bc54056853`.
- Tags de versión observados: `v0.2.5` es el más reciente; apunta a `ec5f9b811c9d7a71178d20edf49b4435286c7659`. No existe tag `v0.3.0-beta.1`.
- GitHub Releases API pública: `v0.2.5` stable, publicada 2026-05-01; `v0.2.4` stable; la beta más reciente es `v0.2.3-beta.3` (2026-04-09). No se consultó GitHub con credenciales ni se cambió una publicación.
- `https://pjgv333.github.io/Chemuson/`, y los `.flatpakref` beta/stable respondieron HTTP 200; ambos canales existen actualmente.
- Se creó la rama `release/v0.3.0-beta.1-prep` desde `origin/main`. HEAD inicial de rama: `8640bed18f0961ef9582f8376e98e9a33dbc3c01`.

## Versionado y pipeline antes de cambios

- `src/chemuson/_version.py`: `0.3.0-dev`.
- `pyproject.toml` usa `dynamic = ["version"]` y `chemuson._version.__version__` como fuente del paquete.
- AppStream metainfo conserva `0.2.3-beta.3` y `0.2.1`, por lo que no refleja la última estable ni la próxima beta.
- Inno Setup ya toma `CHEMUSON_VERSION`, pero conserva un fallback `0.0.0-dev` si falta.
- `release.yml` permite dispatch con versión por defecto obsoleta `0.2.3-beta.3` y canal elegible por separado. Un beta podría seleccionarse como stable; RC se clasifica actualmente como stable. Cada build modifica `_version.py` durante CI, de modo que los bytes del artefacto no corresponden literalmente al contenido fuente del tag.
- `release.yml` no tiene gate de pruebas en el mismo SHA; `test.yml` y release son workflows independientes. El workflow de release dispone de `contents: write` global y puede publicar en `gh-pages`.
- Checksums SHA-256 se generan siempre; HMAC y firma GPG Flatpak son opcionales por secretos. La firma no debe confundirse con el checksum.
- `build_appimage.sh` copia el ejecutable PyInstaller portable y lo nombra `.AppImage`; el reporte de QA anterior confirma que no es un contenedor AppImage Type 2. Se conserva como limitación explícita; no se afirma haber construido un AppImage real.
- Flatpak usa las ramas separadas `beta`/`stable`, App ID estable y runtime KDE 6.10; el manifiesto concede `--share=network` y `--filesystem=home`, sin que esta campaña amplíe permisos.
- No hay `docs/release/` ni política de versionado/aceptación estable.

## Baseline de herramientas y pruebas

- `python -m compileall src tests tools packaging`: exit 0.
- `pytest --collect-only -q`: `2027 tests collected in 0.93s`, exit 0; salida íntegra en `/tmp/chemuson-release-prep-baseline-collect.txt`.
- `timeout 8m pytest -q tests/architecture`: `280 passed in 12.50s`.
- `timeout 8m pytest -q tests/test_release_version_script.py tests/test_version_metadata.py tests/test_update_semver.py tests/test_update_policy.py tests/test_update_core.py tests/test_update_provider.py tests/test_update_security.py tests/test_update_portable.py tests/test_update_windows.py tests/test_update_telemetry.py tests/test_update_ui_text.py`: `52 passed in 0.68s`.
- `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`: exit 1 únicamente por el F401 histórico `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.
- XML AppStream y YAML Flatpak parsean correctamente.
- Herramientas locales: `flatpak`, `wine`, Python/pytest/ruff disponibles; `flatpak-builder`, `appimagetool`, `pyinstaller`, Inno Setup (`iscc`), `actionlint` y `shellcheck` no disponibles. No puede hacerse aquí una compilación real de los artefactos de los tres sistemas.
- `gh release list` no puede autenticarse con `gh`; se usó la API REST pública de GitHub sólo para lectura.

## Suite completa y excepciones históricas

No se repite la suite monolítica. En la revisión inmediatamente anterior, sobre el mismo árbol de producto (el commit actual `8640bed` añade cierre documental a `b4e0e66`), la ejecución baseline `timeout 10m pytest -q` mostró dos marcadores `F` antes de abortar cerca del 63% con SIGSEGV durante teardown Qt (`QUndoStack`/`QWidget`), exit 139 y sin resumen. No se identificaron los dos node IDs a partir de ese proceso; no se atribuyen causas. La deuda `_DescriptorWorker`/Qt sigue abierta. Evidencia previa más detallada: `openspec/changes/stabilize-gui-async-worker-shutdown/validation.md` y `openspec/changes/integrate-ai-reference-structure-resolution/baseline.md`.

Fallos de test conocidos identificados por campañas anteriores, pero no convertidos en skips globales: `tests/test_clean2d_engine_candidates.py::test_generate_candidates_attempts_rdkit_for_cyclic_graphs` (reproducido aisladamente en main); `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend` (fallo histórico de baseline); tres aserciones de `tests/test_smiles_stereo_import.py` bajo RDKit 2026.03.6 no declaradas baseline concluyente. Ninguno se corrige aquí. La política de excepciones exactas y la evidencia se detallan en `docs/release/KNOWN_BASELINE_EXCEPTIONS.json`.

## Baseline del requisito de preview

En el HEAD inicial `8640bed18f0961ef9582f8376e98e9a33dbc3c01` no existía `.github/workflows/build-preview.yml` ni un helper que publicara artifacts preview con procedencia. El único workflow era el release oficial, con trigger de tag y dispatch manual separado; `contents: write` era global. No existía un contrato estático de aislamiento preview. No se ejecutó ninguna build preview en la baseline, no se subieron artefactos y no se consultaron canales con credenciales.

## Baseline del addendum Type 2 AppImage / CI (2026-10-08)

Commit de referencia y HEAD limpio: `f0fde603371255902bf0e630ca1ce032b8f16ad3`; rama `release/v0.3.0-beta.1-prep`. `git status --short --branch` devolvió `## release/v0.3.0-beta.1-prep...origin/main` sin cambios. Se creó checkpoint local `checkpoint/beta1-before-type2-appimage-and-ci` en ese SHA; no es un tag ni una publicación.

- `python -m compileall -q src tests tools packaging`: exit 0.
- `timeout 8m python -m pytest --collect-only -q`: `2072 tests collected in 0.67s`, exit 0; salida en `/tmp/chemuson-beta1-appimage-collect-baseline.txt`.
- `timeout 8m python -m pytest -q tests/architecture`: `280 passed in 13.27s`.
- Suites enfocadas release/version/updater/preview/workflows/Linux existentes: `104 passed in 1.82s`.
- Ruff global (`F401,F811,F821,E722,E741`): exit 1 sólo por el F401 histórico `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.
- No se ejecutó `pytest -q` monolítico por la deuda SIGSEGV de teardown Qt ya documentada y el límite explícito de no repetir suites monolíticas prolongadas. No hubo cambios de código antes de esta captura.
- CI observado en `.github/workflows/test.yml` instala `requirements.txt` y `pytest`, omite `requirements-dev.txt`; `requirements-dev.txt` ya declara `PyYAML`, y `requirements.txt` no lo declara. No se duplicará PyYAML ni se añadirá dependencia runtime.
- Reporte público de Actions consultado por API REST de sólo lectura: run `37826597134`, `f0fde603371255902bf0e630ca1ce032b8f16ad3`, conclusión `success`; cuatro grupos de paquetes y reporte presentes. Los jobs pasaron, pero el artefacto Linux de ese SHA es sólo el PyInstaller renombrado. No se modificó ni descargó/publicó ese run.
- `appimagetool` no estaba instalado localmente; PyInstaller tampoco. `desktop-file-validate` sí existe. Para fijar la fuente se consultó la API oficial `AppImage/AppImageKit`: asset x86_64 ID `98605504`, URL oficial del release `continuous`, commit reportado `5735cc5bed206497cddfbd2a75e1982c2606c35d`, actualizado 2025-07-26. Descarga aislada en `/tmp` validada con SHA-256 `b90f4a8b18967545fda78a445b27680a1642f1ef9488ced28b65398f2be7add2`; su propia firma de AppImage es `AI\\x02`, y el AppRun informa `appimagetool, continuous build (commit 5735cc5), build <local dev build> built on 2023-03-08 22:52:04 UTC`. Se ejecutó desde extracción explícita, sin FUSE.
- `gh auth status`: no autenticado. No se intentó push durante la captura baseline.

## Baseline del defecto P1 de iconos SVG empaquetados (2026-10-08)

Confirmación manual del propietario en artifacts del run `37826597134`: los iconos esenciales están ausentes en Windows portable y Linux portable, también en ambos temas; algunos iconos dinámicos de átomos/esfera sí aparecen. Los casos visuales quedan `FAILED — P1 blocks beta acceptance`; la aceptación no se convierte en PASS por pruebas automatizadas.

Reproducción local sobre el mismo `chemuson.spec`, antes de corregirlo:
- Build real: `timeout 9m bash -lc 'source /tmp/chemuson-appimage-venv/bin/activate && pyinstaller --clean --noconfirm chemuson.spec'`; PyInstaller 6.22.3, Linux x86_64; terminó correctamente.
- El log emitió `collect_data_files - skipping data collection for module 'chemuson' as it is not a package` (también `collect_dynamic_libs`). Los jobs Windows/Linux instalan requirements y PyInstaller, pero no instalan el proyecto como paquete; el spec sólo añade `src` a `pathex`.
- `pyi-archive_viewer -l dist/Chemuson`: **0 SVG estáticos `i-*.svg`**, y `PyQt6/QtSvg.abi3.so` + QtSvgWidgets sí están en el archive.
- `dist/Chemuson --version`, ejecutado desde `/tmp` y `/`, devolvió `0.3.0-beta.1`; confirma que cwd no impide arrancar, no que los iconos estén incluidos.
- `python -m pytest -q tests/test_ui_svg_icons.py`: `41 passed in 1.00s`; confirma rasterizado desde fuente, no empaquetado.
- Causa reproducida: `collect_all("chemuson")` no recolecta los SVG cuando el paquete no está instalado en el entorno del build; el binario compartido queda sin recursos estáticos. Explica el mismo defecto en Windows y Linux. Los iconos dinámicos son generados por código y no requieren esos SVG. `IconProvider._tinted_pixmap` reemplaza explícitamente `currentColor` antes de `QSvgRenderer`; la prueba fuente rasteriza con tinta en ambos temas.
- Verificación posterior a la corrección (Linux real): PyInstaller 6.22.3 archivó exactamente 69 SVG, el executable congelado reportó `QtSvg` y `IconProvider.__file__` bajo `sys._MEIPASS`, y renderizó todos los 69 iconos estáticos + iconos esenciales con píxeles visibles en temas claro/oscuro y DPR 2 desde un cwd `/tmp`. El AppImage Type 2 preview y release repitieron el smoke tras extracción. Windows sigue sin build real en este host; falta el run Actions.

## Verificación final del addendum (commit limpio `46af127e6790d62b713453475eabedf4513f7afb`)

- Entorno Linux local: Python 3.14.7, PyInstaller 6.22.3. `pyinstaller --clean --noconfirm chemuson.spec` terminó; el smoke congelado reportó `frozen=true`, 69/69 SVG visibles, QtSvg, rutas bajo `sys._MEIPASS`, cwd `/tmp`, DPR 2, 11/11 iconos esenciales y 11/11 renders `QToolButton` en claro y oscuro.
- La herramienta fijada por asset ID `98605504`, commit `5735cc5`, SHA-256 `b90f4a8b18967545fda78a445b27680a1642f1ef9488ced28b65398f2be7add2`, generó AppImages Type 2 auténticos. AppImageKit verificó AppStream y `desktop-file-validate`; `appstreamcli` devuelve exit 0 con una advertencia `redundante: 1`.
- Preview local: `dist-preview-final/Chemuson-v0.3.0-beta.1-preview-46af127e-linux-x86_64.AppImage`; SHA-256 `c4d03edffa289118245c45bd165bff8cfe49ba5d171392b4297d05ed1beb21bb`. Extracción sin FUSE, AppDir, versión, Qt/PyInstaller CArchive, smoke congelado y startup headless de 12 s: PASS. Update information y sidecars públicos ausentes.
- Release local (sin publicar): `dist-appimage-final/Chemuson-v0.3.0-beta.1-linux-x86_64.AppImage`; SHA-256 `688ac7bd866a676b8d0a78d2b9cdb3ae45f29d874d087f15e469abf4b25a8d11`. Extracción, smoke, versión y startup headless: PASS. `.updateinfo` 97 B, `.update.json` 475 B, `.zsync` 607312 B; información embebida: `gh-releases-zsync|PJGV333|Chemuson|prerelease|Chemuson-v0.3.0-beta.1-linux-x86_64.AppImage.zsync`.
- Fallo de packaging encontrado y corregido durante esta verificación: AppImageKit rechazó el ID de developer en mayúsculas `io.github.PJGV333` (`developer-id-invalid`) y no detectaba el nombre `.metainfo.xml`. Se normalizó el ID de developer a minúsculas, AppImage usa el nombre `.appdata.xml` esperado, y ambos jobs instalan `appstream`; la validación completa posterior pasó. Se eliminó también la categoría principal duplicada del desktop entry.
- PyInstaller dejó avisos no fatales de `collect_all("chemuson")` (se siguen incluyendo los recursos de iconos por `datas` explícitos) y de módulos Qt opcionales ausentes en este host (`QtStateMachine`, `QtSerialPort`, `QtSensors`, `QtRemoteObjects` y plugins de bases de datos/text-to-speech). Los smoke/version/startup de los binarios requeridos pasaron; esto no acredita las funciones Qt opcionales no ejercitadas ni sustituye la build de Actions.
