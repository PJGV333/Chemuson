# Validación — preparación v0.3.0-beta.1 y builds preview

Fecha: 2026-10-08. Rama de trabajo: `release/v0.3.0-beta.1-prep`. Commits de packaging local: `ca7f89e` (Type 2 + SVG) y `46af127e6790d62b713453475eabedf4513f7afb` (validación AppStream). Commit CI: `ec20a70`. Reproducción baseline y evidencia de paquetes: [`baseline.md`](baseline.md).

## Contratos, tests y análisis estático

- `openspec validate prepare-chemuson-0-3-0-beta-1 --strict`: PASS.
- `timeout 5m python -m pytest -q tests/architecture`: **280 passed in 13.48s**.
- `timeout 8m python -m pytest -q tests/test_appimage_validation.py tests/test_packaged_icon_smoke.py tests/test_linux_distribution_scripts.py tests/test_preview_workflow_contract.py tests/test_release_workflow_contract.py tests/test_release_artifacts.py tests/test_test_workflow_contract.py`: **42 passed in 1.20s**.
- `timeout 4m python -m pytest -q tests/test_ui_svg_icons.py tests/test_version_metadata.py`: **43 passed in 1.03s**.
- Versión/updater: **52 passed in 1.02s**.
- `timeout 5m python -m pytest --collect-only -q`: **2086 tests collected in 0.69s**.
- `python -m compileall -q src tests tools packaging`, `bash -n` de builders/helpers Linux, parseo YAML de los tres workflows y manifiesto Flatpak, parseo XML AppStream y `git diff --check`: PASS.
- `appstreamcli validate packaging/flatpak/io.github.PJGV333.Chemuson.metainfo.xml`: exit 0 (`redundante: 1`).
- Ruff scoped para módulos/tests cambiados (`F401,F811,F821,E722,E741`): PASS. Ruff global vuelve a fallar sólo por el F401 histórico `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no se modificó.
- No se ejecutó la suite monolítica ni se volvió a probar el aborto Qt/SIGSEGV histórico; no se declara resuelta esa deuda. Cada comando de tests respetó el máximo de 10 minutos.

## AppImage Type 2 y P1 de iconos

- Causa reproducida en baseline: `collect_all("chemuson")` omitía los SVG porque ChemUSON no se instala como paquete durante PyInstaller; el binario contenía cero `i-*.svg`, aunque QtSvg estaba incluido. El spec compartido ahora agrega explícitamente los 69 SVG y su licencia; el lookup se comprueba bajo `sys._MEIPASS` y desde un cwd aislado.
- En el commit de packaging `46af127e6790d62b713453475eabedf4513f7afb`, PyInstaller 6.22.3 / Python 3.14.7 generó el ejecutable Linux real. Smoke congelado: 69/69 SVG visibles, QtSvg, ruta `IconProvider` bajo `_MEIPASS`, DPR 2 y 11/11 iconos esenciales dibujados en `QToolButton` en temas claro y oscuro.
- Preview local —no Actions— `Chemuson-v0.3.0-beta.1-preview-46af127e-linux-x86_64.AppImage`, SHA-256 `c4d03edffa289118245c45bd165bff8cfe49ba5d171392b4297d05ed1beb21bb`: ELF x86_64 + `AI\x02`, extracción sin FUSE, AppDir/desktop/AppStream, CArchive, versión y smoke congelado: PASS. 69 SVG, ambos temas, QToolButton 11/11 y startup headless 12 s; sin update-information ni sidecars públicos.
- Release AppImage local (no publicado) `Chemuson-v0.3.0-beta.1-linux-x86_64.AppImage`, SHA-256 `688ac7bd866a676b8d0a78d2b9cdb3ae45f29d874d087f15e469abf4b25a8d11`: los mismos checks Type 2/extracción/iconos/versión/headless: PASS. `.updateinfo` 97 B, `.update.json` 475 B, `.zsync` 607312 B y update-information `gh-releases-zsync|PJGV333|Chemuson|prerelease|Chemuson-v0.3.0-beta.1-linux-x86_64.AppImage.zsync` verificados.
- AppImageKit oficial usado: asset ID `98605504`, commit `5735cc5`, SHA-256 `b90f4a8b18967545fda78a445b27680a1642f1ef9488ced28b65398f2be7add2`. AppImageKit/appstreamcli mostró `redundante: 1` como advertencia, con validación exitosa.
- Durante la primera validación AppImageKit sí falló por `developer-id-invalid io.github.PJGV333` y por no reconocer el sufijo `.metainfo.xml`. Se corrigió el developer ID de AppStream a minúsculas, AppImage usa el nombre `.appdata.xml`, se quitó la categoría desktop principal duplicada y ambos jobs instalan `appstream`; la validación posterior pasó. Sin cambios de publicación Flatpak/remoto.
- PyInstaller aún informa que `collect_all("chemuson")` omite el paquete (los SVG ya se recopilan por `datas` explícitos) y avisa de módulos Qt opcionales ausentes en este host (`QtStateMachine`, `QtSerialPort`, `QtSensors`, `QtRemoteObjects` y algunos plugins no usados). El smoke y startup de los binarios requeridos pasaron; no se afirma que esos módulos opcionales hayan sido probados.

## Aceptación manual y verificación remota

- El propietario confirmó los iconos ausentes en portable Windows/Linux del run `37826597134`; `UI-FAIL-WIN-01` y `UI-FAIL-LINUX-01` permanecen **FAILED — P1 blocks beta acceptance**. UI-07/UI-08 son retests separados y siguen `NOT TESTED`.
- El smoke automatizado, incluido el dibujo QToolButton, **no** cambia los casos manuales a PASS. El propietario debe instalar y revisar visualmente Windows portable y Linux portable/AppImage en ambos temas. Este P1 bloquea explícitamente publicación beta.
- No hay build Windows local, artifact corregido de Actions ni run remoto nuevo. El build Linux local no equivale a un preview verificado por GitHub Actions. Windows, installer/setup, Flatpak real y aceptación gráfica siguen sin verificar.
- `.github/workflows/test.yml` instala `requirements.txt`, `requirements-dev.txt` y el proyecto editable; no duplica PyYAML ni añade dependencias Python runtime. El contrato estático comprueba pytest real sin skips/enmascaramiento, pero aún no se ejecutó el workflow remoto.
- Al comenzar esta tarea, `gh auth status` no tenía una sesión activa y un intento anterior de push había fallado con `could not read Username`. Más tarde el CLI indicó que la sesión HTTPS normal del propietario estaba autenticada; el resultado del push de esta tarea se registra en la sección de validación incremental inferior.
- No se creó tag ni Release, no se modificaron `main`, `gh-pages`, canales públicos ni publicación remota. No se tocaron Clean2D, química, Molecular Assistant ni persistencia `.cmsn`; los cambios de GUI se limitan al onboarding solicitado y las cadenas visibles de marca.

## Dictamen previo a UI-ONBOARDING-001 / BRANDING-001

- `PREVIEW BUILD INFRASTRUCTURE: STATIC CONTRACTS PASS; REMOTE RUN BLOCKED BY MISSING AUTH` (estado al cerrar la validación previa).
- `LOCAL LINUX APPIMAGE TYPE 2: PASS` (preview y release locales, source SHA `46af127e6790d62b713453475eabedf4513f7afb`).
- `WINDOWS FROZEN BINARY / ACTIONS ARTIFACTS: NOT VERIFIED`.
- `MANUAL ICON ACCEPTANCE: FAILED / RETEST PENDING — P1 BLOCKS BETA PUBLICATION`.
- `RELEASE READINESS: NOT READY TO PUBLISH`; `MERGE READINESS: NOT READY`.

## UI-ONBOARDING-001 y BRANDING-001 — validación incremental (2026-10-08)

- `openspec validate prepare-chemuson-0-3-0-beta-1 --strict`: PASS; `git diff --check`: PASS.
- Tests focalizados de GUI, onboarding, marca, instaladores/metadatos, AppImage validator, workflows y actualización: **88 passed in 13.30s**.
- Tests de arquitectura: **280 passed in 10.78s**. Tests de updater y compatibilidad: **36 passed in 0.87s**. Persistencia/autosave/recovery: **34 passed in 0.99s**.
- Geometría del onboarding en escalas Qt simuladas `QT_SCALE_FACTOR=1.0, 1.25, 1.5, 2.0`: **3 tests passed en cada escala** (12 ejecuciones); las pruebas ejercitan 980×600, 1440×900, 1600×900, movimiento, agujeros de objetivos y límites de tarjeta.
- `pytest --collect-only -q`: **2093 tests collected in 0.98s**. No se ejecutó `pytest -q` monolítico por el antecedente Qt/SIGSEGV documentado.
- `ruff check` sobre todos los Python modificados, reglas `F401,F811,F821,E722,E741`: PASS. `python -m compileall src tests tools packaging`: PASS.
- `bash -n` en builders Linux, parseo YAML de Flatpak y workflows de release/preview, AppStream XML y `desktop-file-validate` de los desktop entries y la plantilla AppImage: PASS. AppStream devuelve exit 0 con `redundante: 1`; desktop-file-validate conserva un hint de categorías principales múltiples en el desktop Flatpak (sin error).
- Identidad técnica protegida por tests: versión `0.3.0-beta.1`, package/import, `QApplication.applicationName`, QSettings, `.cmsn`, App ID/repo, rutas/remote updater, ejecutables/nombres de artifact, AppId/instalación existente de Inno. Se agregó además la comprobación ChemUSON al validador AppImage y al índice HTML Flatpak.
- Aceptación manual `UI-ONBOARDING-001` y `BRANDING-001` sigue **NOT TESTED — owner retest required**. Escalas Qt simuladas y validación estática no acreditan paquete/instalación reales Windows/Linux. No se inició ni publicó ningún release/tag.
