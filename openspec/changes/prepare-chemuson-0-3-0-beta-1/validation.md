# Validación — preparación v0.3.0-beta.1 y builds preview

Fecha: 2026-10-08. Rama local: `release/v0.3.0-beta.1-prep`. Base del addendum: `f0fde603371255902bf0e630ca1ce032b8f16ad3`. Evidencia baseline completa y reproducción del defecto: [`baseline.md`](baseline.md).

## Contratos, tests y análisis estático

- `openspec validate prepare-chemuson-0-3-0-beta-1 --strict`: PASS.
- `timeout 5m python -m pytest -q tests/architecture`: **280 passed in 13.48s**.
- Validación enfocada AppImage/Linux/workflows/iconos: **42 passed in 1.28s**.
- `timeout 4m python -m pytest -q tests/test_ui_svg_icons.py`: **41 passed in 1.03s**.
- Versión/updater: **52 passed in 1.02s**.
- `timeout 5m python -m pytest --collect-only -q`: **2086 tests collected in 0.69s**.
- `python -m compileall -q src tests tools packaging`, `bash -n` de los builders/helpers Linux, parseo YAML de los tres workflows y manifiesto Flatpak, parseo XML AppStream y `git diff --check`: PASS.
- Ruff scoped sobre los módulos/tests Python añadidos o editados, con `F401,F811,F821,E722,E741`: PASS. Ruff global conserva únicamente el F401 histórico `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no se modificó.
- No se ejecutó la suite monolítica ni se volvió a probar el aborto Qt/SIGSEGV histórico; no se declara resuelta esa deuda. Cada comando de tests respetó el máximo de 10 minutos.

## AppImage Type 2 y validación de iconos

- La causa P1 se reprodujo en un ejecutable Linux PyInstaller 6.22.3 de baseline: `collect_all("chemuson")` omitió la colección porque el paquete no está instalado en el entorno de build; el archive contenía cero SVG estáticos `i-*.svg` aunque incluía QtSvg. El spec compartido se corrigió para incluir explícitamente los 69 SVG y su licencia, independientemente del cwd o instalación editable.
- Smoke de proceso congelado: valida las rutas bajo `sys._MEIPASS`, inventario completo, QtSvg, render SVG estático y esencial en ambos temas y DPR 2. También compara el dibujo de los iconos esenciales en un `QToolButton` con un control vacío. El smoke está opt-in y el validador ejecuta el binario desde un cwd aislado.
- Un build PyInstaller/AppImage Type 2 Linux previo pasó extracción sin FUSE, verificación ELF/`AI\x02`, inspección del AppDir, versión, sidecars/update-information cuando aplica y smoke headless acotado. Se repetirá el build final desde un commit limpio para validar la versión actual del smoke con `QToolButton`; esa evidencia final aún está pendiente en esta captura.
- El headless/offscreen sólo acredita carga/rasterizado automatizados; **no** acredita funcionamiento gráfico ni acepta el defecto visual.

## P1 y aceptación manual

- El propietario confirmó el defecto en los paquetes Windows y Linux portable del run `37826597134`. La matriz conserva esos casos como `FAILED — P1 blocks beta acceptance`; se añadieron casos de retest separados para Windows portable y Linux portable/AppImage.
- La corrección automatizada no cambia las filas manuales a PASS. **El propietario debe instalar y verificar visualmente los paquetes corregidos en ambos temas y plataformas.** Esa aceptación no se ha realizado y bloquea expresamente la publicación beta.
- Windows no se construyó en este host. El smoke fail-closed está conectado a los jobs Windows/Linux de preview y release, pero aún no se ha ejecutado en runners de Actions. No hay artifacts corregidos de GitHub ni run nuevo.
- Linux/Windows setup, Flatpak real y aceptación gráfica no se declaran verificados por los tests locales. AppImage Type 2 sí cuenta con build local, sujeto al rebuild final limpio indicado arriba.

## Integridad y publicación

- `.github/workflows/test.yml` instala `requirements.txt`, `requirements-dev.txt` y el proyecto editable, sin duplicar PyYAML ni añadir dependencias runtime; la prueba de contrato verifica que pytest real no se omite ni enmascara.
- `gh auth status`: **no autenticado**. No se hizo push ni se inició un Actions run. Si sigue así tras los commits, queda pendiente publicar sólo `release/v0.3.0-beta.1-prep` mediante GitHub Desktop o una sesión autenticada normal; no se intentará autenticación por contraseña.
- No se creó tag ni GitHub Release, no se modificó `main`, `gh-pages`, canales públicos ni remoto Flatpak. No se tocó Clean2D, química, Molecular Assistant, GUI de producto ni persistencia `.cmsn`.
- Pendiente: commits focalizados de packaging y CI; push normal de la rama si hay autorización; comprobar el preview/CI remoto, con Windows y Linux evaluados por separado; retest visual del propietario.

## Dictamen

- `PREVIEW BUILD INFRASTRUCTURE: STATIC CONTRACTS PASS; REMOTE RUN NOT VERIFIED`.
- `LOCAL LINUX APPIMAGE TYPE 2: PASS (build anterior); FINAL CLEAN-COMMIT REBUILD PENDING`.
- `WINDOWS FROZEN BINARY / ACTIONS ARTIFACTS: NOT VERIFIED`.
- `MANUAL ICON ACCEPTANCE: FAILED / RETEST PENDING — P1 BLOCKS BETA PUBLICATION`.
- `RELEASE READINESS: NOT READY TO PUBLISH`; `MERGE READINESS: NOT READY`.
