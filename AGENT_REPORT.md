# AGENT_REPORT — cierre de higiene del repositorio

## Alcance y decisiones

- Rama `maintenance/repository-hygiene-closure`, creada desde `origin/main`
  `1db4f63b52af79247745b3a8a220fb728348218c`, publicada mediante push normal.
- `origin/main` permanece en el mismo SHA. No se hizo merge, force-push,
  rebase, cherry-pick, `git gc --prune` ni reescritura de historia.
- OpenSpec: `openspec/changes/repository-hygiene-closure-2026-10-03/`;
  checklist completo. Memoria y decisiones en
  [`docs/history/CAMPAIGNS.md`](docs/history/CAMPAIGNS.md), política en
  [`docs/history/REPOSITORY_POLICY.md`](docs/history/REPOSITORY_POLICY.md) y
  métricas/refs/prune en
  [`docs/history/REPOSITORY_CLEANUP_2026-10-03.md`](docs/history/REPOSITORY_CLEANUP_2026-10-03.md).

## Cambios

- Se retiraron `src/sys`, `src/repro_v2.png`, dos parches sueltos ya históricos,
  124 outputs orbitales bajo `tests/archive`, prototipo UI ejecutable y
  capturas/generadores intermedios, scripts F7 one-shot y cuatro informes
  redundantes. OpenSpec Markdown/specs, tests y fixtures vigentes, cinco
  capturas UI finales, tres capturas KDE/Wayland, siete evidencias F7 selectas
  y el `theme.py` normativo se conservan.
- Se corrigieron enlaces absolutos `file://` del manual, se acortaron las
  fichas UI/known issues y comentarios de tema que apuntaban a prototipos
  retirados. Ningún archivo test de código ni módulo `src/chemuson/clean2d/`
  fue modificado; tampoco se cambiaron química, geometría, serialización o
  versión.
- Los previews orbitales ahora escriben fuera del árbol por defecto. El
  `render_orbital_family_preview.py` pasó smoke; el `orbital_fit_report.py`
  falla por un bug ya existente (`PiBondingParams.ring` inexistente en
  `_family_metric_strings`), dejando salidas parciales solo en temp. Se
  documenta y no se arregla fuera de alcance.
- Se borraron 23 ramas remotas y tres locales solo tras verificar SHA y
  ancestro de `origin/main`; quedaron intactas ocho ramas únicas/activas,
  incluida `gh-pages` y ambas Clean2D. Inventario final: 10 remotas reales,
  cuatro locales.

## Verificación

- `python -m compileall src tests tools packaging`: PASS.
- `pytest --collect-only -q`: 1816.
- `pytest -q tests/architecture`: 269 passed.
- UI dirigida: 299 passed.
- `pytest -q`: 1760 passed, 55 skipped, 1 fallo de baseline persistente:
  `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend`.
  Fue el único fallo en baseline y en ambos full runs del cierre; no se cambió.
- Ruff scoped `F401,F811,F821,E722,E741`: un F401 histórico en el test
  Clean2D `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py`; no se
  toca fuera de alcance.
- `openspec validate --all --strict`: 43 passed, 0 failed.
- Smoke Qt offscreen 100%/HiDPI 200%: PASS (temas claro/oscuro, 1440×900,
  980×600, DPR 2, tabs, rail/flyout, SidePanel, paleta y onboarding).
- `git diff --check`: PASS; escaneo de 350 Markdown: cero links relativos
  rotos. `git status --short`: limpio después del commit/push final.

Una primera invocación del smoke Qt olvidó `PYTHONPATH=src`; se repitió con
`PYTHONPATH` y `XDG_CONFIG_HOME` temporales y pasó, sin tocar preferencias del
usuario. Un primer `diff --check` encontró tres espacios finales en el baseline
OpenSpec; se corrigieron antes del push.

## Integridad final

Los tamaños, el alcance exacto de archivos, SHAs previos/posteriores y la
estimación de compactación histórica están en el reporte enlazado arriba. Los
objetos Flatpak activos publicados por `gh-pages` se preservan. No se inició
otra fase de producto.

## Continuación ChemUSON — Phase 4.5/5 (2026-10-07)

- Rama `ai/molecular-assistant-foundation`, base `6d4fca961cc7695ed01d78ae8f0a1ad769388384`. `c0233eb` (funcional) y `6d36e3a` (higiene/dependencias) se publicaron con pushes normales. No hubo merge, rebase ni force-push. El commit de cierre de este informe/tasks se publicará también mediante push normal.
- Se cerró el contrato tipado de transformación M23, Insert Variant/Replace undoable, controles y perfiles runtime sin secretos, verificación conservadora de identidad con ChemIO aislado y resumen seguro del evaluador Clean2D. No se modificó `src/chemuson/clean2d/`, persistencia `.cmsn` ni dependencias runtime.
- Tests: colección actual 1946; shards cubrieron 1922 passed, 20 skipped, 4 fallidos. Los fallos son el test Clean2D ya registrado en baseline y tres aserciones de importación estereoquímica bajo RDKit 2026.03.6. Suite de arquitectura 278 passed; OpenSpec estricto 52 passed. Ruff de archivos cambiados PASS; Ruff global conserva solo el F401 `math` histórico.
- No se repitió la suite monolítica: su referencia histórica es 19:26 y el límite activo es 10 minutos. El par CompChem→Assistant pasó 5/5 en ambos órdenes; el bloque ordenado Assistant→transform→UI pasó 99. No se reprodujo SIGSEGV en los shards, pero tampoco se afirma que el aborto monolítico histórico esté descartado ni se encontró un segundo trigger de producción.
- Smoke local únicamente contra Qwen 3.8 27B en `127.0.0.1:8081`: etanol/cafeína pasaron ChemIO e identidad offline; colesterol, vancomicina y eritromicina devolvieron `invalid_json`. No se usó red externa ni se guardó salida cruda. Evidencia detallada: `openspec/changes/stabilize-ai-molecular-assistant-integration/validation.md` y `/tmp/chemuson-qwen-local-smoke/`.

## Cierre documental y autorización de integración Molecular Assistant (2026-10-07)

La rama `ai/molecular-assistant-foundation` parte de `b4e0e661c0062e5e7bb69d5ea0dd712ea7e8e38e`, con `origin/main` en `4068319e9a1deee7dbf239872191a8f509c8b641`, 23 commits ahead y 0 behind; el árbol estaba limpio y el fast-forward era posible. OpenSpec estricto, arquitectura (280), compileall y diff-check pasaron; no hay diff de Clean2D. El propietario confirmó que `origin/main:pyproject.toml` ya declara Pillow y que `6d36e3a` reconcilió `requirements.txt` con esa declaración; conservar Pillow no constituye una dependencia nueva del proyecto. No se modifican dependencias ni código de producto. El bloqueo queda levantado y la microfase continúa con los gates de integración autorizados.

Como captura baseline obligatoria, se ejecutó una sola vez `timeout 10m pytest -q` antes de los cambios documentales: emitió dos marcadores `F` y abortó cerca del 63% en teardown Qt con SIGSEGV/exit 139 (`QUndoStack`/`QWidget`, worker `chemuson-3d_0`), sin resumen final. Se documenta como compatible con la deuda Qt/SIGSEGV histórica ya conocida; no se repitió ni se afirma que esté resuelta. Ruff global volvió a encontrar únicamente el F401 histórico `math` en el test Clean2D. Evidencia y comandos registrados en `openspec/changes/integrate-ai-reference-structure-resolution/baseline.md` y `validation.md`.

## Preparación ChemUSON v0.3.0-beta.1 + Actions preview (2026-10-08)

### Decisiones y cambios

- Rama local `release/v0.3.0-beta.1-prep`, creada desde `origin/main` `8640bed18f0961ef9582f8376e98e9a33dbc3c01`. Se preparó `_version.py` en `0.3.0-beta.1` y AppStream con fecha `2026-10-08`; Pillow y `pyproject.toml` dinámico se conservan. No se creó tag.
- Se endurecieron `release.yml` (tag-only, validación tag/versión/AppStream/SHA/API, gate acotado del mismo SHA, permisos read-only por defecto, `contents: write` sólo en jobs de publicación, no overwrite, colisiones de artifacts bloqueadas y procedencia/checksums). Se mantienen los builders y verificaciones oficiales; `test.yml` no se modificó.
- Nuevo `build-preview.yml` (local) con `workflow_dispatch` y push filtrado a `release/**-prep`. Los jobs verifican el mismo SHA/version canónico y producen cuatro grupos Actions con nombres distintos de release, checksum SHA-256 y `preview-provenance.json`. `contents: read`; sin secrets, tags/releases, `git push`, manifests públicos, remoto Chemuson ni `gh-pages`. El preview Linux omite sidecars públicos; Flatpak usa rama local `preview-<sha>`.
- Inno ahora exige `CHEMUSON_VERSION`; el smoke Windows existente conserva su valor explícito. Los scripts release validan la pareja versión/canal; no se altera lógica química, Clean2D, Molecular Assistant, dependencias ni `.cmsn`.
- Documentación: `docs/release/VERSIONING_POLICY.md`, `0.3.0-beta.1.md`, `manual-acceptance-0.3.0.md` (70 casos `NOT TESTED`), `PREVIEW_BUILDS.md`, catálogo exacto de excepciones, README, guías Linux/Windows/hotfix, OpenSpec y esta validación.

### Verificación

- OpenSpec estricto: `openspec validate prepare-chemuson-0-3-0-beta-1 --strict` PASS.
- Arquitectura: 280 passed. Suites enfocadas release/version/updater/package/preview/workflow: 104 passed en 1.76s. Colecta total: 2072 tests.
- `python -m compileall -q src tests tools packaging`, Ruff scoped, `bash -n`, parseo de AppStream XML/Flatpak YAML/workflows y `git diff --check`: PASS.
- Ruff global sigue fallando sólo por el F401 histórico `math` de `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no se tocó. Suite monolítica no ejecutada. `actionlint`, PyInstaller, flatpak-builder, appimagetool e Inno Setup no están disponibles localmente.
- `gh auth status`: sin sesión. El workflow no está en GitHub; no se inició un run. Builds reales: NOT BUILT; los smoke locales usan un ejecutable stub. Matriz manual: no ejecutada.
- No se realizaron operaciones de tag, Release, canal, remoto Flatpak o `gh-pages`; no hay cambios de Clean2D, dependencias ni persistencia. No hubo commit, push, merge, rebase ni force-push. El propietario debe revisar el árbol local y autorizar por separado cualquier push/ejecución.

### Estado

- `PREVIEW BUILD INFRASTRUCTURE: NOT READY` para ejecución remota: contratos estáticos locales pasan, pero la rama/workflow no está publicada ni probada en runners.
- `PREVIEW ARTIFACTS: NOT BUILT` (ninguno de los cuatro paquetes se compiló realmente).
- `RELEASE READINESS: NOT READY TO PUBLISH`; `MERGE READINESS: NOT READY`. Pendientes: disponibilidad autorizada del workflow en GitHub, builds reales, aceptación manual, revisión Qt/documentación y decisión del propietario.

## Addendum final: AppImage Type 2, CI PyYAML y P1 SVG (2026-10-08)

Este addendum supersede únicamente los estados de “builds reales NOT BUILT” y “test.yml no se modificó” de la captura anterior; los demás límites de publicación/aceptación continúan vigentes.

- Packaging local en `release/v0.3.0-beta.1-prep`: `ca7f89e` (`fix(packaging): build Type 2 AppImages and retain SVG assets`), `ec20a70` (`fix(ci): install development test dependencies`) y `46af127` (`fix(packaging): validate AppStream metadata for AppImage`). El spec común recoge los 69 SVG explícitamente; el smoke congelado valida QtSvg/`sys._MEIPASS`, DPR 2, iconos de ambos temas y dibujo visible en QToolButton. Preview/release Linux son AppImage Type 2 auténticos, no el PyInstaller renombrado. El primero pasó sin updater público; el release local validó update-information, `.updateinfo`, `.update.json` y `.zsync`. No se publicaron.
- Causa del P1 confirmada en el mismo spec usado por Windows/Linux: `collect_all("chemuson")` no detectaba los recursos porque ChemUSON no está instalado como paquete durante PyInstaller; el archive baseline tenía cero SVG. Los casos históricos Windows/Linux permanecen FAILED y el retest visual propietario está pendiente; el smoke no es aceptación gráfica.
- AppImageKit inicialmente falló la validación AppStream por developer ID con mayúsculas; se normalizó el ID y se usó `.appdata.xml`, con validación exitosa posterior. El único aviso AppStream final fue `redundante: 1` (exit 0). PyInstaller también reportó dependencias Qt opcionales ausentes en este host; los smoke/startup requeridos pasaron y no se afirma cobertura de esos módulos.
- CI: `test.yml` instala requirements runtime + dev y el editable, sin PyYAML duplicado ni dependencias Python runtime nuevas. No hay ejecución remota verificada.
- Verificación: 280 arquitectura; 42 pruebas packaging/workflow/CI; 43 UI SVG + versión; 52 versión/updater; 2086 recolectados; compileall, OpenSpec estricto, Ruff scoped, shell/YAML/XML/AppStream y diff checks pasaron. Ruff global conserva sólo el F401 histórico `math` en el test Clean2D. No se corrió la suite monolítica.
- AppImage preview local SHA-256 `c4d03edffa289118245c45bd165bff8cfe49ba5d171392b4297d05ed1beb21bb`; release local SHA-256 `688ac7bd866a676b8d0a78d2b9cdb3ae45f29d874d087f15e469abf4b25a8d11`; ambos se validaron contra source SHA `46af127e6790d62b713453475eabedf4513f7afb`. Evidencia detallada: `openspec/changes/prepare-chemuson-0-3-0-beta-1/{baseline,validation}.md`.
- `gh auth status` sigue sin sesión. El único `git push` normal no interactivo falló antes de actualizar remoto (`could not read Username`); no se usó contraseña. Push por GitHub Desktop/sesión normal autenticada y run de Actions quedan pendientes. Windows, artifacts remotos y retest visual no verificados. No se creó tag/Release, no se modificó `main`, `gh-pages` ni los canales públicos. El P1 bloquea publicación beta; deuda Qt/SIGSEGV no resuelta.
