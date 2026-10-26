# Compilaciones de prueba en GitHub Actions

`Build Preview` compila paquetes para evaluación manual desde ramas `release/**-prep`. El resultado son artifacts temporales de una ejecución de Actions; **no** es una publicación ni una versión disponible en el updater.

## Opción A — Desde GitHub

1. Abrir el repositorio → **Actions** → **Build Preview**.
2. Si aparece **Run workflow**, elegir la rama `release/v0.3.0-beta.4-prep` (o la rama de preparación correspondiente). La versión se obtiene de `src/chemuson/_version.py`; no se escribe manualmente.
3. Ejecutar el workflow y esperar a que la validación de rama/SHA, los dos jobs de ensayo del release (`Release preflight rehearsal` y `Release gate bounded checks rehearsal`) y los tres jobs de paquetes terminen.
4. Windows y Linux ejecutan `validate_packaged_icons.py` y `validate_packaged_rdkit_worker.py` sobre el ejecutable PyInstaller real antes de subirlo. El gate RDKit comprueba que el CArchive no contenga recursos top-level `PIL`/`Pillow`, además de imports nativos dentro del bundle, descriptores conocidos de etanol, SMILES de entrada/canónico y conformero 3D mediante el worker aislado; no admite skip. La metadata de RDKit puede resolver Pillow transitivamente en el entorno de build, pero el spec excluye esos recursos del paquete. En Linux, el job construye AppImage Type 2 y `validate_appimage.py` repite ambos smokes sobre el ejecutable extraído del paquete, además del smoke esencial de flujo real 2D→3D con el worker RDKit congelado (`geometry_flow_smoke`), que exige `ok=true` con posiciones y energía reales. Estos gates automatizados no sustituyen la aceptación manual del propietario.
5. Revisar el resumen final. El run sólo es completo si el ensayo del preflight, el ensayo de los checks acotados del gate, Windows portable/installer, Linux AppImage Type 2 y Linux Flatpak muestran `success`.
6. Abrir la ejecución → sección **Artifacts** y descargar por separado:
   - `chemuson-preview-windows-portable`
   - `chemuson-preview-windows-installer`
   - `chemuson-preview-linux-appimage`
   - `chemuson-preview-linux-flatpak`
   - `chemuson-preview-build-report` (estados por plataforma)
7. En cada paquete, revisar `preview-provenance.json` y verificar los archivos listados en `checksums.sha256`. Confirmar versión, rama, SHA completo, sistema y `publication: false` antes de instalar. El AppImage debe devolver la misma versión con `./Chemuson-...AppImage --appimage-extract` en una carpeta temporal; no se requiere FUSE.
8. Probar en VM/perfil aislado y registrar los resultados en [la matriz manual](manual-acceptance-0.3.0.md). El smoke automatizado sí verifica que Qt rasteriza los iconos empaquetados, pero no sustituye inspección visual interactiva ni la aceptación del propietario. No usar el setup/Flatpak preview sobre un perfil de trabajo sin respaldo: los instaladores comparten el App ID de la aplicación.

**Disponibilidad del dispatch:** GitHub sólo permite `workflow_dispatch` cuando reconoce la definición del workflow en el repositorio. Si la opción no aparece mientras la campaña sigue en una rama de preparación y el workflow aún no está disponible desde la rama por defecto, no hacer merge ni cambiar permisos para forzar el dispatch. El push filtrado a `release/**-prep` es el disparador automático previsto una vez que el workflow se encuentre en GitHub; la rama debe publicarse mediante el procedimiento normal autorizado.

Desde la simplificación de la pipeline, `test.yml` ya no se dispara automáticamente en cada push: la suite completa corre sólo de forma extraordinaria vía `workflow_dispatch` (o localmente con `tools/ci_pytest_shards.py`), según [la política de pruebas](../testing/AGENT_TESTING_POLICY.md). La verificación ordinaria previa al push es local y proporcional al cambio. Los builds de preview nunca dependieron ni dependen de la suite masiva: construyen código ya validado y conservan sus gates ligeros de integridad, incluido el smoke esencial de flujo real 2D→3D en el AppImage. La deuda Qt/teardown histórica sigue documentada.

## Opción B — Desde un agente autenticado

Un agente sólo debe iniciar un run si el propietario autorizó la acción y GitHub informa que el workflow ya está disponible. Primero comprobar `gh auth status` y que la sesión tenga permiso de ejecutar workflows; no asumir que Luna tiene credenciales. Sin autorización, detenerse y dejar los paquetes como `NOT BUILT`.

```bash
gh workflow run build-preview.yml --ref release/v0.3.0-beta.4-prep
gh run list --workflow build-preview.yml --branch release/v0.3.0-beta.4-prep --limit 5
gh run view <run-id> --web
gh run download <run-id> --dir "chemuson-preview-<run-id>"
```

El agente debe devolver la URL del run, SHA/versión y resultado de cada job, además de indicar si se generaron los cuatro grupos. Un run parcial o fallido **no** se informa como build exitoso. No usar API/token con permisos adicionales ni otro workflow para eludir el filtro.

## Aislamiento y alcance real

- `GITHUB_TOKEN` tiene `contents: read`; el preview no recibe secrets de firma/publicación. El ensayo del preflight usa ese mismo token de lectura para la consulta de inexistencia de Release.
- El workflow no crea tags/releases públicos, no hace `git push`, no despliega `gh-pages`, no construye una URL pública ni actualiza canales Flatpak beta/stable o manifests del updater. Única excepción controlada: el job `Release preflight rehearsal` crea un tag anotado temporal `v<versión>` **solo dentro del runner** para ejercitar el validador oficial, verifica vía `git ls-remote` que el tag oficial no existe ni existió, y lo borra con `if: always()` al terminar. Nunca se publica.
- Linux es un AppImage Type 2 auténtico con AppRun/desktop/SVG/AppStream y recursos PyQt6/ChemUSON validados; el preview omite update-information embebido y `.updateinfo`, `.update.json`, `.zsync`.
- El Flatpak es un bundle local de rama `preview-<sha>`; no se produce un remoto público ni un `.flatpakref` actualizable.
- Los tests de contrato inspeccionan permisos, SHA, nombres y ausencia de pasos de publicación. Son evidencia **estática** únicamente; no prueban compilación real.
- La retención esperada de artifacts es 14 días. Descargarlos pronto si se necesita conservar evidencia.

## Ensayo del preflight oficial antes del tag (a partir de beta.3)

El lanzamiento beta.2 demostró que la evidencia estática no basta: el preflight oficial falló en el runner por una ruta mal resuelta que ningún test auxiliar cubría. Desde beta.3, el Build Preview incluye dos jobs de ensayo **fail-closed** que el `preview_report` exige en `success` para declarar el run completo:

- **`Release preflight rehearsal (non-publishing)`**: en `ubuntu-latest` con `fetch-depth: 0` del SHA candidato, comprueba que el tag `v<versión>` no existe en el remoto, crea el tag anotado temporal solo localmente, simula el evento de tag (`GITHUB_REF_NAME`, `GITHUB_REF_TYPE`, `GITHUB_SHA` exacto, `GITHUB_EVENT_PATH`, `GITHUB_OUTPUT`), ejecuta el código real de `packaging/release/validate_release_tag.py` con `RELEASE_PREFLIGHT_TOKEN: github.token` y verifica que termina con exit code 0 y emite `tag/version/channel/sha` esperados. Comprueba además que `REPO_ROOT` coincide con `git rev-parse --show-toplevel` del runner. El tag temporal se borra siempre (`if: always()`).
- **`Release gate bounded checks rehearsal`**: replica en el mismo runner y con el mismo Python 3.11 los pasos acotados del `release_gate` oficial posterior al preflight (catálogo de excepciones, aceptación ChemName, `compileall`, Ruff acotado, tests de arquitectura y contracts acotados). Los contratos comparan esos comandos byte a byte contra `release.yml` para impedir deriva.

Estos ensayos convierten la publicación en una operación ya probada de extremo a extremo: si el ensayo falla o no puede ejecutarse fielmente, el Build Preview se declara incompleto y el candidato no es apto para recibir el tag oficial.

## Guía de aceptación manual — usar la matriz existente

Probar cada paquete por separado en una VM/perfil limpio y registrar resultado, versión, SHA completo, sistema, tester y evidencia **en `manual-acceptance-0.3.0.md`**. No crear otra matriz. Los smokes de Build Preview son evidencia de compilación, no cambian casos manuales de `NOT TESTED` a `PASS`.

### P0 — ejecutar primero; cualquier fallo bloquea

1. Antes de abrir, comprobar `preview-provenance.json`, versión del candidato (actualmente `0.3.0-beta.4`), SHA y `checksums.sha256`; AppImage debe ser Type 2 extraíble. Probar Windows portable e installer, AppImage y Flatpak local por separado, sin modificar canales públicos.
2. Arranque, apertura y cierre repetidos; verificar versión/marca, iconos, menús y herramientas esenciales. Windows/Linux: temas claro y oscuro. Registrar `START-01/05/06`, `UI-07/08` y `BRANDING-001` según corresponda.
3. Dibujar átomos/enlaces/anillos, verificar valencias/conectividad y Undo/Redo; guardar copia `.cmsn`, cerrar, reabrir y comparar estructura/propiedades. Usar `DRAW-*`, `DATA-01/12`.
4. Importar/exportar SMILES y MOL/SDF; validar conectividad, fórmula, carga y **identidad** tras cada paso (un SMILES sintácticamente válido no prueba identidad). Ejecutar explícitamente `DATA-11` para enantiómeros y E/Z; la pérdida estéreo silenciosa es P0.
5. Si hay crash, pérdida/corrupción `.cmsn`, cambio de identidad/estereoquímica o error químico silencioso grave, marcar `FAIL`, conservar el archivo/evidencia y detener esa aceptación.

### P1 — funcionalidad principal

- ChemName con etanol `CCO`, benceno `c1ccccc1` y acetamida `CC(N)=O`; comprobar barra de estado/anotación al editar y Undo/Redo. Un nombre antiguo no debe persistir; salida no confiable debe ser `N/D`. Casos `CHEMNAME-RETEST-01/UPDATE-01`.
- Etanol: descriptores RDKit frente a la referencia del caso; registrar separadamente si faltan. Clean2D en simple/aromática y comparar el grafo antes/después (`RDKIT-*`, `CLEAN-*`).
- Apariencia de paneles/diálogos, onboarding, temas, plantillas, texto, anotaciones, flechas y diagramas (`UI-*`, `UI-ONBOARDING-001`). Retest manual, no el smoke Qt, determina aceptación.

### P2 — capacidades avanzadas y limitaciones

- Assistant local/offline y resolución PubChem: registrar proveedor/procedencia/identidad por separado; una respuesta SMILES válida no demuestra que sea la molécula solicitada. La IA es experimental.
- Nombres IUPAC complejos: sólo aceptar frente a estructura exacta y referencia confiable; la campaña de robustez sigue separada y los casos no soportados pueden ser `N/D`.
- Clean2D de estructuras rígidas/macrociclos y CompChem: funciones avanzadas con límites/entorno propios; no inferir cobertura universal desde casos simples.
- Updater/distribución: confirmar que los previews no alteran beta/estable ni ofrecen actualizaciones públicas (`DIST-*`). Las pruebas automatizadas no sustituyen upgrade/uninstall manual.

**La prioridad de ejecución no reescribe la severidad ni los resultados ya documentados en la matriz.** Mantener incidentes históricos `FAILED`, y todos los casos no ejecutados como `NOT TESTED`.

## Estado de esta campaña

- El run previo [#37826597134](https://github.com/PJGV333/Chemuson/actions/runs/37826597134) pasó con cuatro grupos, pero su Linux `.AppImage` era un ejecutable PyInstaller renombrado. El propietario confirmó en ese run iconos ausentes en Windows y Linux portable, en temas claro y oscuro: **P1 — FAILED, bloquea la aceptación beta**. No se acepta como prueba del nuevo Type 2 ni de iconos corregidos.
- El build remoto [#38008585493](https://github.com/PJGV333/Chemuson/actions/runs/38008585493) es `success`, SHA `f17609f97551aa94d90fac4819c450c6ee5e1005`, versión `0.3.0-beta.1`, branch `release/v0.3.0-beta.1-prep`, `complete=true` y `publication=false`. Source validation, Windows portable/installer, AppImage Type 2, Flatpak y reporte finalizaron correctamente. Los cuatro artifacts, IDs, descargas, checksums y provenance están registrados en la sección post-push del OpenSpec de preparación.
- El test Actions [#38008585666](https://github.com/PJGV333/Chemuson/actions/runs/38008585666) del mismo SHA también terminó `success`: 2.130/2.130 IDs únicos, 2.110 passed, 20 skipped, 0 failed; Windows y Flatpak smoke incluidos. Se verificó AppImage Type 2 real (magic `AI\x02`), los iconos empaquetados y ChemName/RDKit frozen smokes. P1 histórico sigue en la matriz; esos gates automatizados no sustituyen retest del propietario.
- Los artifacts anteriores siguen siendo evidencia histórica, no se mezclan con los actuales. El run nuevo ofrece paquetes únicamente para aceptación; la retención prevista es 14 días.

**Estado: READY FOR MANUAL ACCEPTANCE; NOT READY FOR BETA PUBLICATION.** P0 de identidad/estereoquímica, persistencia, arranque/cierre y los retests P1 de UI/RDKit/ChemName siguen `NOT TESTED` hasta que el propietario pruebe los paquetes exactos y registre resultados en la matriz existente. No anunciar disponibilidad pública; no se crearon tags/releases ni se modificaron canales.

## Addendum — auditoría de CI y preview (2026-10-22)

- CI #422 / Actions [#39782991976](https://github.com/PJGV333/Chemuson/actions/runs/39782991976) pasó para SHA `f0151cf09a7200acc53cb6b371105e79fd3cac1d`: 2.598 node IDs, 2.576 passed, 22 skips y cero fallos/errores; el shard más lento duró 299.398/300 s. CI #422 ejecutó los smoke de Windows/Flatpak del workflow de tests, **no** construyó los cuatro paquetes. El resultado detallado y los skips están en el OpenSpec de estabilización.
- El Build Preview más reciente inspeccionado, [#39679607716](https://github.com/PJGV333/Chemuson/actions/runs/39679607716), es completo para su propio source SHA `61f7647b683db3eca19ea06e5ec82ab2897d640a`; no coincide con el SHA estabilizado `f0151cf…`. Sus cuatro artifacts, IDs y hashes están en la [auditoría OpenSpec](../../openspec/changes/final-stabilization-chemuson-0-3-0-beta-1/validation.md). No etiquetarlos ni usarlos como distribución del SHA estabilizado.
- Entre `61f7647b…` y `f0151cf…` cambió `src/chemuson/gui/controllers/compchem3d_controller.py`; en consecuencia, el set de cuatro paquetes del candidato exacto queda **NOT BUILT**. `build-preview.yml` valida sólo ramas `release/**-prep`, mientras el SHA actual está en `stabilize/v0.3.0-beta.1`. No se creó una rama puente, no se modificó el filtro ni se intentó eludirlo.
- Los smokes automatizados en el preview anterior (69 SVG, ChemName, RDKit/3D, AppImage Type 2/extracción y Flatpak `/app`) son evidencia sólo de `61f7647b…`. Inno se compiló, pero no se instaló; el Flatpak se importó/inspeccionó en OSTree temporal, pero no se instaló/ejecutó en escritorio. La aceptación manual permanece sin cambios y los casos sin evidencia exacta siguen `NOT TESTED`.

## Candidata RC exacta — aceptación manual (2026-10-22)

Este addendum reemplaza la evaluación anterior de candidata sin paquetes exactos para el SHA de estabilización; no reasigna los artifacts de runs anteriores. Build Preview [#39789295451](https://github.com/PJGV333/Chemuson/actions/runs/39789295451) produjo con éxito los cuatro paquetes para SHA `a38fa8045c858f8ff4b6002950f57e50292aea81`, rama `release/v0.3.0-beta.1-rc-prep`, versión `0.3.0-beta.1`, `publication=false`. Los artifact IDs, tamaños, hashes de paquete/manifiestos y verificaciones están en la [validación OpenSpec](../../openspec/changes/final-stabilization-chemuson-0-3-0-beta-1/validation.md). Candidata disponible para **aceptación manual**; esto no autoriza la publicación. Artifacts de Actions con retención de 14 días: descárgalos pronto, extrae cada ZIP por separado y no mezcles artifacts de runs anteriores.

### CachyOS x86_64

1. Para el AppImage, comprueba desde el directorio extraído del artifact: `sha256sum -c checksums.sha256`; el nombre esperado termina `preview-a38fa804-linux-x86_64.AppImage`. Ejecuta `chmod +x ./Chemuson-v0.3.0-beta.1-preview-a38fa804-linux-x86_64.AppImage`, después `./Chemuson-v0.3.0-beta.1-preview-a38fa804-linux-x86_64.AppImage --version` (esperado `0.3.0-beta.1`) y finalmente lánzalo desde el escritorio. Si el montaje/FUSE no está disponible, prueba `--appimage-extract` y ejecuta `./squashfs-root/AppRun`; anota que fue extracción, no instalación.
2. Para Flatpak, verifica su artifact aparte con su propio `sha256sum -c checksums.sha256`. Requiere el runtime KDE 6.10 (`org.kde.Platform//6.10`). Instala sólo el bundle local: `flatpak install --user --bundle ./Chemuson-v0.3.0-beta.1-preview-a38fa804-linux-x86_64.flatpak`; ejecuta `flatpak run --user --branch=preview-a38fa804 io.github.PJGV333.Chemuson`. No añadas un remoto de ChemUSON ni selecciones un canal beta/stable público. Si quieres retirar la prueba tras registrar evidencia: `flatpak uninstall --user io.github.PJGV333.Chemuson//preview-a38fa804`.
3. Prueba AppImage y Flatpak por separado: arranque/versión, tema claro y oscuro, HiDPI, dibujo y selección, etanol `CCO` y benceno, undo/redo, guardado/reapertura `.cmsn`, ChemName y propiedades RDKit. Usa archivos de prueba y conserva logs/capturas con artifact ID y SHA; una prueba de la AppImage no acredita Flatpak.

### Windows x64

1. Descarga los artifacts separados **Windows portable** y **Windows Inno Setup**. Extrae cada ZIP por separado; en cada carpeta ejecuta `Get-FileHash .\<nombre-del-paquete> -Algorithm SHA256` y compara el resultado con la línea correspondiente de `checksums.sha256` (o con la fila exacta del artifact en OpenSpec). No reutilices el portable de otra ejecución.
2. Ejecuta primero el portable en un perfil limpio o VM Windows x64, sin permisos de administrador salvo que sean imprescindibles. Luego prueba el setup en una VM/perfil desechable: registrar instalación, inicio, versión, cierre y desinstalación. No actualices una instalación estable ni habilites un canal público.
3. En ambos, prueba inicio/cierre, dibujo de etanol/benceno/amida, guardar y reabrir `.cmsn`, undo/redo, ChemName, descriptores RDKit, iconos en tema claro/oscuro y DPI. Registra portable y setup como ejecuciones distintas, cada una con artifact ID, SHA, versión de Windows, resultado y evidencia.

Los checks automatizados del run (iconos, ChemName, RDKit/3D, AppImage y Flatpak) no son aceptación de instalación ni GUI del propietario. Mantén las filas manuales sin evidencia en `docs/release/manual-acceptance-0.3.0.md` como `NOT TESTED`; no marques PASS hasta completar y documentar esa ejecución exacta. No se creó tag, Release, merge, PR ni publicación.

## Intento graphics RC — Build Preview #39848401263 (2026-10-22)

El run usó exactamente `release/v0.3.0-beta.1-graphics-rc-prep`, SHA `fa51c896992c2d72d01ebf03deca316d2cce175d`, versión `0.3.0-beta.1`, `publication=false`. Windows portable/installer y Flatpak pasaron, pero AppImage falló en el smoke gráfico Xvfb/XCB. Aunque el CArchive y los smokes de iconos, RDKit y ChemName pasaron, Qt no pudo inicializar el plugin `xcb`; el build no llegó al manifest/upload. El AppImage efímero tuvo digest `47042dbc20e59fe6f9b4bb5753796b4552637b5fc7178d0f2e45d1b79053c28a`, pero no tiene artifact ID ni enlace descargable y **no debe usarse para aceptación**. El reporte del run declara `complete=false`; sigue siendo NO-GO para aceptar/publicar.

Para ese run fallido no hay archivo que abrir; no uses su digest temporal. En un artifact de preview completo, descarga el ZIP `chemuson-preview-linux-appimage`, verifica el SHA/provenance, extráelo en Dolphin, clic derecho al `.AppImage` → **Propiedades → Permisos → Es ejecutable**, cierra Propiedades y haz doble clic; si KDE pregunta, elige **Ejecutar**. Los casos manuales de la matriz siguen `NOT TESTED` hasta ejecutar y documentar el paquete exacto.

El incidente Flatpak se investigó aparte: `test` run #39809759478 fue cancelado mientras apt descargaba 123 paquetes/99 MB durante el límite existente de 10 min. En #39848401263 la instalación tardó cerca de 38 min y luego Flatpak pasó; no se atribuye a código funcional ni se aumentaron timeouts.

## Seguimiento graphics RC — Build Preview #39867831993 (2026-10-22)

El intento posterior usó la misma ref RC en el SHA `7b1c450321567fbe4903393ee530be87caee102c`, versión `0.3.0-beta.1`, `publication=false`. Source validation, Windows portable/installer y Flatpak pasaron. AppImage Type 2 y sus gates de CArchive, iconos, RDKit/3D worker y ChemName pasaron; XCB inició y creó OpenGL 4.5 Core/llvmpipe con Mesa 23.2.1. El gate falló después, cuando la escena compleja quedó vacía mientras el input contenía 24 átomos y 25 enlaces. El reporte mostró renderer/scene con 0 vértices; el run terminó antes de provenance/checksums/upload. **No hay artifact AppImage ni enlace descargable**; el run declara `complete=false` y no sirve para aceptación manual.

La causa se identificó como una carrera del harness: dibujar la fixture 2D activa un timer de sincronización 3D de 180 ms que reemplazaba la escena inyectada por el test. El smoke ya detiene ese timer y suprime el scheduling automático sólo en su ventana de prueba. El binario congelado local pasó bajo Xvfb/Mesa 23.2.1 para escenas simple y compleja, pero esto no sustituye un Build Preview correcto ni la aceptación manual. **Inicio, editor 2D, Clean2D, visor 3D y cierre en CachyOS siguen `NOT TESTED`.**

El AppImage corregido del Build Preview #39869114647 está disponible en la sección siguiente. No uses los archivos temporales ni hashes de logs de jobs fallidos de los intentos anteriores para aceptación. No se creó tag, Release, merge, PR ni publicación.

## AppImage corregido listo para aceptación manual — Build Preview #39869114647

Build Preview [#39869114647](https://github.com/PJGV333/Chemuson/actions/runs/39869114647) completó todos los gates (`complete=true`, `publication=false`) desde `release/v0.3.0-beta.1-graphics-rc-prep`, SHA `7c7a02c51c2af5d32b320839b77d544a0271af1c`, versión `0.3.0-beta.1`.

- **Descarga:** [AppImage artifact 12465108715](https://github.com/PJGV333/Chemuson/actions/runs/39869114647/artifacts/12465108715)
- **Archivo:** `Chemuson-v0.3.0-beta.1-preview-7c7a02c5-linux-x86_64.AppImage` — AppImage Type 2 (`AI\x02`), 175.350.976 bytes.
- **SHA-256:** `6f9e8205820c8662fca316fb9312590d2236ab29d6ba93a5ad6689922c26619e`. El manifest `checksums.sha256` y `preview-provenance.json` descargados fueron verificados contra el archivo; el reporte [12465103049](https://github.com/PJGV333/Chemuson/actions/runs/39869114647/artifacts/12465103049) declara `success` para source, Windows, AppImage y Flatpak.
- **Gate gráfico CI:** XCB/Xvfb + Mesa 23.2.1 llvmpipe; OpenGL 4.5 Core; canvas 2D visible; escenas simple (11.376 vértices) y compleja (31.464); ventana/dock cerraron limpiamente. El artifact descargado además se extrajo sin FUSE, `AppRun --version` devolvió `0.3.0-beta.1` y el smoke gráfico local pasó sobre Xvfb/Mesa 26.

El paquete ya está disponible para probar en CachyOS/KDE siguiendo los pasos de arriba. El smoke automatizado y la extracción local no son aceptación del propietario: inicio interactivo, editor 2D, Clean2D, visor, cierre y demás filas permanecen `NOT TESTED` hasta registrar la ejecución con este SHA.

## Salida de builds fuera del repositorio

Los scripts `packaging/linux/build_appimage.sh`, `build_flatpak.sh` y `build_linux_all.sh` honran `CHEMUSON_OUTPUT_DIR`: si está definido, `dist/`, `dist-appimage/` y `dist-flatpak/` se escriben bajo esa ruta (los argumentos posicionales específicos siguen teniendo prioridad; sin la variable, el comportamiento cambia nada). Así los binarios generados no ensucian el árbol de Git. Ejemplo: `CHEMUSON_OUTPUT_DIR=/tmp/chemuson-out bash packaging/linux/build_linux_all.sh`.

## Preview de la beta integrada en stabilize — run #39881978237 (2026-10-23)

Tras integrar `fix/beta1-appimage-graphics` en `stabilize/v0.3.0-beta.1` (fast-forward al SHA `bb53ba2f…`, más el registro de aceptación y un commit de manifiesto CI), el Build Preview de la beta integrada se ejecutó por el mecanismo oficial RC `-prep`:

- **Run:** [#39881978237](https://github.com/PJGV333/Chemuson/actions/runs/39881978237), rama `release/v0.3.0-beta.1-integration-rc2-prep`, source SHA exacto `7d935669300053ee222fc664444bc00a37e558c9`, versión `0.3.0-beta.1`, `publication=false`, `complete=true`. **Sin suite masiva**: sólo validación de fuente, gates ligeros de paquete y smoke funcional mínimo.
- **Jobs:** Validate preview source, Windows portable and installer, Linux AppImage Type 2, Linux Flatpak bundle y Preview build report: todos `success`. El log del job AppImage registra `Frozen AppImage normal 2D->3D conformer flow smoke passed`; Windows y AppImage muestran el worker RDKit congelado devolviendo el conformero de etanol (`ok=true`, mmff94, energía −1,3368570618), igual que fuente.
- **Artefactos y SHA-256 verificados localmente** (`sha256sum -c checksums.sha256` y provenance con `git_sha=7d935669…`): AppImage [12470444731](https://github.com/PJGV333/Chemuson/actions/runs/39881978237/artifacts/12470444731) `5611d2699ee16f436e39dbda2a777dd22dd04c2f6bac9d605ea8021ffd93e22e` (175.355.072 B, cabecera `AI\x02`); Windows portable [12470718764](https://github.com/PJGV333/Chemuson/actions/runs/39881978237/artifacts/12470718764) `9871633fa5ceae4f074e8a8be5abe57646c76ba6abe02513142f09e127bdd85a`; instalador [12470713964](https://github.com/PJGV333/Chemuson/actions/runs/39881978237/artifacts/12470713964) `1771cc9b5d0485051e9bb5935e381da79c5f7abe0c88ad9c493757ebd51807e0`; Flatpak [12470489974](https://github.com/PJGV333/Chemuson/actions/runs/39881978237/artifacts/12470489974) `458d6e0f151790d51ab55c72168043a8e72d5b7913cb5bae49339a82686e0482`; reporte [12470524536](https://github.com/PJGV333/Chemuson/actions/runs/39881978237/artifacts/12470524536).
- El run anterior [#39881197410](https://github.com/PJGV333/Chemuson/actions/runs/39881197410) (SHA `201a81cb…`) también pasó los cuatro jobs, pero quedó superado por un commit posterior que sólo añade IDs al manifiesto de shards de pruebas; se conserva como historial y no debe usarse para retest.

Estos paquetes son **evidencia de compilación**, no aceptación: el propietario ya validó visualmente apertura/2D/dibujo/visor 3D con cafeína en el build de diagnóstico. Las filas restantes de la matriz (`manual-acceptance-0.3.0.md`) siguen `NOT TESTED` hasta ejecutar estos paquetes exactos. No se creó tag, GitHub Release, merge con `main` ni publicación; sigue pendiente la aceptación final y la auditoría de seguridad.

## Preview corregido de Pillow/Flatpak beta.1 — run #39900808044 (2026-10-23)

Build Preview [#39900808044](https://github.com/PJGV333/Chemuson/actions/runs/39900808044) pasó source validation, Windows portable/installer, AppImage Type 2, Flatpak y reporte para el SHA exacto `e4710b938f0096ff2d02f35010117bb9a6ed8286` de `release/v0.3.0-beta.1-security-prep`; versión `0.3.0-beta.1`, `complete=true`, `publication=false`. Los cuatro paquetes son artifacts temporales de preview; no etiquetar ni publicar este build.

- Windows portable [12479081901](https://github.com/PJGV333/Chemuson/actions/runs/39900808044/artifacts/12479081901): 127.027.888 B, SHA-256 `d107a43dadb64c137dafb6e0278cc38118e2b911559c411654034ca458846754`.
- Windows installer [12478962316](https://github.com/PJGV333/Chemuson/actions/runs/39900808044/artifacts/12478962316): 128.133.031 B, SHA-256 `9859c922b2be0707a38afdb3c2d316b8d7125a0c516c6d4dc25e57d9e3c58f64`.
- AppImage Type 2 [12478253483](https://github.com/PJGV333/Chemuson/actions/runs/39900808044/artifacts/12478253483): 169.608.384 B, SHA-256 `4abb85948b3c4eccd41d51340cded2469cf7e8b5a3d04ce5eca1bbfc53009169`.
- Flatpak [12478758199](https://github.com/PJGV333/Chemuson/actions/runs/39900808044/artifacts/12478758199): 104.897.448 B, SHA-256 `828faf67f74496040594fa396530be0c1e2af92271b9eb9f6ab226e45d2556c3`.
- Build report [12478957557](https://github.com/PJGV333/Chemuson/actions/runs/39900808044/artifacts/12478957557).

Los cuatro ZIP se descargaron; checksums y provenance concuerdan en artifact, source SHA, branch, versión y `publication=false`. Escaneo independiente del CArchive real encontró **0** recursos top-level `PIL`/`Pillow` tanto en Windows portable como en AppImage; el OSTree tree extraído de Flatpak carece de `site-packages/PIL` y metadata Pillow. La resolución transitiva de Pillow 13.0.0 persiste sólo en el entorno PyInstaller por metadata upstream de RDKit; no está en los paquetes congelados. Recursos `rdkit/sping/PIL` son internos legacy de RDKit, no el paquete Pillow. OSV no mostró advisories para 13.0.0 en la consulta del 2026-10-23. No se eliminó `--share=network`: sigue siendo necesario para funciones remotas y está documentado como riesgo de transmisión de datos.

El preview #39899351520 se conserva como historial, pero **no** es candidato: sus CArchives contenían cinco entradas `PIL/*`. El Build Preview más reciente #40144590463 pertenece al SHA `85fa4ed2bd0afe1d31953440228fb73a69b54e7e`; sus paquetes son temporales y no son assets oficiales. El propietario informó el 2026-10-24 una aceptación manual limitada por formato; al no aportar run/artifact/SHA, la matriz no atribuye ese resultado a los previews ni al tag. Los casos no observados —incluido portal Flatpak, persistencia y funciones fuera de su mensaje— conservan `NOT TESTED`. La protección `v*` y los environments beta están configurados.

Artifacts del preview #40144590463, todos con provenance SHA `85fa4ed2bd0afe1d31953440228fb73a69b54e7e`, versión `0.3.0-beta.1`, branch `release/v0.3.0-beta.1-prep`, `complete=true`, `publication=false`; tamaños y hashes binarios cotejados localmente con manifests y `checksums.sha256`:
- Windows portable [12579789321](https://github.com/PJGV333/Chemuson/actions/runs/40144590463/artifacts/12579789321): 127028997 B; SHA-256 `d0fb4c8e17eb645587ba43d14bb5204ef88dcdc3478eb933f766530c287213c5`.
- Windows installer [12579674468](https://github.com/PJGV333/Chemuson/actions/runs/40144590463/artifacts/12579674468): 128135152 B; SHA-256 `02dcedbd84e2ced93f4b1dc17e0af569240f43d85897510de7028813ab8a2e81`.
- Linux AppImage Type 2 [12580420618](https://github.com/PJGV333/Chemuson/actions/runs/40144590463/artifacts/12580420618): 169604288 B; SHA-256 `e5f27656b598b7479e350a1e4d6cbd242120eb9e883146a0e11f574cebeb9d95`.
- Linux Flatpak [12580073877](https://github.com/PJGV333/Chemuson/actions/runs/40144590463/artifacts/12580073877): 104895440 B; SHA-256 `1a71cde8f18cbfa85cc4b8b4f3113c37f5f693abcb85981a567dcf8a00fd8866`.
- Report [12580505627](https://github.com/PJGV333/Chemuson/actions/runs/40144590463/artifacts/12580505627), `complete=true`.

Windows/AppImage `pillow_archive_entries=[]`; AppImage real GUI/2D→3D, icon, ChemName and RDKit checks passed. OSTree checkout of the exact Flatpak contained no `site-packages/PIL`, Pillow package or metadata. **Bloqueo de publicación:** el tag `v0.3.0-beta.1` apunta al SHA verificado, pero el commit lleva `[skip ci]`; no se disparó `release.yml` tag-only y `workflow_dispatch` falla HTTP 422. No se moverá ni recreará el tag; faltan GitHub Release, assets/sidecars oficiales y publicación Flatpak. Consultar el addendum final de `AGENT_REPORT.md`/OpenSpec; requiere decisión de recuperación del propietario.

## Preparación determinista del entorno Flatpak (beta.4, 2026-10-26)

Build Preview #40375931525 (SHA `670a7670235ec1f8833cb03c7d2475bbf5a2589f`) canceló
`Linux Flatpak bundle` a los 45 min. La API pública de Actions muestra el motivo: el
paso `Install Flatpak build tools` duró **43 min 59 s** y dejó sólo 2 min 24 s para el
empaquetado, que necesita ~4 min. En los runs previos del mismo job ese mismo paso tardó
**0,3–0,7 min** (beta.1–beta.3: 0,32/0,35/0,43/0,48/0,52/0,67 min).

No es una regresión de ChemUSON: `apt-get install -y flatpak flatpak-builder ostree`
descarga ~99 MB en ~123 paquetes. En los runners alojados por GitHub el origen apt es
`mirror+file:/etc/apt/apt-mirrors.txt` con `azure.archive.ubuntu.com` como prioridad 1;
cuando ese nodo está degradado, apt permanece en él y baja a ~60 kB/s sin devolver error
(conocido upstream en actions/runner-images [#14823](https://github.com/actions/runner-images/issues/14823),
[#14594](https://github.com/actions/runner-images/issues/14594) y
[#12949](https://github.com/actions/runner-images/issues/12949)). El mismo incidente ya
está registrado arriba para el job Flatpak de `test.yml`.

Contrato actual de preparación (jobs `build_linux_flatpak`, `build_flatpak` y
`publish_flatpak_remote`):

- `packaging/linux/install_flatpak_tools.sh` sustituye al `apt-get` suelto y
  `packaging/flatpak/toolchain.lock` fija imagen del runner, paquetes, versiones mínimas
  admitidas, mirrors oficiales permitidos y límites de tiempo por fase.
- Sondea los cuatro orígenes oficiales de Ubuntu (6 s cada uno), elige el más rápido y
  reordena el `mirror+file` del runner; si el elegido falla, rebota al siguiente. Nunca se
  usan repositorios de terceros ni se desactiva la verificación de firmas.
- Descarga acotada (`--download-only`), publica los `.deb` en
  `${{ runner.temp }}/chemuson-flatpak-apt` y **instala sin red** con `--no-download`; la
  siguiente ejecución restaura esa caché con `actions/cache` (clave derivada del hash del
  lockfile). apt vuelve a verificar cada `.deb` contra el índice firmado, así que una
  caché obsoleta o alterada se vuelve a descargar en lugar de instalarse.
- `packaging/linux/build_flatpak.sh` acota `flatpak-builder` (1500 s),
  `build-update-repo` (300 s) y `build-bundle` (300 s), y registra los tiempos por fase en
  el resumen del paso; si se agota la fase se reporta la causa en lugar de consumir el job.
- `timeout-minutes` se conserva (45 min en los jobs de empaquetado, 25 min en el de
  publicación): el control pasa a estar dentro de cada fase, no inflando el presupuesto.

El script imprime `Flatpak toolchain preparation timings (s): probe=… update=… download=…
install=… verify=… total=…` y el builder `Flatpak packaging timings (s): builder=…
repo_update=… bundle=… packaging_total=…`; con eso se comparan preparaciones reales sin
adjudicar resultados locales a Actions. Los contratos estáticos de este comportamiento
viven en `tests/test_flatpak_toolchain_preparation.py`.
