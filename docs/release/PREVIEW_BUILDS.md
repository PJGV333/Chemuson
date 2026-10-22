# Compilaciones de prueba en GitHub Actions

`Build Preview` compila paquetes para evaluación manual desde ramas `release/**-prep`. El resultado son artifacts temporales de una ejecución de Actions; **no** es una publicación ni una versión disponible en el updater.

## Opción A — Desde GitHub

1. Abrir el repositorio → **Actions** → **Build Preview**.
2. Si aparece **Run workflow**, elegir la rama `release/v0.3.0-beta.1-prep` (o la rama de preparación correspondiente). La versión se obtiene de `src/chemuson/_version.py`; no se escribe manualmente.
3. Ejecutar el workflow y esperar a que la validación de rama/SHA y los tres jobs de paquetes terminen.
4. Windows y Linux ejecutan `validate_packaged_icons.py` y `validate_packaged_rdkit_worker.py` sobre el ejecutable PyInstaller real antes de subirlo. El gate RDKit comprueba imports nativos dentro del bundle, descriptores conocidos de etanol, SMILES de entrada/canónico y conformero 3D mediante el worker aislado; no admite skip. En Linux, el job construye AppImage Type 2 y `validate_appimage.py` repite ambos smokes sobre el ejecutable extraído del paquete. Estos gates automatizados no sustituyen la aceptación manual del propietario.
5. Revisar el resumen final. El run sólo es completo si Windows portable/installer, Linux AppImage Type 2 y Linux Flatpak muestran `success`.
6. Abrir la ejecución → sección **Artifacts** y descargar por separado:
   - `chemuson-preview-windows-portable`
   - `chemuson-preview-windows-installer`
   - `chemuson-preview-linux-appimage`
   - `chemuson-preview-linux-flatpak`
   - `chemuson-preview-build-report` (estados por plataforma)
7. En cada paquete, revisar `preview-provenance.json` y verificar los archivos listados en `checksums.sha256`. Confirmar versión, rama, SHA completo, sistema y `publication: false` antes de instalar. El AppImage debe devolver la misma versión con `./Chemuson-...AppImage --appimage-extract` en una carpeta temporal; no se requiere FUSE.
8. Probar en VM/perfil aislado y registrar los resultados en [la matriz manual](manual-acceptance-0.3.0.md). El smoke automatizado sí verifica que Qt rasteriza los iconos empaquetados, pero no sustituye inspección visual interactiva ni la aceptación del propietario. No usar el setup/Flatpak preview sobre un perfil de trabajo sin respaldo: los instaladores comparten el App ID de la aplicación.

**Disponibilidad del dispatch:** GitHub sólo permite `workflow_dispatch` cuando reconoce la definición del workflow en el repositorio. Si la opción no aparece mientras la campaña sigue en una rama de preparación y el workflow aún no está disponible desde la rama por defecto, no hacer merge ni cambiar permisos para forzar el dispatch. El push filtrado a `release/**-prep` es el disparador automático previsto una vez que el workflow se encuentre en GitHub; la rama debe publicarse mediante el procedimiento normal autorizado.

Un push a la rama también puede disparar el workflow existente `test.yml`, que escucha ramas generales y ejecuta la suite completa. Ese resultado es independiente del preview; no se relaja ni se oculta aquí. La deuda Qt/teardown histórica sigue documentada y no se repitió localmente como parte de esta campaña.

## Opción B — Desde un agente autenticado

Un agente sólo debe iniciar un run si el propietario autorizó la acción y GitHub informa que el workflow ya está disponible. Primero comprobar `gh auth status` y que la sesión tenga permiso de ejecutar workflows; no asumir que Luna tiene credenciales. Sin autorización, detenerse y dejar los paquetes como `NOT BUILT`.

```bash
gh workflow run build-preview.yml --ref release/v0.3.0-beta.1-prep
gh run list --workflow build-preview.yml --branch release/v0.3.0-beta.1-prep --limit 5
gh run view <run-id> --web
gh run download <run-id> --dir "chemuson-preview-<run-id>"
```

El agente debe devolver la URL del run, SHA/versión y resultado de cada job, además de indicar si se generaron los cuatro grupos. Un run parcial o fallido **no** se informa como build exitoso. No usar API/token con permisos adicionales ni otro workflow para eludir el filtro.

## Aislamiento y alcance real

- `GITHUB_TOKEN` tiene `contents: read`; el preview no recibe secrets de firma/publicación.
- El workflow no crea tags/releases, no hace `git push`, no despliega `gh-pages`, no construye una URL pública ni actualiza canales Flatpak beta/stable o manifests del updater.
- Linux es un AppImage Type 2 auténtico con AppRun/desktop/SVG/AppStream y recursos PyQt6/ChemUSON validados; el preview omite update-information embebido y `.updateinfo`, `.update.json`, `.zsync`.
- El Flatpak es un bundle local de rama `preview-<sha>`; no se produce un remoto público ni un `.flatpakref` actualizable.
- Los tests de contrato inspeccionan permisos, SHA, nombres y ausencia de pasos de publicación. Son evidencia **estática** únicamente; no prueban compilación real.
- La retención esperada de artifacts es 14 días. Descargarlos pronto si se necesita conservar evidencia.

## Guía de aceptación manual — usar la matriz existente

Probar cada paquete por separado en una VM/perfil limpio y registrar resultado, versión, SHA completo, sistema, tester y evidencia **en `manual-acceptance-0.3.0.md`**. No crear otra matriz. Los smokes de Build Preview son evidencia de compilación, no cambian casos manuales de `NOT TESTED` a `PASS`.

### P0 — ejecutar primero; cualquier fallo bloquea

1. Antes de abrir, comprobar `preview-provenance.json`, versión `0.3.0-beta.1`, SHA y `checksums.sha256`; AppImage debe ser Type 2 extraíble. Probar Windows portable e installer, AppImage y Flatpak local por separado, sin modificar canales públicos.
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

Cuando haya un artifact AppImage de un Build Preview completo: descarga el ZIP `chemuson-preview-linux-appimage`, extráelo en Dolphin, clic derecho al `.AppImage` → **Propiedades → Permisos → Es ejecutable**, cierra Propiedades y haz doble clic; si KDE pregunta, elige **Ejecutar**. Prueba el archivo exacto cuyo SHA y provenance constan en el artifact. Para este run no hay archivo que abrir; los casos manuales de la matriz permanecen `NOT TESTED`.

El incidente Flatpak se investigó aparte: `test` run #39809759478 fue cancelado mientras apt descargaba 123 paquetes/99 MB durante el límite existente de 10 min. En #39848401263 la instalación tardó cerca de 38 min y luego Flatpak pasó; no se atribuye a código funcional ni se aumentaron timeouts.

## Seguimiento graphics RC — Build Preview #39867831993 (2026-10-22)

El intento posterior usó la misma ref RC en el SHA `7b1c450321567fbe4903393ee530be87caee102c`, versión `0.3.0-beta.1`, `publication=false`. Source validation, Windows portable/installer y Flatpak pasaron. AppImage Type 2 y sus gates de CArchive, iconos, RDKit/3D worker y ChemName pasaron; XCB inició y creó OpenGL 4.5 Core/llvmpipe con Mesa 23.2.1. El gate falló después, cuando la escena compleja quedó vacía mientras el input contenía 24 átomos y 25 enlaces. El reporte mostró renderer/scene con 0 vértices; el run terminó antes de provenance/checksums/upload. **No hay artifact AppImage ni enlace descargable**; el run declara `complete=false` y no sirve para aceptación manual.

La causa se identificó como una carrera del harness: dibujar la fixture 2D activa un timer de sincronización 3D de 180 ms que reemplazaba la escena inyectada por el test. El smoke ya detiene ese timer y suprime el scheduling automático sólo en su ventana de prueba. El binario congelado local pasó bajo Xvfb/Mesa 23.2.1 para escenas simple y compleja, pero esto no sustituye un Build Preview correcto ni la aceptación manual. **Inicio, editor 2D, Clean2D, visor 3D y cierre en CachyOS siguen `NOT TESTED`.**

Cuando exista un artifact AppImage de un Build Preview completo, verifica su SHA/provenance y usa el flujo de Dolphin de arriba. No uses archivos temporales ni hashes de logs de jobs fallidos para aceptación. No se creó tag, Release, merge, PR ni publicación.
