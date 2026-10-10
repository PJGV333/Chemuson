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
- El build remoto [#38006731375](https://github.com/PJGV333/Chemuson/actions/runs/38006731375) es `success`, SHA `57d481d21fee2db6ecacf6a862839caa63a1c56d`, versión `0.3.0-beta.1`, branch `release/v0.3.0-beta.1-prep`, `complete=true` y `publication=false`. Source validation, Windows portable/installer, AppImage Type 2, Flatpak y reporte finalizaron correctamente. Los cuatro artifacts y descargas/checksums/provenance están registrados en la sección post-push del OpenSpec de preparación.
- El test Actions [#38006731537](https://github.com/PJGV333/Chemuson/actions/runs/38006731537) del mismo SHA también terminó `success`: 2.130/2.130 IDs únicos, 2.110 passed, 20 skipped, 0 failed; Windows y Flatpak smoke incluidos. Se verificó AppImage Type 2 real (magic `AI\x02`), los iconos empaquetados y ChemName/RDKit frozen smokes. P1 histórico sigue en la matriz; esos gates automatizados no sustituyen retest del propietario.
- Los artifacts anteriores siguen siendo evidencia histórica, no se mezclan con los actuales. El run nuevo ofrece paquetes únicamente para aceptación; la retención prevista es 14 días.

**Estado: READY FOR MANUAL ACCEPTANCE; NOT READY FOR BETA PUBLICATION.** P0 de identidad/estereoquímica, persistencia, arranque/cierre y los retests P1 de UI/RDKit/ChemName siguen `NOT TESTED` hasta que el propietario pruebe los paquetes exactos y registre resultados en la matriz existente. No anunciar disponibilidad pública; no se crearon tags/releases ni se modificaron canales.
