# Compilaciones de prueba en GitHub Actions

`Build Preview` compila paquetes para evaluación manual desde ramas `release/**-prep`. El resultado son artifacts temporales de una ejecución de Actions; **no** es una publicación ni una versión disponible en el updater.

## Opción A — Desde GitHub

1. Abrir el repositorio → **Actions** → **Build Preview**.
2. Si aparece **Run workflow**, elegir la rama `release/v0.3.0-beta.1-prep` (o la rama de preparación correspondiente). La versión se obtiene de `src/chemuson/_version.py`; no se escribe manualmente.
3. Ejecutar el workflow y esperar a que la validación de rama/SHA y los tres jobs de paquetes terminen.
4. Revisar el resumen final. El run sólo es completo si Windows portable/installer, Linux portable y Linux Flatpak muestran `success`.
5. Abrir la ejecución → sección **Artifacts** y descargar por separado:
   - `chemuson-preview-windows-portable`
   - `chemuson-preview-windows-installer`
   - `chemuson-preview-linux-appimage`
   - `chemuson-preview-linux-flatpak`
   - `chemuson-preview-build-report` (estados por plataforma)
6. En cada paquete, revisar `preview-provenance.json` y verificar los archivos listados en `checksums.sha256`. Confirmar versión, rama, SHA completo, sistema y `publication: false` antes de instalar.
7. Probar en VM/perfil aislado y registrar los resultados en [la matriz manual](manual-acceptance-0.3.0.md). No usar el setup/Flatpak preview sobre un perfil de trabajo sin respaldo: los instaladores comparten el App ID de la aplicación.

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
- El portable Linux conserva `.AppImage` en el nombre por compatibilidad, pero no es un AppImage Type 2; para preview omite `.updateinfo`, `.update.json` y `.zsync`.
- El Flatpak es un bundle local de rama `preview-<sha>`; no se produce un remoto público ni un `.flatpakref` actualizable.
- Los tests de contrato inspeccionan permisos, SHA, nombres y ausencia de pasos de publicación. Son evidencia **estática** únicamente; no prueban compilación real.
- La retención esperada de artifacts es 14 días. Descargarlos pronto si se necesita conservar evidencia.

## Estado de esta campaña

- `PREVIEW BUILD INFRASTRUCTURE: NOT READY` para uso remoto: workflow, builders y contratos están en la rama local y pasan validación estática, pero el código aún no está en GitHub y no se ha probado en runners reales.
- `PREVIEW ARTIFACTS: NOT BUILT`: no hay binarios de Windows portable/setup, Linux portable ni Flatpak de esta campaña.

No anunciar disponibilidad pública a partir de estos artifacts. El estado cambiará sólo después de publicar normalmente la rama autorizada, revisar el aislamiento y completar una ejecución Actions con los cuatro jobs exitosos.
