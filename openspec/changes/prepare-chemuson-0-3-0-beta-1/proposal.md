# Proposal — preparar ChemUSON 0.3.0-beta.1

## Why

ChemUSON tiene ya una línea 0.3.0, builds Windows/Linux/Flatpak, canales beta/stable y un updater, pero el flujo actual permite seleccionar versión y canal independientemente, publica sin gate de tests del mismo SHA y altera la versión fuente durante el build. Hace falta preparar una beta reproducible y documentar una ruta de promoción estable sin crear una release ni alterar canales en esta campaña.

## What Changes

- Establecer SemVer de aplicación y `src/chemuson/_version.py` como fuente canónica; preparar el código para `0.3.0-beta.1` sin tag.
- Sincronizar metadatos derivados y verificar en CI que tag, versión de aplicación/paquete, instalador, AppStream y manifests coinciden.
- Hacer que la publicación se origine únicamente en tags protegidos, derive el canal desde el prerelease y ejecute gates fail-closed sobre el mismo SHA antes de empaquetar/publicar.
- Retirar el dispatch que acepta versión/canal separados; evitar colisión/reutilización de releases, controlar permisos y preservar separación beta/stable.
- Mantener checksums obligatorios, declarar explícitamente el alcance opcional de HMAC/GPG y registrar procedencia verificable de los artefactos.
- Documentar excepciones baseline por test exacto, sin skips generales, y preparar notas beta, política de versionado y matriz manual beta/estable.
- Añadir compilaciones preview independientes para ramas de preparación: portable e instalador Windows, AppImage Linux Type 2 auténtico, y bundle Flatpak. Adjuntar checksums y procedencia sin publicar nada.
- Construir AppImage desde un AppDir verificable con AppRun, desktop entry, icono y bundle PyInstaller; preservar y validar el contrato existente de updater, incluido el update-information embebido.
- Corregir CI Python para instalar `requirements-dev.txt` sin duplicar dependencias runtime y proteger el manifiesto de test mediante contrato estático.
- Aislar por contrato el workflow preview de Releases, tags, canales, manifests públicos y `gh-pages`, con permisos `contents: read` y pruebas estáticas.
- Auditar Windows, Flatpak y el ejecutable portable Linux existente; registrar límites que no puedan probarse en este host.
- Bloquear la aceptación beta por los iconos SVG ausentes en los paquetes Windows/Linux; verificar recursos y rasterizado QtSvg en ejecutables congelados de preview y release.
- Registrar y corregir `UI-ONBOARDING-001`: posponer el recorrido automático hasta que el layout sea visible, conservar sus tres pasos y recalcular geometría ante resize/DPI.
- Registrar y corregir `BRANDING-001`: normalizar las superficies visibles a **ChemUSON**, sin cambiar identidades técnicas, persistencia, nombres de artefactos ni actualización.
- Mantener la aceptación manual de ambos asuntos pendiente hasta el retest del propietario en los nuevos paquetes.
- Registrar `P1 — RDKit isolated backend unavailable in packaged executable`: el propietario observó que Windows portable muestra fórmula/masa/espectros estimados pero no descriptores RDKit. Verificar por separado las extensiones nativas empaquetadas y el worker aislado real en Windows portable, Linux PyInstaller y AppImage.
- Añadir un modo worker interno al ejecutable congelado (sin GUI, sin Python externo, sin import RDKit en el padre, compatible con Windows `console=False`): el AppImage PyInstaller previo reproduce que el ejecutable recibe `_rdkit_worker.py` como argumento CLI no reconocido y termina con código 2. La disponibilidad/import de RDKit se valida por separado, sin atribuir el fallo a ausencia del paquete. Mantener aislamiento, contrato JSON, timeouts y errores controlados.
- Hacer que Preview y el workflow oficial fallen si el RDKit smoke congelado falla o no produce descriptores conocidos/SMILES/3D; no permitir `skip` porque RDKit es dependencia obligatoria.
- Mantener `P1` bloqueante y el retest manual del propietario pendiente; el SIGSEGV de teardown Qt sigue siendo deuda independiente.
- Diagnosticar `P1 — ChemName templates omitted from frozen packages` inspeccionando primero artefactos Windows/AppImage reales y Flatpak por separado; incluir explícitamente los recursos requeridos y comparar Python con nombres del ejecutable congelado sin cambiar reglas de nomenclatura.
- Completar BRANDING-001 con el texto aprobado de Acerca de y el título principal simplificado, preservando atribuciones legales e identidades internas.

## Capabilities

### New Capabilities
- `versioned-release-pipeline`: SemVer, fuente/versiones derivadas, protección de tags y canales, gate de mismo SHA, integridad/procedencia de artefactos y rollback.
- `manual-release-acceptance`: matriz manual reproducible y criterios explícitos de promoción beta→stable.
- `preview-build-pipeline`: paquetes descargables de ramas de preparación, con procedencia/checksums, SHA idéntico en todas las plataformas y aislamiento verificable respecto a publicación oficial.

### Modified Capabilities

- `ui-onboarding`: geometría fiable del recorrido actual de tres pasos, sin alterar el layout ni la semántica QSettings.
- `visible-branding`: nombre presentado como ChemUSON con preservación de identidad técnica.
- `packaged-rdkit-worker`: carga de RDKit nativo y ejecución del worker aislado desde los binarios PyInstaller distribuidos.

El updater conserva sus contratos actuales; estos cambios no cambian canales, rutas ni formatos. La corrección RDKit se limita al protocolo de lanzamiento/diagnóstico del worker, a su verificación de empaquetado y al mensaje de fallo del panel; no altera algoritmos químicos, Clean2D ni persistencia `.cmsn`.

## Impact

Afecta `.github/workflows/release.yml`, `.github/workflows/build-preview.yml`, `.github/workflows/test.yml`, scripts y AppDir de `packaging/linux/`, validadores de `packaging/release/`, `_version.py`, AppStream/Inno metadata, tests de release/versionado/preview/AppImage/CI, `docs/release/`, documentación Linux, README y OpenSpec. No añade dependencias runtime, no toca química, Clean2D, `.cmsn`, `gh-pages` ni canales publicados. No crea tag ni GitHub Release; los previews sólo suben artefactos a la ejecución de Actions.
