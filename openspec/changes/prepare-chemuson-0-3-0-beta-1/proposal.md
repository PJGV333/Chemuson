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
- Añadir compilaciones preview independientes para ramas de preparación: portable e instalador Windows, ejecutable portable Linux con sufijo histórico `.AppImage`, y bundle Flatpak. Adjuntar checksums y procedencia sin publicar nada.
- Aislar por contrato el workflow preview de Releases, tags, canales, manifests públicos y `gh-pages`, con permisos `contents: read` y pruebas estáticas.
- Auditar Windows, Flatpak y el ejecutable portable Linux existente; registrar límites que no puedan probarse en este host.

## Capabilities

### New Capabilities
- `versioned-release-pipeline`: SemVer, fuente/versiones derivadas, protección de tags y canales, gate de mismo SHA, integridad/procedencia de artefactos y rollback.
- `manual-release-acceptance`: matriz manual reproducible y criterios explícitos de promoción beta→stable.
- `preview-build-pipeline`: paquetes descargables de ramas de preparación, con procedencia/checksums, SHA idéntico en todas las plataformas y aislamiento verificable respecto a publicación oficial.

### Modified Capabilities

Ninguna. El updater conserva sus contratos actuales; la campaña endurece el proceso de preparación/publicación.

## Impact

Afecta `.github/workflows/release.yml`, el nuevo `.github/workflows/build-preview.yml`, scripts de `packaging/linux/` y `packaging/release/`, `_version.py`, AppStream/Inno metadata, tests de release/versionado/aislamiento preview, `docs/release/`, README y OpenSpec. No añade dependencias runtime, no toca química, Clean2D, `.cmsn`, `gh-pages` ni canales publicados. No crea tag ni GitHub Release; los previews sólo suben artefactos a la ejecución de Actions.
