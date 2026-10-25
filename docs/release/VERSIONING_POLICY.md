# Política de versionado y publicación de ChemUSON

Esta política aplica a compilaciones oficiales y canales `beta`/`stable`. Una compilación preview de Actions no es una release: no crea tags y no actualiza el updater.

## Fuente canónica

- La única fuente de versión de la aplicación es `src/chemuson/_version.py` (`__version__`). `pyproject.toml` conserva versionado dinámico desde ese valor.
- AppStream (`packaging/flatpak/io.github.PJGV333.Chemuson.metainfo.xml`) y la versión del instalador Inno deben coincidir con la versión canónica.
- `packaging/release/set_version.py` sirve para preparar un commit de versión y sus metadatos; el CI oficial **no** modifica archivos fuente durante la compilación.
- Los builders de Windows, Linux y Flatpak deben empaquetar el mismo commit y la misma versión.

## Línea de desarrollo

- `main` es la línea principal de integración. `fix/*`, `feature/*` y `release/*-prep` son ramas temporales; no se mantienen líneas beta y estable en paralelo.
- Se permite excepcionalmente una rama correctiva temporal para una versión antigua cuando sea necesario mantenerla; después se integra explícitamente a `main`.
- Las ramas de preparación se basan en el candidato exacto. Antes de integrar ramas divergentes, verificar merge-base, commits exclusivos y diff; no forzar merges automáticos ni descartar cambios documentales ya publicados en `main`.

## SemVer con MAJOR 0

Formato: `MAJOR.MINOR.PATCH[-prerelease]`.

- Mientras `MAJOR=0`, correcciones compatibles incrementan `PATCH`; una línea de funcionalidad aditiva incrementa `MINOR`. Cambios incompatibles también incrementan `MINOR` hasta que el propietario apruebe explícitamente `1.0.0`.
- El desarrollo puede usar `X.Y.Z-dev`.
- Secuencia de prepublicación: `X.Y.Z-beta.N`, opcionalmente `X.Y.Z-rc.N`, y después `X.Y.Z` estable. `N` es positivo y monotónico: una beta rechazada se corrige con la siguiente numeración (`beta.2`), nunca se reutiliza `beta.1`.
- La campaña preparó `0.3.0-beta.1`, pero su workflow tag-only no se ejecutó (`[skip ci]`). El push de `0.3.0-beta.2` sí inició el workflow oficial, que falló cerrado en el preflight del gate por un defecto de rutas en el validador; tampoco existe GitHub Release de beta.2. Los tags anotados y sus commits se preservan; no se reutilizan ni mueven. La recuperación prepara `0.3.0-beta.3` en `docs/release/0.3.0-beta.3.md`; el tag beta.3 requiere confirmación explícita del propietario tras revisar la evidencia del ensayo del preflight ejecutado en GitHub Actions.

## Tags y canales oficiales

- El único publicador oficial es `.github/workflows/release.yml`, activado por un tag SemVer `vX.Y.Z`, `vX.Y.Z-beta.N` o `vX.Y.Z-rc.N`. No hay un workflow de publicación independiente para beta y otro para stable.
- `vX.Y.Z` se publica como release final en `stable`. `beta.N` y `rc.N` se publican como prerelease en el canal `beta`; ese metadato de GitHub acompaña al SemVer, pero no sustituye la identidad ni la comparación SemVer del actualizador.
- Las notas públicas siguen `docs/release/<version>-public.md`; el workflow único carga el archivo correspondiente a la versión validada. Un archivo ausente bloquea la publicación.
- Tag, `_version.py`, entrada AppStream vigente, SHA del evento y SHA checkout deben coincidir. Un tag malformado, borrado, ya publicado, con metadatos distintos o sin respuesta verificable de GitHub bloquea el pipeline.
- No se seleccionan versión y canal en inputs independientes; el canal se deriva del tag.
- Nunca mover, borrar, volver a crear ni sobrescribir un tag/release ya publicado. Una corrección usa una versión posterior. Stable `0.3.0` requiere aprobación de aceptación; una regresión estable requiere `0.3.1`.
- Los tags `v*` deben estar protegidos mediante ruleset de GitHub, con creación/actualización/borrado restringidos a mantenedores autorizados. Para la publicación beta.1 se verificó el ruleset activo `Protect official version tags v*` (ID `25801369`); sólo el usuario propietario `PJGV333` tiene bypass para la publicación autorizada. El workflow no prueba que el ruleset exista.

## Pipeline y permisos

1. El workflow oficial se dispara por push de tag `v*`, valida tag/SHA/version y consulta que no exista ya una Release.
2. Un gate acotado corre sobre ese SHA; los builders verifican el mismo SHA. No se consideran suficientes los resultados de otra ejecución de rama/PR.
3. El artefacto contiene SHA-256 y procedencia (versión, tag, canal, SHA fuente). HMAC y GPG son opcionales si no hay secretos y se describen como opcionales.
4. `contents: read` es el permiso por defecto. Sólo el job de Release y el job posterior que publica el remoto Flatpak reciben `contents: write`; proteger los environments de publicación con reviewers y políticas de refs permitidas. Para beta.1, `beta` y `flatpak-beta` están configurados con reviewer propietario, patrón de tag `v*-beta.*` y `can_admins_bypass=false`; el reviewer puede aprobar su propio deploy (`prevent_self_review=false`). Los canales stable no forman parte de este lanzamiento.
5. Flatpak estable y beta mantienen rutas separadas. No promover/renombrar un artefacto beta como stable.
6. El workflow de preview está separado: permisos de sólo lectura, artifacts con retención temporal en Actions, sin tag, Release, `gh-pages`, canal público, firmas de publicación ni manifiesto visible para el updater. Las compilaciones de diagnóstico también permanecen fuera de publicación.

## Actualizador y compatibilidad de canales

- El actualizador integrado considera identidad oficial únicamente un tag `v<SemVer>` válido y ordena releases por precedencia SemVer; `published_at`, el título y el booleano GitHub `prerelease` no deciden cuál release es más nueva ni si el SemVer es estable.
- Una instalación con preferencia `stable` sólo recibe versiones finales más nuevas. Una instalación que optó por `beta` puede recibir betas/RC más nuevas y graduar a una versión final posterior. Las preferencias persistidas y la protección contra downgrade se conservan.
- Los canales instalados no se eliminan en esta simplificación: AppImageUpdate sigue usando `latest` para stable y `prerelease` para beta; los remotos/ramas Flatpak `stable` y `beta` y su actualización con `flatpak update` permanecen; Windows conserva la selección de instalador o portable según el tipo de instalación. Flatpak no consulta el updater interno.
- Los manifests publicados se conservan por compatibilidad mientras no se haya demostrado que no exista consumidor anterior. Descargas siguen sujetas a asset de plataforma y SHA-256/provenance; previews no son fuente de actualización.

## Latest y beta.3

GitHub normalmente reserva el distintivo **Latest** para releases que no son prerelease. No cambiar el atributo `prerelease` de `v0.3.0-beta.3`: AppImageUpdate, rastreo beta de Windows, preferencias de clientes y remotos Flatpak son contratos separados y no se promueven cambiando sólo el metadato de GitHub. Mientras la siguiente versión final no se publique, README y enlaces oficiales deben señalar la beta recomendada y explicar que sigue siendo preliminar. No renombrar ni reconstruir tags, Release o paquetes existentes.

## Beta → estable y recuperación

- La beta sirve para instalar y evaluar artefactos reales. Cada resultado manual debe registrar versión/SHA, plataforma, pasos y evidencia en la matriz de aceptación.
- Stable exige firma del propietario, casos prioritarios completados, cero P0/P1 atribuibles al candidato y decisión documentada sobre P2/P3.
- Un fallo del gate bloquea todos los builders/publicadores. Un build parcialmente fallido no es una release completa.
- Si una publicación se interrumpe después de crear un tag/release, no sobrescribir sus archivos ni volver a usar el tag. Detener la promoción del canal, revisar la ejecución y recuperar mediante el procedimiento aprobado con una versión posterior; una reparación manual debe conservar evidencia y ser explícitamente autorizada.
- La deuda de teardown Qt/SIGSEGV y los límites de verificación de proveedores IA se mantienen abiertos hasta que haya evidencia que los cierre; no se reinterpretan como resueltos por una compilación o por tests aislados.
