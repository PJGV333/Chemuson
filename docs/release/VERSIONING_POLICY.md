# Política de versionado y publicación de ChemUSON

Esta política aplica a compilaciones oficiales y canales `beta`/`stable`. Una compilación preview de Actions no es una release: no crea tags y no actualiza el updater.

## Fuente canónica

- La única fuente de versión de la aplicación es `src/chemuson/_version.py` (`__version__`). `pyproject.toml` conserva versionado dinámico desde ese valor.
- AppStream (`packaging/flatpak/io.github.PJGV333.Chemuson.metainfo.xml`) y la versión del instalador Inno deben coincidir con la versión canónica.
- `packaging/release/set_version.py` sirve para preparar un commit de versión y sus metadatos; el CI oficial **no** modifica archivos fuente durante la compilación.
- Los builders de Windows, Linux y Flatpak deben empaquetar el mismo commit y la misma versión.

## SemVer con MAJOR 0

Formato: `MAJOR.MINOR.PATCH[-prerelease]`.

- Mientras `MAJOR=0`, correcciones compatibles incrementan `PATCH`; una línea de funcionalidad aditiva incrementa `MINOR`. Cambios incompatibles también incrementan `MINOR` hasta que el propietario apruebe explícitamente `1.0.0`.
- El desarrollo puede usar `X.Y.Z-dev`.
- Secuencia de prepublicación: `X.Y.Z-beta.N`, opcionalmente `X.Y.Z-rc.N`, y después `X.Y.Z` estable. `N` es positivo y monotónico: una beta rechazada se corrige con la siguiente numeración (`beta.2`), nunca se reutiliza `beta.1`.
- La campaña actual prepara el valor canónico `0.3.0-beta.1`; a la fecha de este documento no ha creado `v0.3.0-beta.1` ni una GitHub Release.

## Tags y canales oficiales

- El único disparador oficial es un tag SemVer `vX.Y.Z`, `vX.Y.Z-beta.N` o `vX.Y.Z-rc.N`.
- `vX.Y.Z` se publica únicamente en `stable`. `beta.N` y `rc.N` se publican únicamente como prerelease en el canal `beta`.
- Tag, `_version.py`, entrada AppStream vigente, SHA del evento y SHA checkout deben coincidir. Un tag malformado, borrado, ya publicado, con metadatos distintos o sin respuesta verificable de GitHub bloquea el pipeline.
- No se seleccionan versión y canal en inputs independientes; el canal se deriva del tag.
- Nunca mover, borrar, volver a crear ni sobrescribir un tag/release ya publicado. Una corrección usa una versión posterior. Stable `0.3.0` requiere aprobación de aceptación; una regresión estable requiere `0.3.1`.
- Los tags `v*` deben estar protegidos mediante ruleset de GitHub, con creación/actualización/borrado restringidos a mantenedores autorizados. La configuración remota debe ser verificada por el propietario; el workflow no prueba que el ruleset exista.

## Pipeline y permisos

1. El workflow oficial se dispara por push de tag `v*`, valida tag/SHA/version y consulta que no exista ya una Release.
2. Un gate acotado corre sobre ese SHA; los builders verifican el mismo SHA. No se consideran suficientes los resultados de otra ejecución de rama/PR.
3. El artefacto contiene SHA-256 y procedencia (versión, tag, canal, SHA fuente). HMAC y GPG son opcionales si no hay secretos y se describen como opcionales.
4. `contents: read` es el permiso por defecto. Sólo el job de Release y el job posterior que publica el remoto Flatpak reciben `contents: write`; se deben proteger los environments `beta`, `stable`, `flatpak-beta` y `flatpak-stable` y limitar secretos a los jobs que los requieren.
5. Flatpak estable y beta mantienen rutas separadas. No promover/renombrar un artefacto beta como stable.
6. El workflow de preview está separado: permisos de sólo lectura, artifacts con retención temporal en Actions, sin tag, Release, `gh-pages`, canal público, firmas de publicación ni manifiesto visible para el updater.

## Beta → estable y recuperación

- La beta sirve para instalar y evaluar artefactos reales. Cada resultado manual debe registrar versión/SHA, plataforma, pasos y evidencia en la matriz de aceptación.
- Stable exige firma del propietario, casos prioritarios completados, cero P0/P1 atribuibles al candidato y decisión documentada sobre P2/P3.
- Un fallo del gate bloquea todos los builders/publicadores. Un build parcialmente fallido no es una release completa.
- Si una publicación se interrumpe después de crear un tag/release, no sobrescribir sus archivos ni volver a usar el tag. Detener la promoción del canal, revisar la ejecución y recuperar mediante el procedimiento aprobado con una versión posterior; una reparación manual debe conservar evidencia y ser explícitamente autorizada.
- La deuda de teardown Qt/SIGSEGV y los límites de verificación de proveedores IA se mantienen abiertos hasta que haya evidencia que los cierre; no se reinterpretan como resueltos por una compilación o por tests aislados.
