# Diseño: SemVer autoritativo, canales compatibles y publicación versionada

## Fuente de verdad y selección

- GitHub Releases sigue siendo el catálogo público y `.github/workflows/release.yml` el único publicador oficial. No se introduce otro endpoint.
- El `tag_name` oficial `v<SemVer>` identifica la versión; un nombre de Release no rescata un tag inválido. Releases con tags no SemVer se ignoran de forma segura.
- `list_releases()` ordena versiones según SemVer descendente. `published_at` queda como dato de presentación/caché, no determina cuál release es más reciente.
- `latest_for_channel()` usa la clasificación SemVer: `stable` admite sólo versiones sin prerelease; `beta` mantiene el comportamiento existente de aceptar beta/RC y versiones finales. El flag GitHub `prerelease` se conserva como metadato, pero no decide la compatibilidad.
- `AutoUpdateCore` y `UpdateSettings` conservan sus contratos: preferencia persistida beta/estable, veto de downgrade mediante `highest_seen_version`, chequeo de asset por plataforma y SHA-256/firma existentes.
- El parser SemVer rechaza identificadores vacíos y números prerelease con ceros iniciales; evita colapsar metadatos corruptos a una versión estable. Build metadata mantiene la comparación de precedencia SemVer existente.

## Procedimiento de publicación

- El workflow sigue siendo disparado sólo por tags protegidos y conserva gates, builders Windows/AppImage/Flatpak, permisos, validadores, provenance y checksums.
- Las notas se resuelven mediante `docs/release/${{ steps.vars.outputs.version }}-public.md`, de modo que el mismo workflow sirve a beta y stable. La ausencia de notas falla cerrada en la acción de Release.
- No se eliminan manifests `beta.json`/`stable.json`: aunque el updater integrado consulte Releases, consumidores antiguos/desconocidos no se pueden descartar en esta campaña.
- AppImageUpdate conserva `latest` para stable y `prerelease` para beta; Flatpak conserva sus ramas/remotos. Flatpak sigue deshabilitando la búsqueda interna. Windows conserva el flavor installer/portable.
- beta.3 no cambia. GitHub no suele elegir prereleases para el distintivo Latest; se mantienen los enlaces explícitos a la beta recomendada en README y las notas actuales indican que `v0.2.5` conserva Latest. No se modifica el estado de beta.3 ni se reconstruyen archivos.

## Integración documental desde main

`origin/main` está una confirmación documental por encima del merge-base común `f0fde603371255902bf0e630ca1ce032b8f16ad3`; el commit `0f8fd6d7` sólo toca README/manuales y ya describe beta.3. La rama parte del tag/HEAD beta.3 `01652e1c9e79945a050f3733adc307b47866b7bb`. Para preparar la integración futura sin traer ni descartar código, se integrará de manera explícita el cambio documental de `main` en esta rama, se revisará cada archivo resultante y se registrarán el merge-base y los commits. No habrá merge a `main`.

## Riesgos y pruebas

- El orden por SemVer podría diferir de la cronología de publicación; es intencional y los tags oficiales son inmutables.
- Un beta configurado para beta puede avanzar de `0.3.0-beta.3` a beta superior y luego a `0.3.0`; un canal estable no recibirá prereleases aunque GitHub clasifique mal el booleano.
- Probar precedencia/invalidación SemVer, selección por canal con fechas/flags contradictorios, upgrades beta→beta y beta→final, replay/downgrade, ausencia de assets, fallo de red y checksum erróneo.
- Probar contratos YAML del workflow y que la ruta de notas depende de la versión, sin ejecutar publicación.
- No cambiar código de empaquetado de AppImage/Flatpak/Windows; contratos verifican que se conservan sus mecanismos.
