# Hotfix releases

Un hotfix publicado necesita una versión y un tag nuevos. El updater entrega Releases oficiales por canal; no distribuye commits sueltos ni previews de Actions.

## Preparación segura

1. Corregir el defecto en una rama y actualizar `_version.py` junto con AppStream mediante `packaging/release/set_version.py` antes del commit/tag revisado.
2. Construir un preview desde `release/**-prep` con [PREVIEW_BUILDS.md](release/PREVIEW_BUILDS.md) y probar los artifacts en VM/perfil aislado. El preview no es actualización y no modifica beta/stable.
3. Registrar aceptación, versión y SHA. Un hotfix beta incrementa el número prerelease (`X.Y.Z-beta.N`); un hotfix stable posterior a `X.Y.Z` incrementa PATCH (`X.Y.(Z+1)`). Nunca reutilizar un tag ni intentar reemplazar el asset de un release existente.
4. Sólo tras aprobación del propietario, crear el tag protegido `vX.Y.Z-beta.N` o `vX.Y.Z`. `.github/workflows/release.yml` valida que tag, `_version.py`, AppStream y SHA sean idénticos; el canal deriva del tag. No hay dispatch manual de versión/canal.
5. Confirmar resultados del gate y de cada plataforma, checksums, procedencia y estado de publicación. Si se interrumpe después de publicar, detener la promoción: no mover tags, no sobreescribir assets y no forzar una segunda publicación del mismo tag.

## Reglas de canal

- `beta.N` y `rc.N` son prereleases y se enrutan sólo a beta.
- `X.Y.Z` sin prerelease se enruta sólo a stable.
- Los usuarios beta pueden estar sujetos a la política del updater existente; no se debe relabelar un paquete ni cambiar manualmente un manifest para promoverlo.
- La configuración de rulesets y GitHub environments es un prerrequisito externo que debe verificar el propietario.

## Compilaciones sin publicación

Para revisar un hotfix antes de autorizar tag/release, usar la ejecución **Build Preview** y descargar sus cuatro artifacts de Actions. Si el workflow no está disponible en GitHub o el usuario/agente no está autorizado, dejar la prueba pendiente; no improvisar una ruta con permisos ampliados.

La guía completa de SemVer, recuperación e inmutabilidad está en [VERSIONING_POLICY.md](release/VERSIONING_POLICY.md).
