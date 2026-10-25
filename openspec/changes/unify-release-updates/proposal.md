# Propuesta: unificar selección de releases y conservar canales instalados

## Motivación

La publicación oficial ya usa un único workflow tag-only, pero la selección del actualizador mezcla el orden por fecha de GitHub con el booleano `prerelease`, aunque el proyecto publica versiones beta y finales en el mismo catálogo. Además, el workflow fija las notas de beta.3, lo que no sirve para la siguiente versión SemVer sin editar el pipeline. Esa distribución puede esconder una versión SemVer superior o hacer que la clasificación de GitHub contradiga la versión.

## Alcance

- Usar el tag SemVer como identidad y fuente de orden de cada release oficial.
- Seleccionar candidatos compatibles por comparación SemVer, no por `published_at` ni únicamente por `prerelease`.
- Mantener preferencias `stable`/`beta`, rutas AppImageUpdate y remotos/ramas Flatpak existentes, así como checksums, verificación y selección del paquete Windows.
- Hacer versionadas las notas del workflow oficial (`docs/release/<version>-public.md`) sin crear otro procedimiento de publicación.
- Alinear la política y los contratos con la línea principal única y el flujo de integración beta.3 → `main`.
- Incorporar a esta rama de trabajo el commit documental ya publicado en `main`, después de verificar su ancestralidad y que el cambio sólo incluye README/manuales; no integrar esta rama en `main`.

## Fuera de alcance

- Cambiar o reconstruir la Release, el tag, la clasificación `prerelease` o los paquetes públicos de `v0.3.0-beta.3`.
- Eliminar los canales beta/estable, manifests existentes, remotos Flatpak, metadatos AppImageUpdate o preferencias persistidas.
- Cambiar lógica química, GUI, RDKit, ChemName, datos de usuario o privilegios.
- Crear otro servicio, catálogo, workflow publicador, dependencia o automatización de merge.
- Ejecutar la suite general de pytest.

## Módulos probables

`src/chemuson/update/semver.py`, `src/chemuson/update/provider.py`, `tests/test_update_*.py`, contratos de `release.yml`, `.github/workflows/release.yml`, `docs/release/VERSIONING_POLICY.md` y el OpenSpec de esta campaña. Se preservarán literalmente las páginas de beta.3 publicadas en `origin/main`.
