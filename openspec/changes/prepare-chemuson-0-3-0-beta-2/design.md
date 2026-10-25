# Diseño: recuperación beta.2 con tag beta.1 preservado

## Base y aislamiento

La rama de origen es `release/v0.3.0-beta.1-prep`, verificada contra `origin` en `aa2e91c98cf1a16e8299e881955b1614edb2e18d`. La nueva rama de preparación se crea desde ese commit. El tag anotado beta.1 y su commit objetivo `85fa4ed2bd0afe1d31953440228fb73a69b54e7e` se conservan intactos. Las correcciones de seguridad y empaquetado ya integradas permanecen en la historia.

## Cambios permitidos

1. `packaging/release/set_version.py` sincroniza `src/chemuson/_version.py` y la entrada superior de AppStream para `0.3.0-beta.2`; `pyproject.toml` obtiene la versión dinámicamente y no necesita duplicarla.
2. Se agrega una nota pública versionada para beta.2 y se apunta `release.yml` a ella. La nota distingue la evidencia histórica de beta.1 de una aceptación específica de paquetes beta.2.
3. Se documenta que beta.1 quedó sin publicación y que la recuperación utiliza la numeración monotónica beta.2.
4. Se actualizan sólo las expectativas literales de los contratos existentes para la versión de ventana y la ruta de notas beta.2.
5. Se agrega el OpenSpec con baseline y verificaciones de esta preparación.

No se modifica la lógica de producción ni el contrato del pipeline. `release.yml` permanece activado exclusivamente por push de tags `v*`, con validación fail-closed del tag/SHA, tests acotados, builds desde el SHA validado, permisos mínimos, aprobaciones de environments, publicación de assets y validación del remoto Flatpak. Un tag beta.2 sólo puede considerarse después de confirmación del propietario.

## Verificación

- Baseline focal de compilación, colección de tests y Ruff; la suite completa de pytest se omite por instrucción.
- Tests acotados de versionado, semver, metadatos, validación de tags, workflows y artifacts; validaciones OpenSpec y de distribución.
- Revisar configuración activa de ruleset y environments, además de la configuración de los workflows.
- Empujar la rama candidata para ejecutar un Build Preview nuevo. Los cuatro artifacts deben reportar el mismo SHA de origen, provenance verificable y hashes descargados comprobados localmente. Los artifacts permanecen previews temporales y no son assets oficiales.

## Consecuencia

La rama puede avanzar más allá de beta.1 sin mutar su tag. El candidato listo se comunica mediante su SHA completo; no se crea tag en esta tarea, por lo que la publicación oficial y el cambio de numeración aún requieren confirmación del propietario.
