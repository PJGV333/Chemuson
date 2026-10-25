# Tareas

- [x] Confirmar estado limpio, identidad de los tags beta.1/beta.2 (anotados, SHAs objetivo `85fa4ed2…` y `6aecba5a…`) y que beta.1 es ancestro de beta.2; verificar el fallo del run oficial beta.2 (`40160549763`, exit 128 en `git rev-parse HEAD`) como evidencia de la causa raíz.
- [x] Crear rama aislada `release/v0.3.0-beta.3-prep` desde `6aecba5a5585d9560b693186578b76bea9e623b8` y capturar el baseline antes de tocar producción.
- [x] Corregir `REPO_ROOT = RELEASE_DIR.parents[1]` en `validate_release_tag.py` y auditar `parents[` en scripts relacionados (sin defectos equivalentes).
- [x] Añadir regresiones de `REPO_ROOT` vs `git rev-parse --show-toplevel`, ejecución real de `verify_tag_checkout`, `main()` simulado offline y dos subprocess fail-closed; verificar por mutación que reintroducir `parents[2]` vuelve a fallar con exit 128.
- [x] Añadir jobs `preflight_rehearsal` y `release_gate_rehearsal` al Build Preview, integrarlos al `preview_report` (fail-closed) y actualizar/añadir los contracts estáticos correspondientes.
- [x] Cambiar versión canónica/AppStream a `0.3.0-beta.3` con `set_version.py`; crear notas `0.3.0-beta.3-public.md` y documento interno; apuntar `body_path` del workflow a las notas beta.3; actualizar contratos literales de versión.
- [x] Actualizar `VERSIONING_POLICY.md` y `PREVIEW_BUILDS.md` (recuperación beta.3 y ensayos de preflight).
- [x] Ejecutar validaciones focales locales (contracts de workflows, tests de tag, versionado, artefacts, baseline exceptions, ChemName, compileall, Ruff acotado) sin suite masiva.
- [x] Publicar la rama de preparación y ejecutar el Build Preview oficial desde el SHA candidato exacto.
- [x] Verificar: ensayo del preflight aprobado en Actions (exit 0, salidas correctas, tag temporal borrado), ensayo del gate aprobado, cuatro paquetes construidos, checksums/provenance coincidentes, un único SHA fuente, sin errores de versionado ni rutas, repositorio limpio (run 40162819727 sobre `b16ccbcbc…`; evidencia en `validation.md`).
- [x] Descargar los artifacts del Build Preview y validar SHA-256/provenance localmente contra los manifests.
- [ ] Confirmar el Build Preview verde disparado por el commit final de documentación (SHA definitivo del candidato), con la misma evidencia de ensayos, cuatro paquetes y SHA único.
- [ ] Informar al propietario el SHA completo, el enlace del run y la evidencia del ensayo; no crear `v0.3.0-beta.3` sin su confirmación explícita.
