# Tareas — simplificar pruebas y organización de builds

## 1. Workflows
- [x] 1.1 Verificar protección de ramas y checks obligatorios antes de tocar triggers (API: ninguna rama protegida).
- [x] 1.2 `test.yml`: quitar `push`/`pull_request`, conservar `workflow_dispatch` y todos los jobs sin cambios.
- [x] 1.3 Confirmar que `build-preview.yml`/`release.yml` no dependen de `test.yml` y no empaquetan `tests/`.

## 2. Política de pruebas
- [x] 2.1 Crear `docs/testing/AGENT_TESTING_POLICY.md` con las tres actividades, comandos locales y auditoría extraordinaria.
- [x] 2.2 Verificar que ningún contrato de tests existente se debilita (suite de contracts packaging/CI verde).

## 3. Salidas de build
- [x] 3.1 `CHEMUSON_OUTPUT_DIR` en `build_appimage.sh`, `build_flatpak.sh`, `build_linux_all.sh` (+ PyInstaller distpath/workpath).
- [x] 3.2 `.gitignore` ampliado; conservar manifiestos históricos trackeados.
- [x] 3.3 `tests/test_linux_distribution_scripts.py` verde con los scripts modificados.

## 4. Integración
- [x] 4.1 ff-merge a `stabilize/v0.3.0-beta.1` y push.
- [x] 4.2 Ref RC `-prep` nueva desde el SHA integrado; Build Preview sin suite completa; verificar los cuatro paquetes y sus provenance/SHA-256.
- [x] 4.3 Documentar resultados en AGENT_REPORT y docs de release; worktree limpio.
