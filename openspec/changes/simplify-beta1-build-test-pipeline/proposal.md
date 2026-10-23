# Propuesta — simplificar pruebas automatizadas y organización de builds

## Problema

`test.yml` ejecuta la suite completa (2.607 tests en 8 shards más smokes) en
cada push a cualquier rama, consumiendo recursos y generando fallos ajenos al
código de producción. Además, los artefactos de build local se generaban dentro
del árbol del repositorio sin una política explícita de salida.

## Objetivo

Separar tres actividades sin perder capacidad de detectar errores:

- **A. Desarrollo local:** pruebas proporcionales al cambio; los agentes siguen
  usando la colección completa de `tests/` como herramienta; scripts
  exploratorios temporales fuera de producción.
- **B. Verificación previa al push:** pruebas focalizadas + arquitectura +
  compileall + Ruff + OpenSpec + `git diff --check` + smoke GUI cuando aplique.
  La suite completa no es obligatoria por modificación.
- **C. GitHub Actions:** solo produce paquetes validados; la suite completa se
  ejecuta de forma excepcional vía `workflow_dispatch`.

## Alcance

- `test.yml`: eliminar triggers `push`/`pull_request`; conservar `workflow_dispatch`
  y todos los jobs/herramientas de shards exactamente como están.
- Documentar la política en `docs/testing/AGENT_TESTING_POLICY.md`.
- Soporte `CHEMUSON_OUTPUT_DIR` en los scripts Linux de build para sacar los
  artefactos del repositorio; `.gitignore` ampliado para salidas de build.
- Sin cambios en química, GUI, renderer, persistencia, políticas de
  actualización/seguridad ni en los controles de integridad/procedencia de
  `Build Preview`/`release`. No se modifica `main`.

## Fuera de alcance

- Eliminar o reescribir pruebas existentes.
- Cambiar protecciones de rama, tags, Releases o publicación.
- Modificar el contenido empaquetado del binario (ya excluye `tests/`).
