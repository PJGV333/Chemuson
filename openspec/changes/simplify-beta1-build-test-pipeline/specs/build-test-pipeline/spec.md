# BUILD-TEST-PIPELINE-001 — separación de pruebas y construcción

## Purpose

Mantener la detección de errores donde corresponde: pruebas proporcionales en
desarrollo local, validación ligera y fail-closed en la construcción de
paquetes, y la suite completa como auditoría extraordinaria invocable.

## ADDED Requirements

### Requirement: Full pytest suite is opt-in, never automatic per push
El workflow `test.yml` MUST dispararse únicamente mediante `workflow_dispatch`
u otra invocación manual explícita, y MUST conservar sin cambios los jobs de
plan, shards acotados, summary fail-closed y smokes, incluida la verificación
de la colección completa contra el manifiesto de node IDs. NO debe ejecutarse
automáticamente en cada push ordinario.

#### Scenario: Ordinary push to a working branch
- **WHEN** se hace push a una rama de desarrollo o estabilización
- **THEN** no se programa la suite completa de pytest en GitHub Actions
- **AND** Build Preview solo corre en refs `release/**-prep` o dispatch manual.

#### Scenario: Exceptional full audit
- **WHEN** un responsable dispara manualmente el workflow *test* sobre una referencia
- **THEN** se ejecuta la campaña completa con el mismo plan externo, límites de tiempo y verificación fail-closed que antes
- **AND** ningún job de construcción depende de esa ejecución.

### Requirement: Package builds validate without bundling test infrastructure
La construcción de paquetes MUST incluir solo código de producción,
dependencias y recursos necesarios para ejecutar ChemUSON, y MUST ejecutar
después del build las validaciones ligeras independientes: identidad de
código, integridad SHA-256, procedencia, dependencias nativas y smoke
funcional mínimo, incluido el smoke real de generación 2D→3D del ejecutable.
Estas validaciones MUST NOT sustituirse por la suite completa ni debilitarse.

#### Scenario: Build Preview run
- **WHEN** Build Preview construye Windows portable/instalador, AppImage y Flatpak
- **THEN** no ejecuta la suite completa de pytest
- **AND** ejecuta validadores de iconos, worker RDKit, ChemName y (en AppImage) `geometry_flow_smoke` sobre el binario generado
- **AND** produce manifiestos de procedencia y checksums SHA-256 sin publicar Release.

### Requirement: Build outputs live outside the source tree by configuration
Las herramientas de build locales MUST permitir redirigir todas las salidas a
un directorio externo mediante `CHEMUSON_OUTPUT_DIR` sin cambiar el
comportamiento por defecto, MUST mantener los artefactos fuera de Git y MUST
conservar los manifiestos de procedencia confirmados como evidencia histórica.

#### Scenario: External output directory
- **WHEN** se ejecuta `build_linux_all.sh` con `CHEMUSON_OUTPUT_DIR=/salidas`
- **THEN** `dist/`, `dist-appimage/`, `dist-flatpak/` y los workpath de PyInstaller cuelgan de `/salidas`
- **AND** `git status --short` sigue limpio y borrar `/salidas` no afecta a producción.
