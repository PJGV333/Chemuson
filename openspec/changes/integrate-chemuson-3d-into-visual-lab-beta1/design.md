# Diseño — merge aislado 3D + visual lab Beta 1

## Baseline y referencias

- Worktree aislado: `/home/ccachyavgp/Documentos/Chemuson-3D-beta1-integration`.
- Rama local temporal: `integration/chemuson-3d-beta1`, creada exactamente desde `origin/feat/molecular-depiction-visual-lab` (`43a8a27c45cdabfd157bf3f2a52c9216b24a1552`); checkpoint `checkpoint/chemuson3d-beta1-base-43a8a27` en el mismo SHA.
- Fuente 3D: `origin/feature/chemuson-3d` (`ec486626333ccd2e0ae86671861708b30c508c5b`). El árbol de fuente contiene los 36 commits propios del trabajo 3D; se conservarán como ancestros mediante merge.
- Referencia de release observada: `origin/release/v0.3.0-beta.1-prep` (`9a70ee1c908a39e700098bcc0a54d464a606e24c`). No se integrará en ella: sus especificaciones de release/QA conservan gates de pruebas y aceptación del propietario pendientes.
- Baseline detallado y limitaciones: `baseline.md`.

## Estrategia de integración

1. Mantener la rama temporal sobre la base visual exacta; preservar el checkpoint y registrar un commit local de esta especificación antes de modificar código.
2. Usar `git merge --no-ff --no-commit origin/feature/chemuson-3d`; no cherry-pick, squash ni reemplazo de historial.
3. Resolver manualmente cada conflicto. Un `git merge-tree --write-tree --name-only` de solo lectura anticipa seis: `AGENT_REPORT.md`, `architecture/modules.yml`, `src/chemuson/chemio/_rdkit_worker.py`, `src/chemuson/chemio/rdkit_io.py`, `src/chemuson/chemio/rdkit_safe.py` y `src/chemuson/gui/main_window_ui_builder.py`.
4. En `AGENT_REPORT.md`, conservar advertencias y decisiones pertinentes de ambas ramas. En el catálogo, integrar las entradas de módulos/dependencias sin borrar excepciones válidas. En los tres límites RDKit, revisar ambos lados completos y conservar compatibilidad de entradas, protocolos, cancelación/errores y usos 2D/3D; no resolver con “ours/theirs” a ciegas. En el builder de UI, preservar la calibración visual y los controles actuales de Beta 1 además del viewport/dock 3D.
5. Revisar manualmente los merges automáticos que comparten responsabilidades: `canvas_structure.py`, `canvas_view.py`, `docks.py`, `main_window.py`, `shell/assembly.py` y `tests/test_compchem3d_dock.py`. La rama visual es autoridad para sus capacidades 2D; la rama 3D es autoridad para las nuevas capacidades 3D, sin descartar comportamiento del otro lado.
6. Ejecutar primero regresiones de los conflictos y capacidades 2D/3D; luego arquitectura, persistencia/ChemIO y pruebas integradas disponibles. Validar OpenSpec estricto, compilación, colección, Ruff scoped y `git diff --check`. Para la suite extensa, usar el workflow particionado existente y registrar los timeouts/fallos; no declararla verde por cobertura parcial.
7. Sólo crear el merge commit local si el merge compila, conserva ambos historiales, supera las regresiones críticas y el árbol queda verificable. No tocar la rama de release. La aceptación manual de GUI/OpenGL, revisión visual del propietario y otros gates aún pendientes impiden afirmar aceptación completa de Beta 1 aunque los tests locales pasen.

## Riesgos y criterios de parada

- Los conflictos de `chemio` pueden afectar serialización/importación y el worker 2D aunque las pruebas 3D focales pasen; exigir pruebas de ida/vuelta `.cmsn`, importación/exportación y aislamiento del worker.
- Un merge automático en mixins, canvas o ensamblado de shell puede alterar orden de eventos, galería visual, undo o ciclo de vida del dock; revisar el diff contra la base y ejecutar pruebas de ambas familias.
- La rama objetivo documenta fallos/timeout de suite en un orden de Qt y revisiones manuales pendientes. No ocultarlos ni convertirlos en baseline nuevo.
- Si falta una dependencia/fixture, un gate crítico falla sin causa conocida, o la resolución exige cambiar contratos fuera del OpenSpec, dejar el merge detenido sin commit y reportar el bloqueo.