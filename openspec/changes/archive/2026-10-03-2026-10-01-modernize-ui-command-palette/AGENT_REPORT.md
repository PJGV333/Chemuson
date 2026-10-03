# AGENT_REPORT — Fase 6 (paleta de comandos Ctrl+P)

Cambio OpenSpec: `2026-10-01-modernize-ui-command-palette` (strict valid).
Rama `ui/modernization`.

## Estado de publicaciones

- Fase 6 inicial (`Add command palette`, commit `3c8d992`): **ya fue
  publicada por el usuario**. Antes de esta corrección,
  `origin/ui/modernization` estaba en `3c8d99208b3fb77f7019be3eba639bfcc857622e`.
- Esta corrección de UX (`Restore Clean2D shortcut and finalize command
  palette`) se apila sobre `3c8d992` y se publica con push normal (fast-forward,
  sin `--force`).

## Corrección de UX (aprobada tras prueba manual en KDE/Wayland)

El usuario aprobó funcionalidad y apariencia de la paleta con **una** corrección
de atajos:

- `action_clean_2d_full` **recupera** `Ctrl+K` (contexto `WindowShortcut` +
  `window.addAction(...)`), como antes de la Fase 6. `Ctrl+K` vuelve a ejecutar
  *Limpiar 2D (1 paso)*. No se modifica su handler ni ninguna lógica de
  Clean2D.
- `action_command_palette` pasa de `Ctrl+K` a **`Ctrl+P`** (mantiene
  `WindowShortcut`). La `SearchPill` sigue abriendo exactamente la misma
  `action_command_palette`; su badge visible cambia de `Ctrl K` a `Ctrl P`.
- `Ctrl+Shift+K` y `Ctrl+Alt+K` permanecen intactos.

**Verificación previa de conflicto (gate):** antes de asignar `Ctrl+P` se
inspeccionó programáticamente el conjunto completo de `QAction` y `QShortcut`
de la ventana; **ninguna** poseía `Ctrl+P`. No hubo conflicto; no se eligió otro
atajo.

**Alcance:** exclusivamente documental/atajo. No se modifica layout, ranking,
registro, `QAction` existentes ni apariencia de la paleta. No se toca
`src/chemuson/clean2d/` (la campaña de comportamiento geométrico de Clean2D
continúa por separado) ni se inicia Fase 7.

## Inconsistencias documentales corregidas en el mismo commit

- Se añadió `baseline.md` (faltante) con los datos reales registrados
  (suite completa `3c8d992` y baseline pre-Fase 6 `02c8ea8`).
- Se marcaron `tasks.md` 8.2 y 8.3 como completados.
- Este `AGENT_REPORT.md` deja de afirmar que el push está bloqueado (el usuario
  lo realizó; `origin/ui/modernization` estaba en `3c8d992` antes de esta
  corrección).
- `proposal.md`, `design.md` (decisión D4) y `specs/ui-command-palette/spec.md`
  se actualizan: la premisa "migración deliberada de `Ctrl+K`" se sustituye por
  "`Ctrl+K` para Clean2D quick, `Ctrl+P` para la paleta".

## Fallback de suite completa: 5 fallos RDKit (clase preexistente, flaky)

La suite completa (1838 tests = baseline 1811 + 27 nuevos de la paleta) devolvió
**5 failed, 1813 passed, 20 skipped**. Los 5 fallos están en el dominio
químico RDKit/async y **no importan nada de lo que esta fase tocó**:

- `test_smiles_stereo_import.py::test_chiral_smiles_import_creates_wedge_or_hash`
- `test_smiles_stereo_import.py::test_amino_acid_chiral_smiles_import_creates_wedge_or_hash`
- `test_smiles_stereo_import.py::test_tetrandrine_import_preserves_visual_stereo`
- `test_clean2d_engine_candidates.py::test_generate_candidates_attempts_rdkit_for_cyclic_graphs`
- `test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend`

El baseline documentado tenía 4 fallos RDKit de esta clase. El de
`compchem3d` **pasa al re-ejecutarse en aislamiento** (es flaky/async). Ninguno
depende de la paleta, el AppBar ni de los atajos. No se modifica ninguna química
(Clean2D/ChemName/serialización) en esta fase.

## Verificación (toda ejecutada, salida real)

- `compileall` src tests tools packaging → exit 0
- Ruff scoped (F401,F811,F821,E722,E741) sobre los archivos de la fase → All
  checks passed
- `git diff --check` → limpio
- `openspec validate 2026-10-01-modernize-ui-command-palette --strict` → valid
- `tests/architecture/` → 276 passed (con `command_palette.py` y
  `command_registry.py` registrados en M08 de `modules.yml`)
- Tests dirigidos: `test_command_palette.py` (27), `test_ui_app_bar_tabs.py`
  (24), `test_main_window_tabs.py` (15), `test_branch_rotation_shortcuts.py`
  (7), `test_ui_tool_rail.py` (64) → **137 passed**
- Evidencia offscreen: light/dark × 1440×900 (query vacío, "export", "valid")
  y 980×600, en `/tmp/baseline/evidence/`.

### Comprobación de atajos (salida real, ventana real)

- `Ctrl+K` → `['Limpiar 2D (1 paso)']` (restaurado a Clean2D quick)
- `Ctrl+P` → `['Buscar o ejecutar…']` (paleta; única `QAction` global que lo
  posee)
- `Ctrl+Shift+K` → `['Limpiar 2D para publicación']` (intacto)
- `Ctrl+Alt+K` → `['Proponer conformero 2D']` (intacto)
- `action_clean_2d_full.shortcutContext()` == `WindowShortcut`; sigue presente
  en menú *Estructura* y en el registro de la paleta.

## Decisiones de arquitectura

- `command_palette.py` (widget de presentación) no importa dominios; la
  construcción del registro vive en `command_registry.py` (boundary limpio).
- La paleta es hija de `self` (ventana) → overlay a pantalla completa, igual que
  el spike aprobado.
- La `QAction` existente es la fuente de verdad; la paleta indexa y ejecuta con
  `action.trigger()`, sin crear acciones duplicadas.

## Bug corregido durante la fase (evidencia visual)

En `_rebuild()` los headers de sección se limpiaban solo con `deleteLater()`
(asíncrono): al filtrar, cabía que los 9 headers viejos se pintaran junto a la
lista nueva hasta que el event loop procesara la eliminación. Se añadió
`setParent(None)` síncrono antes del `deleteLater()`. Verificado: tras
`"export"` solo aparecen las 5 secciones/8 filas correctas (antes 14 headers).
