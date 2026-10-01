# AGENT_REPORT — Fase 6 (paleta de comandos Ctrl+K)

Cambio OpenSpec: `2026-10-01-modernize-ui-command-palette` (strict valid).
Commit local: `cca441c` (rama `ui/modernization`). Worktree limpio.

## Desviación / bloqueo: push a `origin` no ejecutable

**Situación:** El trabajo está completo, verificado y comprometido localmente
(`cca441c`), pero **no puedo empujar** a `origin/ui/modernization`. El remote es
HTTPS (`https://github.com/PJGV333/Chemuson.git`) y no hay credencial/token
guardado en esta máquina:

```
fatal: could not read Username for 'https://github.com':
       No existe el dispositivo o la dirección
```

Esto es un bloqueo de **autenticación**, no de contenido. Para completarlo,
la persona que tiene acceso debe (una de):

- `git push origin ui/modernization` con su token/credencial, o
- añadir un token (GH CLI / credential helper / key) y repetir el push.

No hay `--force`: `origin/ui/modernization` está en `02c8ea8` y el commit local
`cca441c` lo extiende (1 ahead, 0 behind). Un push normal es fast-forward.

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
depende de la paleta, el AppBar, la migración de Ctrl+K ni de `command_registry`.
No se modifica ninguna química (Clean2D/ChemName/serialización) en esta fase.

## Verificación (toda ejecutada, salida real)

- `compileall` src tests tools packaging → exit 0
- Ruff scoped (F401,F811,F821,E722,E741) sobre los 14 archivos cambiados → All checks passed
- `git diff --check` → limpio
- `openspec validate 2026-10-01-modernize-ui-command-palette --strict` → valid
- `tests/architecture/` → 276 passed (con `command_palette.py` y
  `command_registry.py` registrados en M08 de `modules.yml`)
- Tests dirigidos de la fase → `test_command_palette.py` (27), AppBar, tool rail,
  main window tabs, branch rotation, clean2d safety → **137 passed**
- Evidencia offscreen: light/dark × 1440×900 (query vacío, "export", "valid")
  y 980×600, en `/tmp/baseline/evidence/`.

### Bug corregido durante la fase (evidencia visual)

En `_rebuild()` los headers de sección se limpiaban solo con `deleteLater()`
(asíncrono): al filtrar, cabía que los 9 headers viejos se pintaran junto a la
lista nueva hasta que el event loop procesara la eliminación. Se añadió
`setParent(None)` síncrono antes del `deleteLater()`. Verificado: tras
`"export"` solo aparecen las 5 secciones/8 filas correctas (antes 14 headers).

## Decisiones de arquitectura

- `command_palette.py` (widget de presentación) no importa dominios; la
  construcción del registro vive en `command_registry.py` (boundary limpio).
- La paleta es hija de `self` (ventana) → overlay a pantalla completa, igual que
  el spike aprobado.
- La `QAction` existente es la fuente de verdad; la paleta indexa y ejecuta con
  `action.trigger()`, sin crear acciones duplicadas.
