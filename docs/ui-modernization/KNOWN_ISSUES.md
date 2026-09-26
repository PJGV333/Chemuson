# KNOWN ISSUES — Modernización de UI (ChemUSON)

Deudas funcionales confirmadas durante la validación manual de la UI
modernizada (Fases 1–4). Este archivo registra problemas **no resueltos** que
se documentan para repararse en campañas/tareas separadas; NO bloquean la
modernización visual.

---

## Branch rotation actions do not execute visibly

Estado:
CONFIRMADO POR VALIDACIÓN MANUAL

Afecta:

- Girar rama -60° (`Ctrl+Alt+Left`)
- Girar rama +60° (`Ctrl+Alt+Right`)

Ruta UI:

Editar -> Rotar

Acciones (definidas en `src/chemuson/gui/actions/structure_actions.py`):

- `window.action_branch_rotate_minus`
- `window.action_branch_rotate_plus`

Handlers:

- `window._on_rotate_branch(-BRANCH_ROTATION_STEP_DEG)`
  (`src/chemuson/gui/main_window.py`, método `_on_rotate_branch`)
- `window._on_rotate_branch(+BRANCH_ROTATION_STEP_DEG)`

Comportamiento esperado:
con una rama/selección válida, la rama debe rotar ±60°
(`BRANCH_ROTATION_STEP_DEG = 60`).

Comportamiento observado:
ni el QAction del menú ni su shortcut (`Ctrl+Alt+Left` / `Ctrl+Alt+Right`)
producen una transformación visible en el lienzo.

Alcance:
problema funcional preexistente o independiente del rediseño visual;
NO es una regresión introducida por la modernización de UI;
NO bloquear Fase 5 de modernización de UI.

Acciones relacionadas que deben verificarse cuando se investigue este issue
(NO marcadas como rotas; sin evidencia):

- Invertir rama (180°) — `Ctrl+Alt+I` (`action_branch_invert` → `_on_invert_branch`)
- Autoacomodar rama — `Ctrl+Alt+A` (`action_branch_auto_arrange`)
- Fragmento con pivote (rotación de fragmento, `FRAGMENT_ROTATION_STEP_DEG`)

Pendiente:
investigar en campaña/tarea separada después de estabilizar el port visual.

---

## Escopo: Clean2D queda fuera de la campaña de UI

Clean2D (`src/chemuson/clean2d/`) **NO se tratará en esta campaña de
modernización de UI**. Existe una campaña independiente activa para Clean2D
(política de candidados y layout); durante la modernización de UI no se
abre ni modifica `clean2d/`.
