# KNOWN ISSUES — Modernización de UI (ChemUSON)

Deudas funcionales confirmadas durante la validación manual de la UI
modernizada (Fases 1–4). Este archivo registra problemas y su estado; los
resueltos quedan documentados con su commit de referencia; los abiertos se
reparan en campañas/tareas separadas y NO bloquean la modernización visual
ni la Fase 5.

---

## Branch rotation shortcuts were inactive with the menu closed

Estado:
**RESUELTO** — commit `68c0051` (`Fix branch rotation window shortcuts`,
rama `ui/modernization`)

Afectaba (originalmente confirmado por validación manual):

- Girar rama -60° (`Ctrl+Alt+Left`)
- Girar rama +60° (`Ctrl+Alt+Right`)
- Invertir rama (180°) — `Ctrl+Alt+I` (`action_branch_invert`)
- Autoacomodar rama — `Ctrl+Alt+A` (`action_branch_auto_arrange`)

Ruta UI:

Editar -> Rotar

Acciones (definidas en `src/chemuson/gui/actions/structure_actions.py`):

- `window.action_branch_rotate_minus`
- `window.action_branch_rotate_plus`
- `window.action_branch_invert`
- `window.action_branch_auto_arrange`

Handlers (funcionaban; no se modificaron):

- `window._on_rotate_branch(±BRANCH_ROTATION_STEP_DEG)`
- `window._on_invert_branch()`
- `window._on_auto_arrange_branch()`

Diagnóstico (demostrado con experimento y tests, no asumido):

- La ruta de menú funcionaba: el mismo QAction vivo en
  ``Editar -> Rotar`` rotaba la rama al activarse.
- El defecto era **solo de wiring del shortcut**: las QActions estaban
  asociadas únicamente al QMenu de la barra de menús, y la UI moderna deja
  la barra oculta (`menuBar().setVisible(False)`). Con la barra oculta, el
  shortcut map de esos QMenus NO está activo cuando el foco está en el
  lienzo (o en cualquier hijo de la ventana), por lo que
  ``Ctrl+Alt+Left`` / ``Ctrl+Alt+Right`` caían en el nudge de flechas del
  canvas (traslación de 1 px, imperceptible) y ``Ctrl+Alt+I`` /
  ``Ctrl+Alt+A`` caían en el vacío.
- Las acciones de Clean2D (``Ctrl+K`` etc.) funcionan globalmente porque
  además se registran en la ventana con
  ``setShortcutContext(Qt.ShortcutContext.WindowShortcut)`` +
  ``window.addAction(...)``.

Arreglo (patrón Qt correcto, sin duplicar acciones ni conexiones):

- A las cuatro QActions históricas (las mismas que viven en el menú
  ``Rotar``) se les aplicó
  ``setShortcutContext(Qt.ShortcutContext.WindowShortcut)`` +
  ``window.addAction(action)``.
- Los shortcuts ahora se activan siempre que ``ChemUSONWindow`` tiene foco
  (lienzo o cualquier hijo), sin necesidad de abrir el menú.
- La ruta de menú sigue funcionando con los mismos QActions.
- La política existente de supresión de shortcuts se mantiene: las
  combinaciones Ctrl+Alt+* no escriben en editores de texto y, por
  semántica estándar de ``WindowShortcut``, siguen activas con foco en un
  editor no modal (mismo comportamiento que el ``Ctrl+K`` de Clean2D);
  la supresión por foco en editor cubre las letras de herramienta vía
  ``_tool_shortcuts_suppressed``.

Verificación:

- `tests/test_branch_rotation_shortcuts.py` (nuevo): wiring en la ventana,
  E2E con `ChemusonWindow` real de los cuatro shortcuts con el menú
  cerrado (rotación ±60°, inversión 180°, autoarrange con obstáculo),
  undo/redo, ausencia de doble ejecución (1 step de undo por pulsación),
  ruta de menú intacta y semántica de foco en editor no modal.
- `tests/test_branch_reorientation.py` (canvas level): sin cambios, en verde.

Acción relacionada verificada en el mismo ejercicio:

- Fragmento con pivote (rotación de fragmento,
  `FRAGMENT_ROTATION_STEP_DEG`): fuera de este defecto; no se tocó.

Pendiente:
ninguno para este issue.

---

## Escopo: Clean2D queda fuera de la campaña de UI

Clean2D (`src/chemuson/clean2d/`) **NO se tratará en esta campaña de
modernización de UI**. Existe una campaña independiente activa para Clean2D
(política de candidados y layout); durante la modernización de UI no se
abre ni modifica `clean2d/`.
