# Design: Paleta de comandos Ctrl+K (Fase 6)

## Resumen de decisiones

### D1 — La `QAction` existente es la fuente de verdad (identidad, no texto)

La paleta **no** crea `QAction` por comando. Cada entrada del registro apunta a
una `QAction` real ya viva en la ventana (menús, toolbars ocultos, app bar).
Ejecutar una entrada llama a `action.trigger()` (o `action.activate`), de modo
que identidad, `enabled`, `checkable`/`checked`, icono, shortcut mostrado,
handlers y semántica son exactamente los de la acción histórica. La paleta
solo añade *metadata de presentación* (`section`, `keywords`, `icon`) que no
copia estado funcional.

### D2 — Deduplicación por identidad de `QAction`

`CommandRegistry.register(action, section, keywords)` indexa por
`id(action)` (identidad de objeto), no por texto. Registrar la misma `QAction`
dos veces (p. ej. desde dos fuentes) NO crea dos entradas: la segunda
inscripción es no-op (conserva la primera sección). Un test verifica que no
existe ninguna `QAction` representada dos veces y que reabrir la paleta no
duplica registros (el registro se construye una sola vez al montar la ventana).

### D3 — "Último comando usado" local a la sesión

El "último resultado usado como atajo rápido" se implementa **solo en memoria**
de la instancia de `CommandPalette` (no se persiste a `platform.settings`). Con
query vacío, la fila 0 es el último comando ejecutado en la sesión si existe,
seguido del catálogo ordenado. No se agregan nuevas preferencias: no hay razón
arquitectónica para persistir una preferencia que solo afecta la ordenación con
query vacío y cuyo valor es trivialmente recalculable.

### D4 — Migración deliberada de `Ctrl+K` (conflicto real resuelto)

- `action_clean_2d_full` **deja** de poseer `Ctrl+K`: se elimina
  `setShortcut(QKeySequence("Ctrl+K"))` y `window.addAction(...)`. Su `QAction`,
  su conexión `triggered → _on_clean_2d_full` y su item de menú
  *Estructura → Limpiar 2D (1 paso)* permanecen intactos. **No** se le asigna
  otro shortcut en esta fase.
- Se crea **una única** `QAction` de apertura `action_command_palette` con
  `QKeySequence("Ctrl+K")`, `Qt.ShortcutContext.WindowShortcut`, y
  `window.addAction(action_command_palette)`. Su `triggered` abre la paleta.
- La `SearchPill` del AppBar no registra su propio `QShortcut`; su clic emite
  `activated` que la ventana conecta a `action_command_palette.trigger()`.
  Así **Ctrl+K y la píldora comparten un mismo camino de apertura** (una misma
  `QAction`), cumpliendo "no implementes dos caminos".
- `Ctrl+Shift+K` (`action_clean_2d_publication`) y `Ctrl+Alt+K`
  (`action_clean_2d_propose`) no se tocan.
- Es una migración deliberada de atajo, no una regresión. Se actualizan
  únicamente los tests/documentación que afirmaban `Ctrl+K → Clean2D quick`.

### D5 — Plantillas vía contrato `QAction` existente (menú dinámico)

`TemplateBrowserService.refresh_templates_menu` ya crea una `QAction` real por
plantilla (conexión `triggered → start_template_insert_by_id(tid)`) en
`window.templates_menu`. La paleta **reutiliza esas mismas `QAction`**: al
refrescar sus entradas de plantilla lee `templates_menu.actions()` (recursivo
por submenús) y las registra con sección "Plantillas". Así no se inventa un
segundo sistema de inserción química; la paleta dispara la `QAction` existente,
que a su vez llama al contrato público `start_template_insert_by_id`.

Adaptador mínimo (solo si el menú dinámico no expone un contrato limpio): la
paleta expone un `refresh_templates()` que, en ausencia de `QAction` por
plantilla, crea entradas ligeras **presentacionales** que conectan su `triggered`
a `window._start_template_insert_by_id(tid)` leyendo la biblioteca
(`template_library.grouped_templates()`) — pero esto NO duplica la inserción:
delega en el mismo `TemplateController.start_template_insert_by_id`. Se prefiere
la vía `templates_menu` (D5 principal) y el adaptador queda documentado como
fallback. En esta implementación se usa la vía `templates_menu`.

### D6 — Overlay centrado sobre la ventana (fiel al spike), no diálogo modal

Como el spike aprobado (`pyqt6-spike/palette.py`), la paleta es un `QWidget`
hijo que cubre la ventana con una tarjeta centrada, capturando el teclado vía
`keyPressEvent` y un `QEventFilter` para el Esc/clic. Ancho objetivo
`paletteW` (560 px) limitado a `min(560, 92% de la ventana)` para no salirse en
ventanas estrechas; altura acotada con scroll. No es `QDialog` (evita el diálogo
del SO y es fiel al mockup "backdrop sobre la app").

## Ranking y filtro

- **Matching**: case-insensitive. Para cada entrada se evalúa el `title` y los
  `keywords` (y la `section`).
- **Prefijo > substring**: una entrada cuyo `title` (o keyword) **empieza por**
  el query puntúa más alto que una que solo lo contiene como substring. Ranking
  simple y determinista:
  - `tier 0` = prefijo del título;
  - `tier 1` = prefijo de keyword/section;
  - `tier 2` = substring en título;
  - `tier 3` = substring en keyword/section.
  - Dentro de cada tier, orden estable por (sección declarada, título) —
    determinista y testeable (test "prefijo antes que substring").
- Sin fuzzy; cero dependencias.

## Teclado / interacción

- `open()`: centra sobre la ventana, foco inmediato al `QLineEdit`, input vacío
  (o último comando en fila 0).
- `↑`/`↓` mueven la selección (envolvente), resaltada con `accentSoft`/`accent`.
- `Enter` ejecuta la fila seleccionada **una vez** y cierra la paleta.
- `Esc` cierra sin ejecutar.
- Clic en una fila ejecuta esa fila y cierra.
- Fila de `QAction` disabled: mostrada atenuada, **no ejecutable** (Enter/clic la
  saltan o la ignoran); la selección visible refleja el estado disabled.
- `checkable` conserva semántica: ejecutar una acción checkable alterna su
  `checked` (por delegación a `trigger()`, no por copiar el estado).
- Los atajos de herramienta de letra simple (`ToolShortcutDispatcher`) se
  suprimen cuando el foco está en el `QLineEdit` (comportamiento ya existente:
  el dispatcher cede ante widgets de entrada de texto).

## Arquitectura / fronteras

- `command_palette.py` (nuevo, dentro de M08 `gui`) depende de `QAction`,
  metadatos de presentación, `theme` (tokens/QSS/IconProvider) y el shell
  (ventana). **NO importa** `chemuson.clean2d`, `chemuson.chemname`,
  `chemuson.chemio.persistence`, ni internals de `gui.canvas` ni controllers
  químicos concretos (contrato AST, ver spec).
- Registro de comandos alimentado desde `main_window`/`shell` (que sí conoce
  todas las `QAction`); la paleta recibe las `QAction` ya creadas.
