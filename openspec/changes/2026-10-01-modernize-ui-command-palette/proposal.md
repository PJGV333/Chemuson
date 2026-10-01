# Proposal: Paleta de comandos Ctrl+K (Fase 6)

## Why

Tras las fases 1–5, la región superior y lateral de la UI está consolidada
(app bar, pestañas, rail, side panel, barra de estado). Falta la mayor ganancia
de "práctico" por esfuerzo según `docs/ui-modernization/PLAN.md` §2.4/§Fase 6:
una **paleta de comandos** (Ctrl+K) que haga buscables todas las acciones de la
aplicación (archivos, edición, vista, estructura, análisis, exportaciones,
paneles laterales, plantillas, preferencias/tema) desde un único punto, con la
píldora de búsqueda del AppBar como entrada visible.

La píldora de búsqueda del AppBar es hoy un placeholder deliberado (Fase 3)
que no registra atajo alguno; su badge `Ctrl K` está oculto porque `Ctrl+K`
pertenece a `action_clean_2d_full`. Esta fase convierte la píldora en entrada
real a la paleta y resuelve de forma deliberada ese conflicto de atajo.

## What Changes

- Nuevo `src/chemuson/gui/command_palette.py`: `CommandEntry` (presentación:
  sección, keywords, icono), `CommandRegistry` (registro + deduplicación por
  identidad de `QAction`) y `CommandPalette` (overlay centrado sobre la
  ventana: input + lista filtrada + navegación ↑/↓ + Enter/Esc + clic).
- **La `QAction` existente es la fuente de verdad**: la paleta presenta, filtra
  y ejecuta `QAction` ya existentes (mantiene identidad, enabled,
  checkable/checked, shortcut mostrado, icono, `triggered` y handlers). No
  genera una segunda `QAction` por comando ni duplica handlers ni lógica
  química.
- **Migración deliberada de `Ctrl+K`**: `action_clean_2d_full` conserva su
  `QAction`, su handler y su función, pero **deja de poseer** `Ctrl+K` (se
  quita el shortcut; la acción sigue accesible desde
  *Estructura → Limpiar 2D (1 paso)* y desde la paleta). `Ctrl+K` pasa a ser el
  shortcut global de una única `action_command_palette` (`QKeySequence("Ctrl+K")`,
  `WindowShortcut`, `window.addAction(...)`), que también es el camino de
  apertura de la píldora de búsqueda (mismos dos caminos → una misma acción).
  `Ctrl+Shift+K` (publicación) y `Ctrl+Alt+K` (conformero) quedan intactos.
  No se inventa otro shortcut para `action_clean_2d_full` en esta fase.
- `shell/assembly.py` / `main_window.py`: construir una única
  `CommandPalette` (registrando las fuentes reales de comandos) y la única
  `QAction` de apertura `action_command_palette` (Ctrl+K) conectada a la paleta
  y a `app_bar.search_pill`.
- `app_bar.py` (mínimo): la `SearchPill` deja de ser placeholder — el badge
  `Ctrl K` se muestra y el tooltip deja de decir "(próximamente)"; su clic emite
  `activated` que la ventana conecta al MISMO camino de apertura (la `QAction`
  de Ctrl+K). La SearchPill no se convierte en editor permanente: el campo
  editable pertenece a `CommandPalette`.
- **Fuente de plantillas**: se reutiliza el contrato `QAction` ya existente del
  menú dinámico de plantillas (`window.templates_menu`, que
  `TemplateBrowserService.refresh_templates_menu` reconstruye). La paleta
  refresca sus entradas de plantilla desde ese menú cuando la biblioteca
  cambia; si no existe contrato limpio, se usa un adaptador mínimo apoyado en
  `template_browser_service` (ver design D5). No se modifica el contenido
  químico de las plantillas.
- `theme/tokens.py` + `theme/qss.py`: métrica de ancho (`paletteW`) y estilos
  tokenizados del overlay/tarjeta/input/sections/rows/selected (light + dark).
- `architecture/modules.yml`: registro de `command_palette.py` en M08 (gui).
- Tests dirigidos en `tests/test_command_palette.py` y ajuste **solo** de los
  tests/documentación que afirmaban que `Ctrl+K` ejecuta Clean2D quick.

## No Changes

- No se modifica `src/chemuson/clean2d/`, `src/chemuson/chemname/`,
  `src/chemuson/chemio/persistence.py`, `src/chemuson/gui/canvas/` ni
  `src/chemuson/gui/editor2d/` (guardarraíles §13).
- No se duplica Clean2D ni su lógica; la paleta no llama a Clean2D directamente.
- No se altera geometría molecular, selección, hit-testing, orbitales,
  diagramas, `tool_id`, comportamiento de SidePanel ni la estructura visual
  aprobada de AppBar/rail/canvas. No se rediseña el AppBar.
- No se agregan nuevas dependencias externas (fuzzy-search, etc.).
- No se agregan nuevas preferencias de persistencia para "último comando"
  (queda local a la sesión, ver design D3).
- No se inicia Fase 7; no se reabren las fases 1–5. Los defectos conocidos de
  Plantillas quedan diferidos a polish/hardening posterior.
