# Proposal: Paleta de comandos Ctrl+P (Fase 6)

## Why

Tras las fases 1–5, la región superior y lateral de la UI está consolidada
(app bar, pestañas, rail, side panel, barra de estado). Falta la mayor ganancia
de "práctico" por esfuerzo según `docs/ui-modernization/PLAN.md` §2.4/§Fase 6:
una **paleta de comandos** (Ctrl+P) que haga buscables todas las acciones de la
aplicación (archivos, edición, vista, estructura, análisis, exportaciones,
paneles laterales, plantillas, preferencias/tema) desde un único punto, con la
píldora de búsqueda del AppBar como entrada visible.

La píldora de búsqueda del AppBar es hoy un placeholder deliberado (Fase 3)
que no registra atajo alguno. Esta fase convierte la píldora en entrada real a
la paleta. La paleta usa `Ctrl+P` (verificado libre de conflictos en
producción); `Ctrl+K` se conserva como el atajo histórico de
`action_clean_2d_full` (limpia 2D, 1 paso).

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
- **Atajo de la paleta: `Ctrl+P`; `Ctrl+K` se conserva para Clean2D quick**:
  `action_command_palette` es la única `QAction` global que posee
  `QKeySequence("Ctrl+P")` (verificado libre de conflictos en producción
  antes de asignarlo), con `WindowShortcut` y `window.addAction(...)`. La
  píldora de búsqueda comparte ese mismo camino de apertura. `Ctrl+K` es el
  atajo histórico de `action_clean_2d_full` (limpia 2D, 1 paso), sin tocar su
  handler ni su lógica: se conserva `setShortcut("Ctrl+K")` +
  `window.addAction(...)`. `Ctrl+Shift+K` (publicación) y `Ctrl+Alt+K`
  (conformero) quedan intactos. (Nota: la implementación original de la Fase 6
  migró `Ctrl+K` a la paleta; la corrección de UX lo restaura a Clean2D quick
  y asigna `Ctrl+P` a la paleta.)
- `shell/assembly.py` / `main_window.py`: construir una única
  `CommandPalette` (registrando las fuentes reales de comandos) y la única
  `QAction` de apertura `action_command_palette` (Ctrl+P) conectada a la paleta
  y a `app_bar.search_pill`.
- `app_bar.py` (mínimo): la `SearchPill` deja de ser placeholder — el badge
  `Ctrl P` se muestra y el tooltip deja de decir "(próximamente)"; su clic emite
  `activated` que la ventana conecta al MISMO camino de apertura (la `QAction`
  de Ctrl+P). La SearchPill no se convierte en editor permanente: el campo
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
- Tests dirigidos en `tests/test_command_palette.py` (y ajuste de los tests de
  AppBar/shortcuts) que demuestran: Ctrl+K ejecuta `action_clean_2d_full` y no
  abre la paleta; Ctrl+P abre la paleta; la SearchPill abre la misma ruta; solo
  una `QAction` global posee Ctrl+P; Ctrl+Shift+K/Ctrl+Alt+K intactos; Clean2D
  quick sigue accesible por menú y por la paleta.

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
