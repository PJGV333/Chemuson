# Baseline — X de pestañas y visibilidad de barra de herramientas

Fecha: 2026-10-15. Rama solicitada: `fix/beta-save-shortcut-menu-access`. Toda captura siguiente precede a cambios de esta corrección.

## Refs y árbol

- Base/HEAD inicial: `0765dcba7abf999e5f3c56efa5634d55bc2d2348`; `origin/fix/beta-save-shortcut-menu-access` coincide. Rama beta preservada: `origin/release/v0.3.0-beta.1-prep` = `9bc8192db3835262191ddb5247a8ff87f4d9169f`.
- `git status --short --branch`: `## fix/beta-save-shortcut-menu-access...origin/fix/beta-save-shortcut-menu-access` (sin cambios).
- `python -m compileall src tests tools packaging`: exit 0.
- `pytest --collect-only -q`: `2169 tests collected in 0.57s`, exit 0.
- `pytest -q` bajo `timeout --signal=TERM --kill-after=10s 480s`: alcanzó el límite externo, exit 124 tras 480 s (no PASS). Antes de interrumpirse emitió un `F` en el node ID 1113, `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend`; ejecutar sólo ese node inmediatamente después dio `1 passed in 0.18s`. Es un fallo transitorio/baseline fuera del alcance, observado antes de los cambios de UI y no se modifica.
- `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`: exit 1 por el F401 histórico `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no relacionado ni corregido.
- Focus UI antes de cambios: `pytest -q tests/test_ui_app_bar_tabs.py tests/test_beta_save_shortcut_menu_access.py tests/test_platform_settings.py tests/test_ui_tool_rail.py` bajo límite externo 300 s: `115 passed in 105.90s`.

## Reproducción de X ausente

Sonda local Qt offscreen, ventana mostrada, ejecución acotada a 120 s:

```text
initial icon_null= False icon_size= 11 11 button_size= 18 18 pixmap_null= False
new-tab icon_null= True pixmap_null= True toolTip= Cerrar pestaña menu= False toolbar_visible= False rail_visible= True side_panel_visible= False
after-dark-refresh new-tab-icon-null= False
EXIT:0
```

`DocumentTabBar.__init__()` refresca antes de que haya tabs; el primer tab recibe el `refresh_icons()` global que ejecuta `_apply_theme()` al final del shell. Luego `_add_tab()` instancia cada `_TabSideButton` sin fijar el icono de cierre. Un cambio posterior de tema refresca todas las side-buttons y hace aparecer la X. Causa reproducible: inicialización incompleta en `_add_tab()`/`_TabSideButton`, no el flujo de cierre. Además la regla QSS general de `QToolButton` da padding 8 px y la regla `tabClose` no lo sobreescribe; el botón era 18×18 con icono 11×11, dejando poco espacio efectivo.

`main_toolbar` se añade con sus QActions originales y al final de `gui/shell/assembly.py` se fuerza `setVisible(False)`. La acción actual “Mostrar copiar/pegar/zoom en barra superior” sólo inserta/remueve esas acciones auxiliares y no cambia la visibilidad del toolbar. `UiPreferences` sólo guarda `ui/show_menu_bar`, cuyo default es false; panel lateral y rail tienen controles/estado separados.

## Archivos de logs completos de baseline local

Las salidas crudas de captura de este entorno, anteriores a cambios, se conservaron durante la tarea en `/tmp/chemuson-interface-final-baseline/` (`compile-collect.log`, `pytest.log`, `ruff.log`, `ui-focused.log`, `icon-probe.log`, `failing-node.log`). Los comandos automáticos usaron límite externo ≤480 s; ninguna prueba Qt interactiva queda sin timeout.
