# Diseño — visibilidad de la fila de menús

## Causa y wiring actual

`_build_view_menu()` crea `action_show_menu_bar` («Mostrar barra de menús») y la conecta a `_on_toggle_menu_bar`, que muestra/oculta la `QMenuBar` y persiste `ui/show_menu_bar`. Separadamente, `build_main_toolbar()` añade a Ver `main_toolbar.toggleViewAction()` bajo el texto «Mostrar barra de herramientas»; esa QAction se conecta a `_on_toggle_main_toolbar` y persiste `ui/show_main_toolbar`. De ahí la discrepancia de superficie. La AppBar ya expone una hamburguesa que abre `_menu_aggregator`, construido con los mismos seis `QMenu`, de modo que Ver y su acción checkable siguen disponibles mientras la barra clásica está oculta.

## Decisión mínima

- Reusar `action_show_menu_bar` y su callback/QSettings. Cambiar sólo su texto visible a «Mostrar barra de herramientas» y retirar el toggle nativo de la toolbar del menú Ver; no crear otra QAction ni duplicar handlers.
- En la construcción del toolbar antiguo, mantenerla montada para compatibilidad de sus acciones pero ocultarla incondicionalmente y deshabilitar su `toggleViewAction()` nativo, de modo que valores históricos `ui/show_main_toolbar=true` no la resuciten ni desde Ver ni desde el menú contextual de QMainWindow.
- No cambiar `ui/show_menu_bar`, defaults (ya `True`), `SHELL_DEFAULTS_VERSION` ni el wiring de `_on_toggle_menu_bar`; esta es la preferencia persistida que ya representa la fila que se oculta.
- Actualizar comentarios/docstrings y pruebas que afirmaban que la QAction visible controlaba la QToolBar. La AppBar, los QMenus/QActions de sus comandos, atajos y panel permanecen intactos.

## Archivos previstos

`src/chemuson/gui/main_window_ui_builder.py`, documentación de preferencias en `src/chemuson/platform/settings.py` si hace falta, tests UI/preferences relevantes, y este OpenSpec. No tocar Clean2D ni capas químicas.
