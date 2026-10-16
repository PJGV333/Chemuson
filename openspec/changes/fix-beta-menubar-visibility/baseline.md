# Baseline — opción de visibilidad de la barra superior

- Rama: `fix/beta-final-ui-clean2d`.
- HEAD local/remoto: `94ffdd80080c9904a4a2d472b06a5174f65d8630`.
- `git status --short`: sin salida; árbol limpio antes de editar.
- `python -m compileall src tests tools packaging` (timeout 300 s): pasó.
- `pytest --collect-only -q` (timeout 600 s): `2236 tests collected in 0.62s`.
- Baseline completo `pytest -q` (timeout externo 600 s): terminó por límite (`timeout` exit 124) con progreso hasta 38%, sin fallo de test visible; no se repitió el monolítico. No quedaron procesos pytest/RDKit en ejecución. Suite se verificará por shards bounded.
- Regresión UI/preferencias baseline: `pytest -q tests/test_beta_compact_shell_and_tabs.py tests/test_ui_app_bar_tabs.py tests/test_ui_tool_rail.py tests/test_platform_settings.py`: `127 passed in 165.55s`.
- `pytest -q tests/architecture`: `280 passed in 9.59s`.
- Ruff focal `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`: un único finding conocido ajeno `F401 math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.
- Wires investigados: `action_show_menu_bar` controla `QMenuBar` y persiste `ui/show_menu_bar`; el native `main_toolbar.toggleViewAction()` se presenta con el mismo label solicitado y controla `QToolBar`/`ui/show_main_toolbar`; AppBar hamburguesa abre `_menu_aggregator` con los mismos seis `QMenu`. `show_menu_bar` default es true; la toolbar default false, pero QSettings v2 puede conservar `show_main_toolbar=true`.
