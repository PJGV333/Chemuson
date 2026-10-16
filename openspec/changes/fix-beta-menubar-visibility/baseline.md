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

## Baseline de seguimiento: geometría del popup Ver (HEAD `9c4d08283894b2ebbb1bc9793a77ded46b6c0656`)

- Rama `fix/beta-final-ui-clean2d`; `git status --short` vacío. Desde el commit de implementación previo `d4d2f67` hasta este HEAD sólo se habían actualizado estos documentos OpenSpec; producción/tests seguían idénticos.
- `python -m compileall src tests tools packaging`: pasó.
- `pytest --collect-only -q`: `2235 tests collected in 0.66s`.
- `pytest -q tests/architecture`: `280 passed in 12.30s`.
- Baseline focal previa sobre el código idéntico: `pytest -q tests/test_beta_compact_shell_and_tabs.py tests/test_beta_save_shortcut_menu_access.py`: `45 passed in 55.68s`.
- Baseline heredada de la suite monolítica `pytest -q`: timeout externo 600 s al 38% sin fallo visible; no se repitió por el límite ya observado. El SHA fuente/test es el mismo que el de la CI previa verde; se validará la modificación nueva por grupos acotados y shards CI.
- Ruff global con `F401,F811,F821,E722,E741`: único finding de baseline, `F401 math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.
- Prueba diagnóstica Qt offscreen antes del cambio: `QScreen.availableGeometry=(0,0,800,800)`; el popup Ver medía `778×786` px, exponía 22 opciones directas en dos columnas y casi consumía toda la altura. La captura de diagnóstico quedó fuera del repositorio en `/tmp/chemuson-view-menu-before.png`. El área QScreen offscreen es fija en este entorno; los casos de 1600×900, 1280×720 y 1024×768 se parametrizarán como tamaños de ventana, verificando las coordenadas reales del popup contra el área disponible.
- Una prueba ad hoc de mouse ocultó/restauró la fila correctamente, pero el primer proceso salió 139 al destruirse con los QMenus emergentes aún abiertos. Repetir cerrando explícitamente `view_menu`/`_menu_aggregator` y drenando deletes diferidos terminó limpio; las regresiones deben cerrar popups en `finally`. No hubo cambio de código de lifecycle ni fallo en la suite pytest.
