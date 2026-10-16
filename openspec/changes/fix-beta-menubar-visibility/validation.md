# Validación — visibilidad de la fila de menús

## Causa raíz y wiring

En `Ver` coexistían dos acciones checkable: `action_show_menu_bar` («Mostrar barra de menús») que controlaba `QMenuBar` y guardaba `ui/show_menu_bar`, y el `toggleViewAction()` de `main_toolbar` («Mostrar barra de herramientas») que mostraba Nuevo/Abrir/Guardar y guardaba `ui/show_main_toolbar`. El label que el propietario quería usar estaba conectado a la superficie equivocada. La hamburguesa ya abría `_menu_aggregator` con los mismos seis QMenu originales.

La corrección reaprovecha `action_show_menu_bar` y su persistencia existente, cambia su label al solicitado, y retira del menú Ver el toggle de la toolbar. La QToolBar legada se fuerza oculta y se deshabilita su toggle nativo, incluso si el INI mantiene `ui/show_main_toolbar=true`. No se cambia la preferencia de menú, los defaults (`show_menu_bar=True`) ni su versión de migración.

## Baseline

- Rama/HEAD: `fix/beta-final-ui-clean2d` / `94ffdd80080c9904a4a2d472b06a5174f65d8630`, igual a `origin`; árbol limpio antes del cambio.
- `python -m compileall src tests tools packaging` pasó.
- `pytest --collect-only -q`: `2236 tests collected in 0.62s` antes del cambio.
- Baseline `pytest -q`, con timeout externo 600 s: exit 124 al alcanzar 38% de progreso; no se mostró fallo de test. No se repitió el monolítico y no quedaron procesos pytest/RDKit activos.
- Baseline focal UI/settings: 127 passed in 165.55s (`test_beta_compact_shell_and_tabs`, `test_ui_app_bar_tabs`, `test_ui_tool_rail`, `test_platform_settings`). Arquitectura: 280 passed in 9.59s.
- Ruff global acotado: único hallazgo baseline conocido `F401 math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.
- Mapa previo: `action_show_menu_bar` → `QMenuBar`/`ui/show_menu_bar`; `main_toolbar.toggleViewAction()` → `QToolBar`/`ui/show_main_toolbar`; AppBar hamburguesa → agregador de los seis menús originales.

## Verificación local tras el cambio

- Primera pasada de regresión encontró cuatro aserciones antiguas que esperaban que «Mostrar barra de herramientas» mostrara la QToolBar. Se actualizaron para el nuevo contrato; no fue un fallo de runtime ni de acciones del documento.
- `pytest -q tests/test_beta_compact_shell_and_tabs.py`: 25 passed in 18.82s.
- Suite focalizada final, timeout 600 s: `pytest -q tests/test_beta_compact_shell_and_tabs.py tests/test_ui_app_bar_tabs.py tests/test_ui_tool_rail.py tests/test_platform_settings.py tests/test_beta_save_shortcut_menu_access.py tests/test_side_panel.py`: **161 passed in 309.35s**.
- `pytest -q tests/architecture`: 280 passed in 9.29s. `pytest -q tests/test_test_workflow_contract.py`: 4 passed in 0.06s.
- Ruff sobre todos los archivos GUI/tests modificados: pasó. Ruff global mantiene sólo el F401 baseline indicado arriba.
- `python -m compileall src tests tools packaging`: pasó. `openspec validate fix-beta-menubar-visibility --strict --no-interactive`: válido. `git diff --check`: pasó.
- Colecta final/plan: **2235 IDs**, 8 shards con distribución `[280,280,280,279,279,279,279,279]`, SHA-256 `7356386ef65f4686c066facb08da25d8a5aee2cb118c8e50fecc2e66bc9b6aad`. Se retiraron/reemplazaron node IDs de pruebas que afirmaban la visibilidad de la toolbar; el manifiesto fue actualizado con el updater tras revisar la diferencia.

## CI/aceptación

Pendiente de push y GitHub Actions por shards. La comprobación manual del propietario sigue pendiente; no integrar en `release/v0.3.0-beta.1-prep`, no modificar `main`, ni publicar tags/releases/paquetes.
