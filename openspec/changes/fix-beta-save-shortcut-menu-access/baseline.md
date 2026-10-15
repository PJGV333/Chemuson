# Baseline — Ctrl+S, cierre de pestañas y menús visibles

Captura previa a cambios de producto. Rama solicitada: `fix/beta-save-shortcut-menu-access`, creada desde `release/v0.3.0-beta.1-prep`.

## Refs y árbol

- Rama base local/remota: `release/v0.3.0-beta.1-prep` = `origin/release/v0.3.0-beta.1-prep` = `9bc8192db3835262191ddb5247a8ff87f4d9169f`.
- Rama correctiva recién creada: `fix/beta-save-shortcut-menu-access` en el mismo SHA `9bc8192db3835262191ddb5247a8ff87f4d9169f`; no existía una rama correctiva remota.
- `git status --short --branch` antes de la rama: `## release/v0.3.0-beta.1-prep...origin/release/v0.3.0-beta.1-prep` (limpio). Después: `## fix/beta-save-shortcut-menu-access` (limpio).
- No se modificaron beta, `main`, stable ni tags.

## Baseline de herramientas y pruebas

Ejecuciones locales con timeout externo; ninguna invocación de pruebas excedió diez minutos:

| Comando | Resultado exacto resumido | Estado |
|---|---|---|
| `python -m compileall src tests tools packaging` | `exit=0` | PASS |
| `python -m pytest --collect-only -q` | `2154 tests collected in 0.80s` | PASS |
| `python -m pytest -q` bajo `timeout --signal=TERM --kill-after=10s 480s` | `exit=124` al alcanzar 480 s, con progreso hasta ~90%; un `F` en el nodo 1098 y el shard focalizado UI posterior pasa. No se considera suite completa ni PASS. | INCOMPLETA / TIMEOUT |
| `python -m pytest -q tests/test_ui_polish.py tests/test_main_window_tabs.py tests/test_ui_app_bar_tabs.py` bajo timeout 240 s | `78 passed in 40.37s` | PASS |
| `python -m pytest -q tests/architecture` bajo timeout 180 s | `280 passed in 10.01s` | PASS |
| `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` | `Found 1 error`: F401 `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3` | FAIL baseline conocido |
| `openspec validate prepare-chemuson-0-3-0-beta-1 --strict` | `Change 'prepare-chemuson-0-3-0-beta-1' is valid` | PASS |

La suite completa produjo el nodo histórico `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend` (posición 1098/2154), que se ejecutó aisladamente y pasó (`1 passed in 0.16s`); no se atribuye una causa nueva. Llegó al nodo 1945, `tests/test_ui_polish.py::test_onboarding_mask_is_path_subtraction_and_hole_is_transparent`, antes de agotar el límite; ese nodo aislado pasó (`1 passed in 0.37s`) y el shard de UI completo pasó. No se repitió la suite monolítica.

## Diagnóstico funcional reproducible previo al fix

La prueba se hizo con una ventana Qt mostrada/activa, el canvas con estado dirty, un archivo `.cmsn` existente en directorio temporal y un átomo añadido al contenido. Se leyó el JSON guardado, además de observar la `QAction` y `QUndoStack`:

| Interacción | Acción invocada | Átomos serializados | Estado dirty |
|---|---|---:|---|
| `Ctrl+S` (QTest sobre canvas enfocado; menú compacto/oculto) | ninguna | 0 (el archivo inicial queda intacto) | permanece dirty |
| `Ctrl+G` (misma configuración) | ninguna | 0 | permanece dirty |
| Ejecutar `window.action_save` desde Archivo/acción original | `window.action_save` | 1 | clean |
| Repetir `Ctrl+S` tras asociar **ese mismo** `QAction` también a la ventana mediante `QMainWindow.addAction` | `window.action_save` | 1 | clean |
| `Ctrl+S` con `QMenuBar` visible | `window.action_save` | 1 | clean |

La interacción se ejecutó con `QApplication.setActiveWindow(window)` y foco explícito en el canvas, por lo que no depende de la limitación de activación de ventanas del plugin offscreen. La hipótesis queda confirmada: el atajo pertenece a una `QAction` alojada en menús/barra que permanecen ocultos en el modo compacto; el action manager no registra la combinación en la ventana activa. La acción sí conecta al handler existente `_on_file_save`, cuyo `FileController.save_file` conserva el flujo CMSN. No se necesita un segundo manejador ni otro atajo. `Ctrl+G` está presente en la misma acción, pero también queda inoperante por la misma asociación oculta; se conservará.

## Cierre y alcance observado

- `DocumentController.confirm_discard_changes` consulta el `QUndoStack` del canvas solicitado. Para Save activa temporalmente ese documento, ejecuta el guardado existente y sólo autoriza cerrar si `canvas.undo_stack.isClean()`; un cancel/error que deje dirty devuelve falso. Discard requiere `QMessageBox.StandardButton.Discard`; Cancel mantiene la pestaña. La ruta `.mol` está separada del guardado nativo.
- `DocumentTabBar` crea el botón X de 18×18 en cada pestaña y emite `closeRequested` hacia el flujo existente. El botón tiene test unitario de señal, pero no había cobertura integrada de Save/Discard/Cancel ni prueba de contenido CMSN al usar Ctrl+S.
- El `QMenuBar` y las seis `QMenu` existentes se crean en `MainWindowUiBuilder`; `assembly.py` oculta el menubar y el botón hamburguesa abre el agregador de esos menús. El menú Ver no tenía preferencia de visibilidad persistida. `platform.settings.UiPreferences` ya es la política existente para preferencias UI (`ui/theme`).

## Restricciones de implementación

Sólo adaptadores GUI/settings y tests/OpenSpec. No tocar ChemName, Clean2D, ChemIO, química/estereoquímica, persistencia interna `.cmsn`, command palette, rail de herramientas ni reglas de publicación/release.