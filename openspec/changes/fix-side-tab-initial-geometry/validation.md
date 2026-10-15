# Validación — geometría inicial de SideTabRow

## Baseline

Rama `fix/beta-properties-panel-polish`, HEAD inicial `81f7d3763a520a23d0128cf2539660f975eacf58`, igual al remoto; árbol limpio antes de registrar este delta. Baseline focalizado: `tests/test_side_panel.py` 13 passed (5.22 s comando), compileall PASS (0.03 s), Ruff focalizado PASS, 2152 tests recogidos (0.82 s comando). No se ejecutó la suite completa, según el alcance pedido.

## Reproducción y causa raíz

- Se añadió primero la aserción de regresión, que falló antes del fix en la apertura desde **Ver**, sin click de tab: el botón Inspector medía 52 px y su `QToolButton.sizeHint()` con el stylesheet efectivo requería 62 px (`QFontMetrics` del texto 46 px + padding declarado 6 px). El strip inicial era 295 px, viewport 308 px y rango de scroll 0. Esto reproduce la geometría insuficiente del primer paint.
- En una sonda del mismo test se hizo `QTest.mouseClick` sobre Validación y se midió antes/después: ambas geometrías siguieron iguales. El cambio específico “un click arregla el texto” de las capturas del propietario no se reprodujo en este entorno Qt offscreen. Se registra explícitamente; la causa geométrica (ancho fijo inferior al sizeHint efectivo) sí es automática y la prueba nueva falla antes del fix.
- Causa: los anchos se fijaban en `__init__` con `QFontMetrics(text)+2*sideTabPadX`, antes del stylesheet/QStyle efectivos; `_strip.adjustSize()` y `minimumWidth` también se calculaban una sola vez, mientras el widget permanece oculto. El ancho calculado ignoraba los márgenes horizontales adicionales del estilo Qt. `set_active()` no es un hook de geometría fiable.

## Cambio y pruebas

`SideTabRow._refresh_tab_geometry()` mide cada botón tras `ensurePolished()` usando el máximo de `QFontMetrics+padding` y `QToolButton.sizeHint()`, actualiza anchos fijos y recalcula layout/tamaño/mínimo del strip en conjunto. Sólo cambia la geometría cuando el tamaño efectivo difiere, evitando transacciones de layout vacías durante cambios de tema. Se invoca al construir, mostrar, redimensionar y recibir cambios de estilo/fuente. Si el strip supera el viewport, se conserva el scroll existente y `ensureWidgetVisible`; no se altera QSplitter, preferencias ni el estado del panel.

| Comando | Resultado | Tiempo total |
|---|---:|---:|
| `pytest -q tests/test_side_panel.py` | 15 passed | 5.36 s |
| `pytest -q tests/test_ui_polish.py` | 38 passed | 10.32 s |
| `pytest -q tests/test_spectroscopy_dock.py` | 6 passed | 1.08 s |
| `pytest -q tests/test_ui_app_bar_tabs.py` | 24 passed | 5.01 s |
| `pytest -q tests/test_main_window_tabs.py` | 16 passed | 6.68 s |
| `pytest -q tests/architecture` | 280 passed | 10.44 s |
| `pytest -q tests/test_test_workflow_contract.py` | 4 passed | 0.21 s |
| `pytest --collect-only -q` | 2154 collected | 0.82 s |

Cada pytest tuvo timeout externo ≤180 s. Cobertura añadida: primera apertura desde Ver antes de interactuar con tabs; prueba de scroll/viewport insuficiente y todos los nombres visibles con ancho suficiente; hide/show, resize, light/dark y geometría del canvas sin segundo salto.

- `python -m compileall -q src tests tools packaging`: PASS (0.04 s).
- Ruff focalizado `F401,F811,F821,E722,E741` en Python cambiado: PASS.
- OpenSpec estricto para `fix-side-tab-initial-geometry` y `fix-beta-properties-panel-polish`: PASS.
- Contrato `tests/test_test_workflow_contract.py`: PASS; el manifiesto se actualizó sólo con los dos node IDs nuevos y confirma 2154 IDs.
- `git diff --check`: PASS.

## Riesgos / aceptación

La prueba confirma anchuras legibles por el estilo Qt y que, si no caben todos los botones, el tab elegido se desplaza completamente al viewport; al ensanchar el panel caben los cinco sin scroll. La estabilidad del canvas se mide contra la geometría tras la apertura inicial: el refresco no causa un segundo desplazamiento. La verificación automatizada no reemplaza el retest visual del propietario en el nuevo AppImage; ese estado manual sigue pendiente.

## Entrega remota

El primer commit `b639a1c` se publicó en `fix/beta-properties-panel-polish`. CI run [38672067445](https://github.com/PJGV333/Chemuson/actions/runs/38672067445) falló dos veces en el mismo test fuera de los archivos de este cambio: `tests/test_spectroscopy_dock.py::test_chemical_properties_scientific_display_is_presentation_only_and_theme_safe`, assertion `table.item(0, 1).isSelected()` después de `setCurrentCell(0, 1)`. Los otros siete shards, el plan completo, Windows smoke y Flatpak smoke pasaron en ambas ejecuciones. La prueba exacta pasó en 10 repeticiones consecutivas, en ejecuciones standalone adicionales y dentro del archivo `tests/test_spectroscopy_dock.py` completo; causa de la diferencia aún no determinada. Se redujo la geometría no-op durante `StyleChange` para mantener la rutina estrictamente idempotente; validar el commit de seguimiento en CI antes de cerrar la tarea.

No integrar a `release/v0.3.0-beta.1-prep` ni crear tag.
