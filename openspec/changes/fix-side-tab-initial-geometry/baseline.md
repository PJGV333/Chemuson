# Baseline — geometría de las pestañas al revelar SidePanel

- Rama: `fix/beta-properties-panel-polish`.
- HEAD inicial: `81f7d3763a520a23d0128cf2539660f975eacf58`; igual a `origin/fix/beta-properties-panel-polish`.
- `git status --short` inicial: limpio (antes de registrar este cambio OpenSpec).
- `python -m compileall -q src tests tools packaging` (timeout externo 180 s): exit 0, `ELAPSED=0.03 s`.
- `pytest --collect-only -q` (timeout externo 180 s): exit 0, 2152 tests collected; pytest 0.60 s / comando 0.82 s.
- `pytest -q tests/test_side_panel.py` (timeout externo 180 s): 13 passed in 5.05 s; comando 5.22 s. No existe todavía una aserción que mida el primer paint tras abrir Ver antes de la interacción con tabs; la prueba existente de legibilidad ejecuta `set_active()` antes de medir.
- `ruff check src/chemuson/gui/side_panel.py tests/test_side_panel.py --select F401,F811,F821,E722,E741` (timeout externo 180 s): PASS.
- La CI del HEAD de partida había verificado 2152 node IDs en 8 shards y las pruebas de panel; este baseline local se limita intencionalmente a pruebas focalizadas, según la solicitud.

## Inspección inicial del orden Qt

- `SideTabRow.__init__` calcula anchos fijos con `QFontMetrics(tab_font)` y `METRICS["sideTabPadX"]`.
- Antes de asignar `_strip` al `QScrollArea`, ejecuta `adjustSize()` y copia una vez `strip_layout.sizeHint().width()` a `minimumWidth`.
- Después, durante el ensamblaje, `ChemusonWindow._apply_theme()` aplica el stylesheet/QPalette al `QMainWindow`; el SidePanel queda oculto hasta una acción Ver.
- La hoja `QToolButton[sideTab="true"]` define padding, `font-size` y peso; `QScrollArea` se crea no-resizable y con scrollbars no visibles.
- `set_active()` vuelve a unpolish/polish cada botón y llama `ensureWidgetVisible()` para la pestaña activa. La prueba existente de cinco tabs mide después de invocar `set_active()` para cada una, por lo que no certifica la primera presentación tras show.
