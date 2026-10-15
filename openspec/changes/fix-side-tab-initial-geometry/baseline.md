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

## Baseline — bloqueo de selección en CI (HEAD `29e3594618afd131ec8f22999d664f6e43185c69`)

Antes del cambio de prueba, el árbol estaba limpio en `fix/beta-properties-panel-polish`; el remoto apuntaba al mismo HEAD.

- `python -m compileall -q src tests tools packaging`: PASS.
- `pytest --collect-only -q`: 2154 node IDs.
- Ruff en `side_panel.py` y `test_spectroscopy_dock.py`, reglas de errores/imports seleccionadas: PASS.
- El test aislado pasó en Python 3.14.7 / PyQt 6.11.0; el shard 1 fallaba en CI Ubuntu 24.04 / Python 3.11.17 / pytest 9.1.1.
- Se creó `/tmp/chemuson-ci-py311` (fuera del repo) instalando `requirements.txt`, `requirements-dev.txt` y el editable del proyecto; resolvió PyQt6 6.11.0 / Qt 6.11.2, pytest 9.1.1, NumPy 2.4.6 y RDKit 2026.9.1. La colección reprodujo el plan CI exacto: 2154 IDs, 8 shards `[270,270,269,269,269,269,269,269]`, SHA `516534cabacbc81975501ab9ea0e5a7420e95affc070d15a04d61be68f446a93`.
- En ese entorno equivalente, el test aislado y el archivo completo pasan; el shard 1 sin cambios falla reproduciblemente: 267 passed, 2 skipped, 1 failed en 89.24 s (`isSelected()` falso). El par de tests `test_ctrl_k_shortcut_triggers_clean_2d_full` y `test_ctrl_alt_k_shortcut_triggers_clean_2d_propose_action`, seguido del test fallido, reproduce el caso en 1.18 s.
- Instrumentación temporal, eliminada antes del fix: al inicio del tema light, `currentRow/Column=-1/-1`, modifier `ControlModifier`; `setCurrentCell` selecciona. Para dark, `currentRow/Column=0/1` y la celda ya está seleccionada, pero persiste `ControlModifier`; repetir `setCurrentCell` conmuta la selección a False. El estado compartido procede de los eventos sintéticos Ctrl de tests previos. No se toca producción.
- Con liberación explícita de modificadores al inicio del test y aserción de `NoModifier`, la secuencia de tres tests pasa y el shard 1 completo pasa: 268 passed, 2 skipped en 90.84 s.
