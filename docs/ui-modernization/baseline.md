# Baseline — campaña de modernización de UI

**Fecha:** 2026-09-22
**Checkout de referencia:** `346b22924b0501cf2239d36a66c4c1ce7fc2fed7` ("Archive Clean2D Campaign 3 safety closure")
**Rama de trabajo:** `ui/modernization` (worktree `/home/unison-pjgv/Documentos/GitHub/Chemuson-ui-modernization`)
**Alcance de este trabajo:** solo `docs/ui-modernization/` (documento + maqueta). **No se toca código de producción.**

> Nota de aislamiento: la campaña Clean2D continúa en el checkout principal
> (`clean2d/campaign-implementation`). Este worktree se creó con
> `git worktree add -b ui/modernization … 346b229` sin cambiar de rama el
> checkout principal. Los cambios sin commit del agente anterior
> (`docs/ui-modernization/PLAN.md`, `docs/ui-modernization/mockup-ui.html`)
> se trasladaron íntegros y exclusivamente a este worktree.

## Comandos de baseline (AGENTS.md §1.2)

### `git status --short --branch` (checkout principal, antes de mover `docs/ui-modernization/`)

```
## clean2d/campaign-implementation...origin/clean2d/campaign-implementation
?? docs/ui-modernization/
```

(Tras trasladar el directorio no versionado a este worktree, el checkout
principal queda con status vacío — sin cambios propios de Clean2D pendientes.)

### `git worktree list`

```
/home/unison-pjgv/Documentos/GitHub/Chemuson                  346b229 [clean2d/campaign-implementation]
/home/unison-pjgv/Documentos/GitHub/Chemuson-ui-modernization 346b229 [ui/modernization]
```

### `python -m compileall src tests tools packaging`

```
(exit 0, sin salida — compila limpio)
```

### `pytest --collect-only -q` / `pytest -q`

**No ejecutables en este entorno:** `pytest` no está instalado ni en
`.venv/bin/python` ni en el Python del sistema (`No module named pytest`).
Se documenta como limitación del entorno; el trabajo de esta etapa es solo
documentación (`.md`/`.html`), por lo que no afecta la suite.

### `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`

**No ejecutable en este entorno:** `ruff` no está instalado
(`No module named ruff`). Misma justificación que pytest.

## Verificación de integridad del plan (auditoría de PLAN.md vs. repo)

Comprobado contra `346b229`:

- Conteo de líneas §1.2: `main_window.py` 2164, `main_window_ui_builder.py` 536,
  `toolbar.py` 1517, `docks.py` 940, `styles.py` 773, `icons.py` 1320,
  `shell/assembly.py` 271, `tab_manager.py` 265 → **todo exacto**.
- `PyQt6` presente en `requirements.txt` (línea 1) y `pyproject.toml`.
- 7 docks (`PlantillasDock`, `InspectorDock`, `ChemicalPropertiesDock`,
  `SpectroscopyDock`, `CompChemDock`, `ValidationDock`,
  `AppearanceDock`) en `docks.py`; los 7 ocultos por defecto en `shell/assembly.py`.
- Menús: Archivo, Editar, Ver, Estructura, Reacción (placeholder "Próximamente"), Ayuda.
- Señal `tool_changed = pyqtSignal(str)` en `toolbar.py` (líneas 67 y 834).
- `architecture/modules.yml`: M08 = GUI, M09 = canvas, M21 = `platform.settings` ✓.
- APIs citadas existen: `styles.get_main_stylesheet(theme_name)`,
  `icons.draw_generic_icon/draw_atom_icon/…`, `template_browser_service.py`.
- **Corrección aplicada (supuesto obsoleto):** §1.1 afirmaba "No hay ningún
  asset gráfico en el repo (ni SVG ni PNG)". El repo contiene PNG de baselines
  de regresión (`assets/baseline/orbitals_palette.png`, `src/repro_v2.png`,
  `tests/archive/…`) pero **no** assets SVG ni assets de iconos de UI. La frase
  se precisó en consecuencia.
