# Baseline — 2026-09-25-modernize-ui-tool-rail-flyouts

Capturada el 2026-09-25 en la rama `ui/modernization` (HEAD `d5405d9`),
**antes** de tocar la implementación de la Fase 4.

## Entorno de validación

`.venv` del repositorio (CPython 3.12.12, creado con `uv`; `chemuson` editable).
Dependencias de validación: `pytest==9.1.1`, `ruff==0.16.9`, `PyYAML`. Sin
dependencias de producción nuevas.

```
.venv/bin/python -m pytest -q
.venv/bin/python -m ruff check src tests tools packaging --select F401,F811,F821,E722,E741
.venv/bin/python -m compileall src tests tools packaging
```

## git status --short

```
(clean)
```

## python -m compileall src tests tools packaging

```
compileall rc=0
```

## pytest --collect-only -q

```
1701 tests collected in 0.43s
```

## pytest -q

```
4 failed, 1677 passed, 20 skipped in 252.87s (0:04:12)
```

Fallos (conocidos, previos a la Fase 4, no introducidos por el venv; ver
`AGENT_REPORT.md` del setup del entorno — son de código, reproducibles con
rdkit 2023.09.2 y 2026.3.6):

- `tests/test_clean2d_engine_candidates.py::test_generate_candidates_attempts_rdkit_for_cyclic_graphs`
- `tests/test_smiles_stereo_import.py::test_chiral_smiles_import_creates_wedge_or_hash`
- `tests/test_smiles_stereo_import.py::test_amino_acid_chiral_smiles_import_creates_wedge_or_hash`
- `tests/test_smiles_stereo_import.py::test_tetrandrine_import_preserves_visual_stereo`

## ruff check src tests tools packaging --select F401,F811,F821,E722,E741

```
F401 [*] `math` imported but unused
 --> tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3:8
Found 1 error. (rc=1)
```

(Aviso preexistente; fuera del alcance de esta fase: no se toca por la regla
de no-refactor oportunista.)

## Inventario 1:1 de tool_ids / acciones / señales (auditoría Fase 4)

### `ChemusonToolbar` (izquierda, `toolbar.py:61`)

Señales: `tool_changed(str)`, `bond_palette_changed(object)`,
`ring_palette_changed(object)`, `element_palette_changed(str)`,
`periodic_table_requested()`. `action_group` (QActionGroup exclusiva,
compartida con `SymbolPaletteToolbar`).

| Botón (QAction objectName) | Tipo | tool_ids de la paleta |
|---|---|---|
| `tool_select` (select_button) | paleta (2, cols=2) | `tool_select`, `tool_select_lasso` |
| `tool_bond` (bond_button) | paleta (11, cols=3) + `bond_palette_changed` | `tool_bond` (specs: sencillo, bold, doble, triple, aromático, wedge, hashed, wavy, flexible, interacción, coordinativo) |
| `tool_chain` (chain_action) | acción simple | `tool_chain` |
| `tool_ring` (ring_button) | paleta (11, cols=4) + `ring_palette_changed` + footer "Tamaño personalizado..." (QInputDialog) | `tool_ring` (benceno + anillos 3–12) |
| `tool_atom` (label_button) | paleta (10, cols=5) + `element_palette_changed` + footer "Tabla periódica..." | `tool_atom` (C, N, O, S, P, F, Cl, Br, I, H) |
| `tool_coordination_center` (coord_action) | acción simple | `tool_coordination_center` |
| `tool_rotate_3d_precise` (rotate_3d_precise_action) | acción simple | `tool_rotate_3d_precise` |

### `SymbolPaletteToolbar` (derecha, `toolbar.py:829`)

Señales: `tool_changed(str)`, `atomic_diagram_requested()`,
`diatomic_mo_diagram_requested()`, `ligand_field_diagram_requested()`,
`electronic_diagram_preset_requested(str)`.

| Botón (QAction objectName) | Tipo | tool_ids de la paleta |
|---|---|---|
| `tool_text` (text_button) | paleta (menú de acciones de formato + colores, vía `set_text_menu`) | `tool_text` |
| `tool_brackets` (bracket_button) | paleta (10, cols=2) | `tool_brackets_square`, `_square_left`, `_square_right`, `_corner`, `_curly`, `_curly_left`, `_curly_right`, `_frame`, `_frame_rounded`, `_round` |
| `tool_annotation` (annotation_button) | paleta (16, cols=4) | `tool_arrow_line`, `_line_dashed`, `_forward`, `_forward_open`, `_forward_dashed`, `_retro`, `_retro_open`, `_retro_dashed`, `_both`, `_both_open`, `_both_dashed`, `_equilibrium`, `_equilibrium_dashed`, `_retrosynthetic`, `_curved`, `_curved_fishhook` |
| `tool_plates` (plate_button) | paleta (2, cols=2) | `tool_tlc`, `tool_electrophoresis` |
| `tool_symbol_palette` (symbol_button) | paleta (12, cols=4) | `tool_charge_plus`, `tool_charge_minus`, `tool_charge`, `tool_symbol_plus`, `tool_symbol_minus`, `tool_symbol_radical`, `tool_symbol_lone_pair`, `tool_symbol_wavy_anchor`, `tool_symbol_radical_cation`, `tool_symbol_radical_anion`, `tool_symbol_partial_plus`, `tool_symbol_partial_minus` |
| `tool_energy_diagrams` (energy_diagram_button) | paleta (8, cols=2) + submenús (3 diálogos + 25 presets) | `tool_energy_diagram_{sublevel_s, custom_level, sublevel_p, sublevel_d, sublevel_f, hybrid_sp, hybrid_sp2, hybrid_sp3}` |
| `tool_orbitals` (orbital_button) | paleta (23 celdas, grid 4×7) | `tool_orbital_{kind}` para los 23 kinds de `ORBITAL_PALETTE_MODEL` |

### Consumidores (ventana/canvas)

- `main_window.py`: `toolbar.tool_changed → _on_tool_changed` (+ `_update_status`);
  `bond/ring/element_palette_changed → _handle_*_palette` (→
  `canvas.set_active_bond/ring/element`); `periodic_table_requested →
  _show_periodic_table`; `symbols_toolbar.tool_changed → _on_tool_changed`
  (+ `_update_status`); 4 señales de diagramas electrónicos → diálogos/presets.
- `_clear_active_tool_selection()` (cambio de pestaña): `tool_none` +
  `toolbar.clear_tool_selection()`.
- `canvas/canvas_selection_input.py::set_current_tool(tool_id)`: normaliza
  `atom_*`/`coord_*`/`bond_*`/`tool_brackets_*`/… y actualiza
  `state.active_tool`.
- `main_window_ui_builder.tool_status_label` mapea los mismos ids a texto.

### Atajos: auditoría de conflictos (letra simple)

Atajos existentes sin modificadores: **ninguno** (todos los `QKeySequence`
existentes llevan Ctrl/Shift/Alt o son F-keys; el canvas solo usa `A`/`E`
con Ctrl/Meta en `_handle_select_all`). Por lo tanto `V, A, L, B, R, C, T,
N, G, E, O` no colisionan con nada existente; se habilitan los 9 propuestos
en la tarea + `A` (lasso) y `L` (cadena) del mockup.

## Decisión `orbitals.py`

`draw_orbital_icon` usa `QPainter`/`QPainterPath`/`QPixmap` (renderer propio,
gradientes y caminos complejos). Migrarlo a la infraestructura SVG de la
Fase 2 implicaría redibujar 23 orbitales con riesgo de cambiar el aspecto
exacto y sobrecargar la fase. **Decisión**: el flyout reutiliza
`draw_orbital_icon(kind)` (QPainter) tal cual; se documenta como residuo
para una fase posterior. Sin cambio de comportamiento.
