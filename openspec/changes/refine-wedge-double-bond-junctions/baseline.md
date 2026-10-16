# Baseline — refinamiento de wedges y dobles enlaces

## Git y entorno

- Rama: `feat/molecular-depiction-excellence`
- HEAD inicial: `bd597dffc15688615aa7ceb3a519c920e499a5e6` (`fix(gui): improve canvas edge and Haworth drawing`)
- Upstream: `origin/feat/molecular-depiction-excellence` en el mismo SHA; árbol limpio (`git status --short` sin salida).
- Qt ejecutado con `QT_QPA_PLATFORM=offscreen` para renders/pruebas.

## Comandos de baseline

| Comando | Resultado |
|---|---|
| `python -m compileall src tests tools packaging` | Exit 0; todos los árboles compilaron. |
| `pytest --collect-only -q` | Exit 0; 2275 node IDs collected in 0.57s. |
| `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` | Exit 1 por el hallazgo histórico único `F401 math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no corregido. |
| `QT_QPA_PLATFORM=offscreen timeout 480s pytest -q` | Exit 124 (timeout) a 480s; progreso ≈34%, sin fallo de aserción visible. Log exacto: `/tmp/chemuson-depiction-refinement-baseline/pytest-baseline.log`. |
| Render existente `tools/visual/render_bond_junction_gallery.py --theme both` | Exit 0. PNG baseline en `/tmp/chemuson-depiction-refinement-baseline/before-existing/`; SHA256 dark `b6f84ef8ae69a91c032051b9912db215a7408a94d95b5fa8b7b5db22dc948889`, light `3c65230cc8b149a2602b2ea132162e8e89976e2cae626c6eadd198cc738c7e83`. |

## Primera lectura técnica (hipótesis, no conclusiones)

- La galería base ya contiene casos simples/bold/dobles/cuñas y 20 casos de clearance, pero no los casos dedicados de etiqueta O en C=C ni Haworth pentágono/hexágono wedge↔bold solicitados en esta iteración.
- `canvas_structure._set_bond_item_join_context()` retira contexto si `bond.style != WEDGE`; por tanto las marcas HASHED no adaptan su endpoint con el doble vecino. `BondItem.update_positions()` genera hashes transversales sin clipping de vecinos.
- `canvas_text._update_bond_label_shrinks()` mantiene clearance fijo de 1.5 para WEDGE y un valor ligado al half-width para HASHED. Aún hay que medir caps/paths en casos concretos antes de cambiar la política.
- `_label_open_direction()` cuenta direcciones de enlace, no orden/anchura; `AtomItem._center_label()` usa advance/bounding box. Aún hay que separar el efecto óptico de font metrics de un posible conflicto con sigma/pi.
- Los joins sólidos ya reciben vecinos plain/bold y usan miter/fork; hay que aislar el caso Haworth antes de alterar parámetros.
