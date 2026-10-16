# Validación — refinamiento de wedges y dobles enlaces

## Alcance implementado

- `canvas_structure.py`: detecta continuaciones dobles en extremos hash y contexto cíclico bold en el extremo ancho de cuñas.
- `items.py`: ajuste terminal limitado de la última marca hash, integración de cuña acotada en la unión cíclica, y centrado óptico por contorno de glifo cuando la etiqueta está junto a enlaces múltiples.
- `canvas_text.py`: centrado condicional, clearance de cuña sólida/hash medido desde la geometría visible y ancho efectivo de dobles.
- Se añadió cobertura geométrica y de invariantes en `tests/test_wedge_double_label_junctions.py`, se reforzaron límites en `tests/test_wedge_geometry.py` y se actualizaron los 10 node IDs nuevos del manifest CI.
- Galería de 15 casos en claro/oscuro con comparación idéntica de escala. Evidencia versionada: `docs/visual/wedge-double-labels-before-after.png` (SHA256 `19e51722ed7836a95b121dcff0dbef24622d19410ece5c9a14160f37e52b11f2`). Render after exploratorio en `/tmp/chemuson-depiction-refinement-after/` (dark `0a37b70ccdcecfb4781d0bdb3ad2c1109116511d25c129c8455c82335d45da7b`, light `010071f2fbf314962d84a8065b594483ecf5c48a8d3f73f5b858e1ee8cd7ec7f`).

Los cambios no alteran coordenadas ni topología: las pruebas comparan snapshots de identidad del grafo, enlaces y estéreo. No se modificaron Clean2D, ChemName, persistencia ni dependencias.

## Baseline

OpenSpec comenzó en `bd597dffc15688615aa7ceb3a519c920e499a5e6`; resultados exactos de baseline en `baseline.md`. `compileall` pasó, se recolectaron 2275 node IDs, Ruff selectivo reportó el único F401 histórico (`math`, `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`) y el pytest monolítico expiró a 480 s (~34%), sin aserción visible. El árbol de baseline estaba limpio y el upstream coincidía con ese SHA.

## Verificación del estado actual

| Verificación | Resultado |
|---|---|
| `pytest --collect-only -q` | **2285 tests collected** en 0.59 s; manifest actualizado con los 10 node IDs añadidos. |
| Pruebas focales (`test_wedge_double_label_junctions.py`, `test_wedge_geometry.py`, `test_canvas_drawing_edges.py`, `test_template_library.py`, `test_selection_tightness.py`) | **54 passed** en 1.32 s. |
| `pytest -q tests/architecture` | **280 passed** en 9.31 s. |
| `pytest -q tests/test_test_workflow_contract.py` | **4 passed** en 0.05 s. |
| `python -m compileall src tests tools packaging` | Exit 0. |
| Ruff selectivo sobre los archivos modificados (`F401,F811,F821,E722,E741`) | Exit 0. |
| Ruff global (`ruff check src tests tools packaging --select F401,F811,F821,E722,E741`) | Exit 1 por el mismo F401 histórico en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no se modificó. |
| `openspec validate refine-wedge-double-bond-junctions --strict` | Válido. |
| `git diff --check` | Sin errores. |

## Shards de la colección de 2285 tests

Plan de 8 shards generado en `/tmp/chemuson-depiction-refinement-final-plan.json`; distribución `[286, 286, 286, 286, 286, 285, 285, 285]`, SHA256 del manifiesto esperado `c35acb729dc75fc3813e04bb2b4412004e851b0a2c28da2551e404a39247bd`.

| Shard | Resultado |
|---|---|
| 1 | 286/286; exit 0, sin timeout (46.166 s). |
| 2 | 286/286; exit 0, sin timeout (22.527 s). |
| 3 | 280 passed, 6 skipped; exit 0, sin timeout (95.845 s). |
| 4 | **Bloqueado:** en dos ejecuciones completas falló `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend`; 282 passed, 3 skipped, 1 failed (≈203 s por ejecución). El log sólo registró `worker.finished`, con trabajo todavía activo y resultado pendiente; no fue timeout. |
| 5–8 | No ejecutados después del fallo inesperado del shard 4. |

Diagnóstico acotado: el archivo CompChem completo pasa aislado (4 passed/0.98 s) y el nodo aislado pasa (1 passed/0.18 s). Sin embargo, el nodo falla también tras los 166 tests que lo preceden en el orden del shard (166 passed, 1 failed/137.99 s). Es un fallo dependiente del orden/estado Qt; su causa y relación con este cambio no están establecidas. No se ocultó, omitió ni alteró el test. Los resultados de los shards 1–4 están en `/tmp/chemuson-depiction-refinement-final-artifacts/`; el segundo intento del shard 4, en `/tmp/chemuson-depiction-refinement-final-retry/`.

## Estado de aceptación

**NOT READY.** El shard 4 conserva un fallo reproducible bajo el orden del shard y los shards 5–8 no se ejecutaron. La suite no se declara verde. La comparación visual es automatizada; queda pendiente la revisión visual independiente del propietario. Sin merge a ramas protegidas, tag o release.
