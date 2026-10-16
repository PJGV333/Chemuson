# Validación — uniones de enlaces moleculares

## Alcance y estado

- Rama experimental: `feat/molecular-depiction-excellence`, basada en `origin/release/v0.3.0-beta.1-prep` (`9a70ee1c908a39e700098bcc0a54d464a606e24c`). El reporte OpenSpec/baseline se publicó antes de producción en `0856981`.
- No se modificaron prep, `main`, tags, releases, Clean2D, persistencia `.cmsn`, presets ni la semántica química. No se añadió `ChemUSON Publication`.
- La única modificación de producción centra las marcas hash en intervalos iguales; orden/orientación, número de marcas, grosores interpolados y `BondItem` API se conservan. El snapshot de prueba confirma que refrescar el item no modifica los átomos/enlace.
- `READY FOR MANUAL VISUAL REVIEW`; no es aceptación visual del propietario.

## Pruebas y verificaciones

| Verificación | Resultado |
|---|---|
| Baseline antes del código | 2238 IDs; ocho shards completados: 2218 passed, 20 skipped, 0 fallos. El F401 `math` indicado abajo ya existía. |
| `pytest --collect-only -q` final | 2267 tests collected; 29 IDs nuevos del test de geometría/galería. El updater añadió sólo esos IDs al manifiesto. |
| Focal de renderizado/geometría + wedge + doble + export + templates/Haworth + aromáticos | 74 passed. |
| `pytest -q tests/architecture tests/test_test_workflow_contract.py` | 284 passed. |
| Ocho shards (300 s externos) | Índices 0, 1, 2, 3, 5, 6 y 7 finalizaron: respectivamente 284/0, 284/0, 280/4, 278/5, 283/0, 274/9 y 282/1 passed/skipped. |
| Recuperación del shard 4 | La ejecución local simultánea agotó 300 s sin aserción fallida, durante `test_theme_refresh_rebuilds_rail_and_flyouts`. Se dividió ese assignment (sin repetir el shard combinado): 80 IDs Clean2D = 80 passed; `tests/test_ui_tool_rail.py` aislado = 64 passed; los 139 IDs restantes = 138 passed/1 skipped. Así los 2267 IDs quedaron cubiertos: 2247 passed, 20 skipped, 0 fallos. La herramienta `verify` no se declara verde porque el primer proceso del shard 4 fue timeout; el timeout y la recuperación quedan explícitos. |
| `python -m compileall src tests tools packaging` | Exit 0. |
| Ruff focalizado `F401,F811,F821,E722,E741` | Sin hallazgos en archivos cambiados. |
| Ruff global `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` | Exit 1 por el único F401 histórico: `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3:8`, `math` no usado. No se modificó. |
| `openspec validate improve-molecular-depiction-junctions --strict` | Cambio válido. |
| `git diff --check` | Exit 0. |

## Evaluación visual y límites

- Galería fuente: `tools/visual/render_bond_junction_gallery.py`; render Qt determinista en claro/oscuro, variantes de numeración y composición de la referencia. Las imágenes exploratorias están en `/tmp`.
- Evidencia comparativa versionada: `docs/visual/bond-junctions-before-after.png` (3360×2408, before/after por columnas, claro/oscuro por filas).
- Mejora confirmada: marcas hash uniformes y centradas; se reduce el saliente transversal que coincidía exactamente con el átomo compartido. Ver la lámina en el tamaño de uso/zoom real para decidir si el margen terminal resulta perceptible.
- Sin cambios deliberados a sigma/pi ni al contorno de cuña sólida. La supuesta deficiencia de perfiles sigma/pi no se confirmó; permanece como investigación posible, no como defecto corregido. El preset Publication se aplaza porque no se validó ese alcance ni sus contratos de preferencias/documento/exportación.
- La revisión visual del propietario sigue pendiente; no se presume aprobación por passing tests o por la evaluación del agente.
