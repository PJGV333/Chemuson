# Validation — fix-3d-dock-tools-and-conformers

**Fecha:** 2026-10-17 (MST, UTC−07:00). **Rama:** `feature/chemuson-3d`.

## Diagnóstico verificado

1. El estado inicial era contradictorio: `MolecularViewportWidget.edit_mode` estaba apagado mientras el botón Move aparecía marcado. La prueba de regresión fallaba por este estado.
2. Ante un resultado sin coordenadas, `_on_compchem_job_finished()` dejaba el detalle solo en `CompChemDock`; el viewport se había vaciado al comenzar el trabajo y no recibía un motivo ni un reintento. Por eso un fallo del worker podía parecer simplemente una vista vacía. Se corregirá en el hito C.
3. La generación aislada sí pudo generar la tetrandrina de referencia, así que “ETKDG no puede generar tetrandrina” no es una explicación general. No se encontró en el repo el `.cmsn` del dibujo del usuario; las pruebas locales existentes solo son “tetrandrine-like”. La causa del archivo real sigue pendiente hasta probar ese documento.

Referencia exacta utilizada para una regresión reproducible (no equivale al documento 2D original del usuario): PubChem CID 73078, fórmula C38H42N2O6, isomeric SMILES:

```text
CN1CCC2=CC(=C3C=C2[C@@H]1CC4=CC=C(C=C4)OC5=C(C=CC(=C5)C[C@H]6C7=C(O3)C(=C(C=C7CCN6C)OC)OC)OC)OC
```

El baseline produjo 46 átomos pesados, 52 enlaces, 46 coordenadas documentales, dos centros `S` preservados, método `mmff`, energía 243.577994 kcal/mol y sin error en aproximadamente 0.327 s. PubChem identifica el CID, la fórmula y la cadena isomérica usada ([CID 73078](https://pubchem.ncbi.nlm.nih.gov/compound/tetrandrine)).

## Hito A

- Se añadió tema claro/oscuro scoped a `MolecularViewportDock`, paleta desde tokens, estilos explícitos para estados y re-aplicación en `topLevelChanged`, `dockLocationChanged` y `_apply_theme()`.
- El dock inicia en navegación sin botón de edición marcado; seleccionar herramienta activa modo edición; desactivarlo limpia el check visual sin perder la preferencia de herramienta.
- Tests incluidos: tema claro/oscuro, paleta, señal de flotado/acoplado, preservación de escena/cámara/selección/tool; tema principal llega al dock; contexto renderer retiene escena al reinicializar.
- `pytest -q tests/test_viewport3d_widget.py tests/test_compchem3d_dock.py tests/test_ui_theme_foundation.py` → **80 passed, 1 skipped in 76.77s**.
- `pytest -q tests/architecture` → **280 passed in 10.42s**.
- `compileall` → exit 0; Ruff scoped → All checks passed; OpenSpec strict → válido; `git diff --check` → limpio.
- **No es validación visual GPU/compositor:** no se inició la GUI; la comprobación de legibilidad en KDE y recuperación del contexto real siguen manuales.

## Hito B

- La selección múltiple mantiene orden de pick; la malla resalta átomos, H derivados pickeados y enlaces, y Ctrl selecciona el componente conectado.
- El cambio de elemento usa `ChangeAtomCommand` con validación de valencia. Construcción/enlace y cambio de orden se validan contra un grafo candidato y se confirman como una sola acción undoable.
- Movimiento incluye átomos/fragmentos, H derivados y profundidad con Alt+arrastre o botón central; la ayuda contextual documenta los gestos.
- El dock ofrece handles guiados de longitud, ángulo y dihedro. Las pruebas verifican mediciones en Å/grados, una transacción undoable, Escape, conservación de H derivado y rechazo de enlaces cíclicos/geometría degenerada.
- La prueba parametrizada de interacción de viewport cubre grafos C–C (etano), C–C–O (etanol) y C–C–C–C (butano), sin entrada manual de IDs.
- `pytest -q tests/test_viewport3d_camera.py tests/test_viewport3d_geometry.py tests/test_viewport3d_widget.py tests/test_compchem3d_dock.py tests/test_geometry3d_editing.py tests/test_geometry3d_editing_metrics.py` → **102 passed, 1 skipped in 37.73s**.
- Ruff scoped sobre archivos B → **All checks passed**; `git diff --check` → limpio.
- **No es aceptación visual:** no se inició GUI; contraste real KDE/OpenGL/compositor permanece pendiente.

## Hito C — implementación y validación sin aceptación visual

- Se mantuvo el contrato `metadata["stereochemistry"] = {"atoms": ..., "bonds": ...}` y se separaron centros potenciales/no asignados en `stereochemistry_diagnostics`; esto conserva compatibilidad con consumidores y pruebas existentes.
- El worker RDKit aislado registra transiciones reales, tiempos por etapa, conteos/estado químico (carga formal, valencia, aromaticidad y anillos), versión/código de retorno RDKit, semillas/intentos, estado de fuerza de campo y estereoquímica especificada/no asignada. ETKDG prueba hasta cuatro combinaciones finitas con semillas distintas y parámetros macro cíclicos solo cuando hay anillos de tamaño ≥8; valida coordenadas y estereoquímica antes y después de minimizar, y preserva el embedding válido como `unoptimized` ante fallo de minimización.
- La ruta NDJSON propaga progreso, timeout y cancelación con terminación del proceso aislado; servicio/backend/controlador reenvían diagnóstico y progreso. La UI muestra progreso, error/reintento y descarta resultados obsoletos con código `stale_result`.
- La escena expone estado `scene_ready`, upload y render solicitados; fallos de construcción o display se identifican como `scene_failed`/`display_failed` y llegan al viewport y al panel CompChem. `render_requested` confirma únicamente que el renderer recibió la solicitud; **no es aceptación visual GPU**.
- La regresión de referencia PubChem CID 73078 verifica 46 átomos pesados, 52 enlaces, coordenadas finitas y dos centros S, además de agua, etano, etanol, butano, benceno, ciclohexano, eteno, etino y adamantano. El SMILES de referencia no sustituye el `.cmsn` original del usuario.
- `pytest -q tests/test_rdkit_live_stream.py tests/test_geometry3d_service.py tests/test_geometry3d_stereochemistry.py tests/test_compchem3d_dock.py tests/test_viewport3d_camera.py tests/test_viewport3d_geometry.py tests/test_viewport3d_widget.py tests/test_geometry3d_editing.py tests/test_geometry3d_editing_metrics.py` → **166 passed, 2 skipped in 80.18s**.
- `pytest -q tests/architecture` → **280 passed in 12.54s**; `compileall` → exit 0; Ruff scoped en archivos del cambio (`F401,F811,F821,E722,E741`) → **All checks passed**; OpenSpec strict → válido; `git diff --check` → limpio.
- Ruff de repositorio con las reglas baseline sigue con el mismo F401 preexistente y fuera de alcance: `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3` (`import math`). No se repitió la suite completa de pytest por el SIGSEGV Qt/SIP ya observado en baseline.

## Pendiente

No se inició GUI: la aceptación KDE/OpenGL/compositor permanece manual. Tampoco se encontró/probó el `.cmsn` original de tetrandrina; la referencia CID 73078 es una regresión separada y no se atribuye al documento real.
