# Validation — fix-3d-dock-tools-and-conformers

**Fecha:** 2026-10-16 (MST, UTC−07:00). **Rama:** `feature/chemuson-3d`.

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

## Pendiente

El Hito C, su commit local y regresiones de generación/diagnóstico están pendientes. No se encontró el `.cmsn` original de tetrandrina dibujado por el usuario; la referencia PubChem es una regresión separada y no se atribuye al documento real. La aceptación gráfica manual del Hito A también sigue pendiente por las restricciones de esta campaña.
