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

## Pendiente

Los hitos B y C, sus regresiones y commits están pendientes. El documento original tetrandrina no está disponible localmente. No se atribuye el resultado exitoso de PubChem al dibujo 2D real.
