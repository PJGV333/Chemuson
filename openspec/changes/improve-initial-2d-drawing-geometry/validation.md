# Validación — geometría inicial 2D

## Registro antes/después

| Caso | Ángulo previsto | Antes | Después | Longitud prevista/obtenida | Cruces/excepción |
|---|---:|---:|---:|---|---|
| Benceno aislado | 120° | 120° (máx. error medido `1.2e-13°`) | 120° (máx. error `6e-14°`) | 40.0 / `39.99999999999997–40.00000000000003` | 0; renderer sin cambios |
| Benceno con sustituyente, vértice girado 17° | 120° local y radial | 137°/103°, radial error 17°; el clic corto no creaba el enlace | 120°/120° y radial error < `1e-6°`; el clic corto añade átomo/enlace | 40.0 / 40.0 | 0; los 6 átomos y enlaces aromáticos previos quedan idénticos |
| Sugerencia sp3 y cadena acíclica | 120° | 109.5° | 120° | 40.0 / `39.99999999999997–40.00000000000003` | 0 |
| Centro carbonilo | 120° trigonal | 119.69–120.32° con eventos dirigidos cuantizados | 120° en las tres separaciones locales | 40.0 / 40.0 | 0 |
| 1-fenilbutan-2-ona construida por clics cortos | C10H12O y geometría 2D coherente | 11 átomos/11 enlaces al arrastrar; local 119.69–120.32° | 11 átomos/11 enlaces; SMILES `CCC(=O)Cc1ccccc1`; ángulos 120° | 40.0 / `39.99999999999997–40.00000000000003` | 0 cruces no adyacentes |
| Anillo anclado a enlace C–O / anclado a átomo | hexágono regular | No había regresión reproducida | 6 lados de 40.0, 120°; la arista compartida se reutiliza y los átomos existentes no se mueven | 40.0 / 40.0 | 0; O se conserva |

Los cuatro ángulos internos del caso completo son 120° dentro de `6e-14°`; los trazos/aromaticidad siguen siendo los de producción. La constante tetraédrica de compatibilidad permanece en 109.5°; sólo las nuevas sugerencias interactivas sp3 suministran 120°.

## Baseline

Ver `baseline.md`, logs completos en `baseline/` y capturas/modelo en `evidence/before/`. Resultados importantes:

- Colección basal: **2,567 tests**.
- `pytest -q` basal: excedió el límite externo de **600 s**.
- `pytest -x` basal: **366 passed** antes de fallar en `tests/test_beta_compact_shell_and_tabs.py::test_view_menu_submenus_reuse_actions_and_remove_obsolete_auxiliary_option`; el menú incluye el visor 3D.
- Ruff basal: un único F401 preexistente (`math`) en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.

El nodo del menú se repitió aislado tras el cambio y mantiene el mismo fallo de visor 3D (`QOpenGLWidget is not supported on this platform`); no se modificó ese subsistema ni se ocultó/actualizó la aserción.

## Validación posterior

- `timeout 120s pytest -q tests/test_initial_2d_drawing_geometry.py`: **11 passed** (incluye polígonos n=3–6, zoom/DPR, anclaje a átomo/enlace con O, colisiones, retícula de 30°, `Alt`/`Shift`, preview/commit, undo/redo, identidad química y round-trip `.cmsn`).
- Familias focales 2D, Clean2D sin cambios, persistencia, Visual Lab/proyecciones, CI node IDs y arquitectura: **469 passed, 2 skipped in 40.10s**. Incluye `tests/architecture` y `tests/test_test_workflow_contract.py`.
- `pytest --collect-only -q`: **2,578 IDs**, iguales al manifiesto `tests/ci/expected_pytest_nodeids.txt` (sin faltantes, extras ni duplicados; comparación normalizada por orden lexicográfico).
- `python -m compileall -q src tests tools packaging`: exit 0.
- Ruff scoped en los cinco módulos de producción cambiados, la prueba y el generador posterior: **All checks passed** (`F401,F811,F821,E722,E741`). El comando global conserva únicamente el F401 basal citado arriba.
- `openspec validate improve-initial-2d-drawing-geometry --strict --no-interactive`: **valid**.
- `git diff --check`: exit 0.

La suite monolítica no se repitió: su baseline ya documenta timeout y el fallo previo aislado sigue reproduciéndose. Los grupos acotados posteriores relevantes sí pasan; no se declara la suite completa verde.

## Evidencia posterior reproducible

Comando, desde la raíz del repo:

```bash
PYTHONPATH=src QT_QPA_PLATFORM=offscreen \
INITIAL_2D_OUT=openspec/changes/improve-initial-2d-drawing-geometry/evidence/after \
python openspec/changes/improve-initial-2d-drawing-geometry/evidence/reproduce_initial_2d_after.py
```

El script crea un anillo girado, construye la secuencia con clics cortos reales de Qt, guarda 12 imágenes del renderer actual (1×/1.25×, fondos de workspace claro/oscuro) y verifica la identidad RDKit, coordenadas y cruces. Resulta: SMILES `CCC(=O)Cc1ccccc1`, 11 átomos, 11 enlaces, 0 cruces no adyacentes, lados entre `39.99999999999997` y `40.00000000000003`, ángulos del anillo entre `119.99999999999994°` y `120.00000000000006°`. Métricas y modelo están junto a las imágenes en `evidence/after/`; `SHA256SUMS` registra los 14 PNG/JSON/TXT y `sha256sum -c SHA256SUMS` los valida.

Las variantes “claro/oscuro” son fondos de escena del renderer (`#F1F5F9` / `#111827`); el papel de ChemUSON sigue blanco en ambos. Son comprobaciones Qt/QImage offscreen, no una captura del escritorio ni aceptación visual manual. El propietario debe hacer esa revisión en CachyOS antes de integrar.

## Revisión manual pendiente en CachyOS

1. Abrir una hoja vacía, crear benceno girado, y hacer clic corto sobre un vértice: debe aparecer un enlace radial de 40 unidades sin cambiar ningún átomo/enlace del anillo.
2. Construir `Ph–CH2–C(=O)–CH2–CH3` con clics cortos y orden 2 para C=O; comprobar que la vista previa coincide con el commit y que cada longitud es uniforme.
3. Añadir un anillo a un enlace C–O y otro a un átomo ya sustituido; comprobar cierre, fusión, heteroátomo y extremos compartidos inmóviles.
4. Probar arrastres dirigidos, ángulos/pasos/longitudes fijos, `Alt`, `Shift`, aromaticidad, undo/redo y guardar/abrir `.cmsn`.
5. Revisar capturas del caso a 1×/1.25× en los temas claro/oscuro disponibles y confirmar que se conservan los refinamientos visuales beta.

No se hicieron push, merge, PR, tags ni publicación. El prompt de campaña continúa sin seguimiento y fuera de los commits.
