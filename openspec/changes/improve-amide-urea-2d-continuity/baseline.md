# Baseline — continuidad 2D de amidas y ureas

## Estado y comandos obligatorios

La rama es `fix/initial-2d-geometry-beta1`, HEAD previo `ef241d4a06434fe719a1d36752138eb2806cd53c` (trabajo anterior conservado, 3 commits locales sobre `3e0ef606`). No hay cambios tracked ni staged. `git status --short` sólo muestra el prompt anterior `?? ChemUSON_Campana_Geometria_Inicial_2D_Beta1.md`, que permanece sin seguimiento y fuera de commits.

| Comando | Baseline |
|---|---|
| `python -m compileall src tests tools packaging` | Exit 0. |
| `pytest --collect-only -q` | 2,578 tests collected en 0.97 s. |
| `timeout --signal=TERM 600s pytest -q` | Exit 124 al límite de 600 s, progreso ~30%; ya mostraba un fallo y skips, sin resumen de pytest porque fue terminado. |
| `timeout --signal=TERM 180s pytest -x -q` | 366 passed, 1 failed en 35.10 s: `tests/test_beta_compact_shell_and_tabs.py::test_view_menu_submenus_reuse_actions_and_remove_obsolete_auxiliary_option`; el menú incluye `Visor molecular 3D` y stderr indica `QOpenGLWidget is not supported on this platform`. Es el mismo fallo basal documentado en la campaña anterior. |
| `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` | Exit 1 por el F401 histórico `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no se modifica. |
| `pytest -q tests/test_functional_group_hydrogen_adjustment.py tests/test_sp3_bond_placement.py tests/test_canvas_drawing_edges.py tests/test_wedge_geometry.py tests/test_wedge_double_label_junctions.py tests/test_persistence_contract.py tests/test_depiction_lab_renderer.py tests/test_depiction_lab_parameters.py` | 90 passed en 2.18 s. |

Logs exactos: `baseline/01-git-status.log` a `baseline/05-ruff.log`, `baseline/06-pytest-first-failure.log.gz` (gzip del stdout exacto; conserva los espacios finales de traceback), y `baseline/07`–`09`. El JSON/QImage before proviene de eventos `QTest` y está en `evidence/before/`.

## Reproducción de `Ar–NH–C(=O)–NH–CH3`

El benceno usa press/move/release de la herramienta de anillo; los enlaces `Ar–N`, `N–C(=O)`, `C–N` y `N–CH3` se construyen con clic corto; `C=O` y `N–CH3` también se contrastan con arrastres dirigidos. Cada clic guarda el endpoint de preview previo al release. Los puntos generados coinciden con el commit. El caso contiene el anillo (6 átomos/enlaces) y la ruta N/C/O (5 átomos, 5 enlaces); no hay corrección de Clean2D.

| Enlace | Distancia centro-a-centro | Longitud configurada | Trim inicio/fin | Longitud visible del path principal |
|---|---:|---:|---:|---:|
| Ar–N1H | 40.000000000000000 | 40.0 | 0 / 5.631875 | 34.368125 |
| N1–C(=O) | 39.999999999999986 | 40.0 | 15.620745 / 0 | 24.379255 |
| C=O (2 trazos) | 39.999999999999964 | 40.0 | 0 / 8.204943 | 31.795057 por trazo |
| C(=O)–N2H | 39.999999999999986 | 40.0 | 0 / 15.743313 | 24.256687 |
| N2H–CH3 | 39.999999999999986 | 40.0 | 5.679473 / 0 | 34.320527 |

La separación de centros nunca se reduce para resolver el defecto. Los trims de los enlaces no aromáticos de la referencia se deben al label-shrink; trims de endpoint son cero. Los ángulos en clic-only son 120°; con drag el rango medido es 119.485807°–120.514193°.

### Glifos y líneas de enlace

- N1: centro `(479.000, 499.000)`; label `NH` 11 pt; N bbox escena `(476.9025,497.4637,8.375,10.3125)`; H `(487.7306,497.4637,8.375,10.3125)`.
- N2: centro `(479.000,429.718)`; N bbox `(476.8897,420.9344,8.375,10.3125)`; H `(487.7179,420.9344,8.375,10.3125)`.
- Anchaje óptico actual: el centro del glyph N queda en el centro del átomo más el offset del grupo; el H inline se ordena al lado de N. El reporte JSON conserva posiciones de label/glyph completos.
- N1→C: proyección actual `13.320745` + clearance `2.30` = trim `15.620745`; pero al cruzar el `shape()` real con el corredor de stroke 2.6 + 2×clearance 2.3 la extensión de intersección termina en `7.623205` (exceso `7.997540`).
- C→N2: proyección `13.443313` + `2.30` = `15.743313`; corredor real hasta `7.795145` (exceso `7.948168`).

La causa se clasifica como **recorte visual excesivo por proyección sin intersección**, no error de coordenadas. El N óptico y el H visible están desplazados localmente; no deben compensarse moviendo centros ni acortando la referencia de 40. La orientación del caso mínimo ya es trigonal, pero falta comprobar/encapsular la clasificación química para los casos de amidas/ureas y contrastarla con aminas/aromáticos.
