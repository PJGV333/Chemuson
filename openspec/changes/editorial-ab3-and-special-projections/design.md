# Diseño — Preset AB3 y geometría de proyecciones

## Baseline y checkpoint

Rama `feat/molecular-depiction-visual-lab`, HEAD inicial `70efde4bca68adc6ffc55a10bd18727fb6fa1cec`, árbol limpio. Ese commit es el checkpoint de recuperación previo a cambios. La salida completa está en `baseline-command-output.txt`; resumen en `baseline.md`. Baseline: compileall PASS; 2,326 tests recogidos; `pytest -q` con 2,305 passed, 20 skipped y un fallo heredado/dependiente del orden en `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend`; Ruff focaliza sólo el F401 histórico `math` de `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py`.

## Preset y elegibilidad editorial

Un módulo de presets nombrados del laboratorio construye un `DepictionStyleParameters` inmutable desde el vector editorial aprobado; UI y CLI pueden activarlo explícitamente. `Current`, el constructor por defecto, `CHEMDOODLE_LIKE` y la inicialización normal del editor permanecen intactos. Las exportaciones identifican el preset y conservan los valores exactos. Los casos llevan metadatos explícitos de elegibilidad/rol: sólo A1, A2, A3, B1, B2 y B3 forman el conjunto de selección aprobado; A4 y B4 quedan excluidos, y B4 sigue en la galería y renderer como `geometric-fixture`, `non-editorial`, `excluded from preset selection / style calibration`. C1–C4 quedan como casos diagnósticos de proyección, no fuente de selección del preset. Los reportes/contact sheets enumeran estos roles para impedir contaminación de futuros resúmenes.

## Geometría de proyección

Crear `src/chemuson/gui/projections/` como código puro de geometría bajo el módulo catalogado M08, sin imports Qt, dialogs, canvas, items ni `DrawingStyle`. `common.py` define descripciones geométricas inmutables (átomos, enlaces y `ProjectionGeometry`) y conversión determinista a `MolGraph`. `haworth.py`, `fischer.py` y `chair.py` definen parámetros y constructores propios. La geometría puede etiquetar aristas frontales con `BondStyle.BOLD`, que es semántica local de la proyección; el grosor, cap/join y demás estilo editorial global siguen exclusivamente en el renderer.

Haworth usa vértices geométricos regulares, con orientación y tamaño explícitos, el O del anillo en posición canónica, política fija para el arco frontal y vectores de sustitución que respetan arriba/abajo y despejan el exterior del anillo. `HaworthParameters` expone tamaño/familia, longitud, aspecto/skew, arco frontal, anómero α/β, configuración D/L y sustituyentes opcionales por posición. El preset histórico beta-D-glucosa produce el mismo inventario básico de 11 átomos/11 enlaces, sin `BondStereo` sintético. `build_haworth_template(bond_length, anomeric_up, bold_front)` queda como adaptador compatible a la nueva API.

Fischer y silla dejan constructores utilizables y parámetros extensibles para backbone/cruces/terminaciones y conformación/vectores axial-ecuatorial, respectivamente; no asignan stereo metadata. `templates.py` mantiene los nombres públicos antiguos como adaptadores.

## Validación y riesgos

Añadir pruebas unitarias deterministas para el preset, defaults, selección editorial y reportes; pruebas geométricas de conteos, coordenadas, frontales, ausencia de stereo inventado y colisiones con el renderer real; pruebas de parámetros/constructores Fischer y silla; auditoría arquitectónica estática de imports puros. Actualizar el catálogo M08 para el nuevo subdirectorio y revisar los node IDs nuevos antes de actualizar el manifiesto CI. Generar evidencia Haworth antes/después con el mismo renderer.

El riesgo principal es el layout de etiquetas Qt, que depende del render efectivo y no sólo de coordenadas: se validará con `AtomItem`/`BondItem` y se exportará evidencia. La suite monolítica basal tarda más de 20 minutos y tiene un fallo conocido; la verificación final usará comandos acotados, cubrirá colección completa y reconciliará shards sin ocultar el fallo.
