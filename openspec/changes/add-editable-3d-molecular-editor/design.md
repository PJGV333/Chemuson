# Design — editor molecular 3D y mecánica molecular en vivo

## Estado actual comprobado

- `MolGraph` es la única autoridad química y ya se comparte con `ChemusonCanvas`; el canvas posee su propio `QUndoStack`.
- `CoordinateSet3D` es un modelo inmutable con coordenadas por ID de átomo. El estado actual de CompChem, `_compchem_coordset`, vive en la ventana y se descarta al cambiar documento; no es suficiente para conformaciones independientes por pestaña.
- `CanvasStructureMixin.get_persistence_data()` se integra en `PersistenceManager.save_to_dict()` bajo la clave raíz `canvas`. Un subcampo opcional puede añadirse sin cambiar los datos de `model` ni el número raíz `version` (`0.1.0`). Al cargar un archivo anterior sin el subcampo, el documento debe seguir siendo válido.
- `AddBondCommand` ya crea un átomo opcional y un enlace en una sola operación undoable; `AddChainCommand` y `AddRingCommand` construyen topologías 2D. La edición 3D debe usar esos comandos o comandos compuestos, no duplicar otro grafo.
- `CompChem3DController` crea un `QThread`, pero `rdkit_backend.optimize_iter()` ejecuta `optimize()` y sólo itera los frames que `optimize_3d_isolated()` devuelve al concluir `_run_worker()`/`subprocess.run()`. Para Live MM hace falta streaming real, no replay de la lista final.
- `scene_molecule_from_graph()` ignora coordenadas ausentes y ya trata los H/enlaces generados como metadata efímera. Esa propiedad se conserva.

## Límites de arquitectura

1. **Modelo único.** `MolGraph` sigue siendo autoridad de conectividad, estereoquímica y propiedades; `CoordinateSet3D` contiene sólo la geometría 3D. Cámara, selección y modo de interacción quedan en el widget/dock. `geometry3d` no importa Qt ni la GUI; el worker sólo recibe/retorna datos primitivos serializables.
2. **Propiedad por documento.** Cada canvas mantiene su coordenada 3D y sus revisiones. La ventana renderiza la coordenada del canvas activo; cambiar pestaña no copia ni descarta la geometría de otro documento. La coordenada de un worker se asocia a token de documento, revisión del grafo y revisión de coordenadas.
3. **Edición aislada de navegación.** `Edición 3D` inicia desactivada. En modo edición, el picking existente elige primitivas y los gestos geométricos no cambian la cámara. Con edición desactivada se preservan los gestos de navegación actuales. `Fit` y proyección sólo cambian cámara.
4. **Operaciones geométricas puras.** Añadir helpers en `geometry3d` para vector/rayo-plano, distancia, ángulo, torsión, traslación rígida y rotación de Rodrigues. Cada helper valida IDs, geometría finita y degeneraciones antes de producir posiciones nuevas. No se modifica la coordenada hasta confirmar el gesto.
5. **Comandos undoables.** Una pulsación-arrastre-confirmación es un solo `QUndoCommand` con snapshot de coordenadas antes/después. Distancia, ángulo, torsión y movimiento de fragmento son una transacción cada uno. La construcción es un comando compuesto/macro que agrupa la operación `MolGraph` existente con la geometría 3D. Frames MM son preview efímero y no se empujan individualmente al stack. Al confirmar el gesto o terminar una relajación se guarda el snapshot final con un comando; Escape o desactivar Live MM cancela y restaura el snapshot committed.
6. **Ciclos y partes rígidas.** Un dihedro sólo rota un lado si, al retirar el enlace central del grafo, sus extremos quedan desconectados. Se exige enlace covalente simple, no aromático, no query y no restringido. Si los extremos siguen conectados —por ejemplo, un enlace de anillo— se deshabilita la operación y se explica el motivo; no se corta el grafo ni se retuerce accidentalmente todo el ciclo. Selección de componente se obtiene por DFS/BFS sobre `MolGraph`; una traslación rígida aplica el mismo vector a todos sus IDs.
7. **Sincronización 2D/3D.** La señal/ciclo de cambios de `undo_stack` compara una firma química/topológica que omite `Atom.x/y`. Un cambio sólo de layout 2D conserva la geometría. Si cambia conectividad o propiedades químicas, la geometría existente se conserva en el canvas pero se marca obsoleta, se oculta y se deshabilitan optimización/exportación; undo/redo que devuelve la firma compatible vuelve a mostrar exactamente ese snapshot. Átomos 2D nuevos no reciben coordenadas inventadas: se solicita construirlos en 3D o regenerar el conformero. Los H RDKit siguen siendo metadata efímera y se recalculan en el siguiente trabajo RDKit.

## Persistencia `.cmsn`

Añadir, cuando existan coordenadas, un objeto opcional `canvas.geometry3d`:

- `schema_version: 1`
- `source`, `method`, `forcefield` y `energy_kcal_mol` (energía puede ser `null`)
- `positions`: lista ordenada de `{atom_id, xyz}` para IDs presentes en `MolGraph`; XYZ son números finitos en Å

No se persiste cámara, selección, flags de edición/live, resultados parciales, H generados ni coordenadas de IDs ya eliminados. Al cargar, se aceptan archivos sin `geometry3d`; se filtran/rechazan entradas malformadas, no finitas, duplicadas o de IDs desconocidos con diagnóstico controlado, sin fallar la carga 2D. Se conserva `PersistenceManager.VERSION == "0.1.0"`: el nuevo subcampo es opcional y autoversionado; una aplicación antigua puede ignorarlo y perderlo al volver a guardar, pero no debe perder la química 2D.

## Interacción y herramientas

- Barra del dock: activar/desactivar edición; herramienta (`Mover`, `Construir`); elemento; orden del nuevo enlace; controles numéricos para distancia (Å), ángulo (°) y dihedro (°); campo Live MM (MMFF94/MMFF94s/UFF) y toggle de optimización.
- Clic selecciona átomo localmente; `Ctrl+clic` selecciona su componente conectado y `Shift+clic` extiende/togglea la selección. Arrastrar la selección mueve rígidamente sus IDs en un plano paralelo a cámara a la profundidad de inicio.
- La distancia, el ángulo y el dihedro se editan numéricamente con IDs explícitos; el dihedro requiere cuatro átomos consecutivos y se rechaza en ciclos o enlaces no rotables. Los controles indican Å/grados y rechazan valores no finitos o geometría degenerada. Los gestos directos de estiramiento/ángulo quedan como extensión futura.
- En construcción se elige elemento y orden; clic en espacio crea un átomo libre, o clic sobre un átomo y después en espacio crea un átomo enlazado. La operación usa comandos de canvas y el validador existente. Se mantienen IDs y no se crean H implícitos documentales automáticamente.

## Live Molecular Mechanics — flujo real

1. Reutilizar M04 (`rdkit_safe.py` y `_rdkit_worker.py`) y el modelo `OptimizationSettings`; no ejecutar RDKit ni esperar procesos desde el hilo GUI.
2. Añadir modo streaming NDJSON al worker aislado: emitir y `flush()` un evento `frame` tras cada bloque de iteraciones y un evento terminal `done/error`; cada frame lleva secuencia, energía kcal/mol y posiciones identificadas por IDs. El worker usa MMFF94 como preferencia, permite MMFF94s y puede pasar a UFF sólo si la política automática está activa y lo informa. Si ningún campo asigna parámetros, el estado es `no_forcefield` y no se simula progreso.
3. Para el átomo o fragmento sujeto al puntero, añadir puntos fijos al force field RDKit; se verificó `ForceField.AddFixedPoint` en la versión RDKit local y la API oficial también describe `AddFixedPoint`, `CalcGrad` y la energía en kcal/mol. Quitar la restricción al soltar. La minimización reduce la energía de la estructura restante; no aplica una dinámica temporal ni un termostato.
4. La GUI aplica cada frame sólo si `job_id`, token de documento y `document_revision` siguen vigentes. El primer movimiento inicia un snapshot tras 80 ms; si el usuario vuelve a mover mientras hay un worker activo, éste se cancela y sólo se conserva el preview más reciente. El snapshot pendiente se inicia después de confirmar la terminación del worker anterior, por lo que nunca hay dos procesos Live MM simultáneos.
5. El supervisor debe terminar el proceso al cancelar/cerrar, esperar un plazo acotado, matar si no termina, cerrar pipes y evitar eventos después de destruir el dock. No se conecta `QThread.requestInterruption()` como única cancelación de un `subprocess.run()` bloqueante.
6. La energía se etiqueta siempre `kcal/mol`. Mostrar fuerza de campo, iteración, estado (`Preparando`, `Relajando`, `Convergió`, `Cancelado`, `Sin parámetros`, `Error`) y energía disponible. Los frames se muestran en preview y no se escriben en undo/disco durante el drag; confirmar el gesto crea un único comando de coordenadas y cancelar restaura el snapshot committed.

## Riesgos

- **Estructura 2D no parametrizable.** Conservar el documento y los datos 3D previos; indicar incompatibilidad del force field y permitir UFF si el usuario/política lo acepta.
- **Worker nativo se cuelga o falla.** Subproceso separado, timeout, terminación acotada y sin bloqueo de GUI.
- **Respuesta obsoleta.** IDs de trabajo y revisiones verificadas en el hilo UI justo antes de aplicar cada frame.
- **Undo/redo de topología.** Mantener IDs estables; una firma compatible restaura el snapshot 3D tras undo, mientras una topología divergente lo oculta sin destruirlo y exige regeneración antes de operar.
- **Anillo o conformación degenerada.** No intentar cortar enlaces cíclicos; devolver error estable en vez de posiciones parciales.
- **Regresión persistencia/renderer.** Campo opcional con valor por defecto; conservar viewport read-only y formato existente; no modificar cámara desde comandos químicos.
- **SIGSEGV suite.** Baseline de esta rama termina `pytest -q` con código `-11` y dos marcadores `F` antes del crash. No alterar ni ocultar ese baseline; reportar validación focalizada y volver a registrar la suite al final, sin convertir esta campaña en una investigación general del apagado.

## Archivos previstos

- Existentes: `src/chemuson/geometry3d/model.py`, `service.py`, `rdkit_backend.py`, `chemio/rdkit_safe.py`, `chemio/_rdkit_worker.py`, `gui/viewport3d/{widget.py,dock.py}`, `gui/canvas/{canvas_view.py,canvas_structure.py}`, `gui/main_window.py`, `gui/shell/assembly.py`, `chemio/persistence.py`.
- Nuevos dentro de paquetes existentes: helpers puros de edición geométrica en `geometry3d`; comandos y controller del editor dentro de `gui/commands` y `gui/viewport3d` (sin nuevo paquete raíz).
- Pruebas: `tests/test_geometry3d_*`, `tests/test_viewport3d_*`, `tests/test_compchem3d_dock.py`, persistencia y nuevos casos de integración; actualizar `architecture/modules.yml` únicamente si cambian imports/ownership.
- Evidencia: `avogadro-1.2.0-analysis.md`, `baseline.md` y `baseline/*.log`.

## Verificación necesaria

Pruebas RED/GREEN por operación pura, comando undo, persistencia, worker streaming/cancelación/revisión, widget interactivo y sincronización 2D; prueba real en display con etano/butano/benceno; registrar para butano geometría/energía a través de varios frames antes de convergencia; `compileall`, Ruff selectivo, tests focalizados, `openspec validate ... --strict` y `git diff --check`. Ejecutar y reportar `pytest -q` una vez al final; si repite SIGSEGV se informa con el baseline, no se declara aprobado global.
