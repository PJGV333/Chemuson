# Informe técnico — Avogadro 1.2.0 y mecánica molecular interactiva

Fecha de inspección: 2026-10-16. Este documento separa evidencia de código fuente, observaciones de Ghidra y conclusiones de diseño. No es una prueba de ejecución.

## Método, identidad y licencia

- Checkout estático: tag `1.2.0`, commit `aaa67d012828b35314049a9463b77f70415c29c3` (2016-06-08); árbol de fuente limpio y detached en scratch.
- AppImage extraído en `/home/ccachyavgp/.hermes/cache/scratch/avogadro-1.2.0-appimage-root`; no se ejecutó Avogadro ni sus bibliotecas. `COPYING` contiene GNU GPL versión 2; los encabezados de las herramientas confirman GPLv2 o GPLv2-or-later.
- SHA-256 de binarios extraídos usados en el análisis: `autoopttool.so` `bd18a6d7559d68527469e4d81f2e7a9884a53cda72c91d73d45783928a0488c8`; `manipulatetool.so` `38408987a3cd74636dd974e235c84200697b952c4fbc9c0daf9d4e434c32e999`; `bondcentrictool.so` `160842b2d7d20acb97805cd1e86aa47ae4ac1e6344cb4a166a33bcb7d253db23`; `libavogadro.so.1.2.0` `f9c1f3c1d9e5c8aa3f5aa26d4eb188f2a4a64af0ec25b6b4a4e7b2e3483b37ba`.
- Ghidra 12.1.4/REA se usaron sobre ELF extraídos. No se modificaron bytes ni se copió código GPL. Las conclusiones funcionales abajo proceden principalmente de lectura estática de la fuente taggeada; una descompilación no se presenta como ejecución.

## Arquitectura funcional observada

### AutoOptTool: interacción y bucle de optimización

**Demostrado por fuente** (`libavogadro/src/tools/autoopttool.cpp`):

- El panel deja escoger un force field de Open Babel, algoritmo y “Steps per Update” (1–100, por defecto 4). El selector inicial usa UFF; los algoritmos disponibles son Steepest Descent, Conjugate Gradients y Molecular Dynamics a 300, 600 o 900 K (líneas 267–333, 541–556).
- Al activar la herramienta inicia un `AutoOptThread`, empuja un `AutoOptCommand` a `QUndoStack` y usa un timer Qt de 50 ms para solicitar actualizaciones (353–424, 467–471).
- Con optimización activa, arrastrar un átomo fija temporalmente ese punto mediante `SetFixAtom`; al soltar llama `UnsetFixAtom`. El movimiento transforma el átomo/selección y los vecinos responden a la siguiente minimización (133–209, 104–131, 173–185).
- Al acabar cada actualización, la UI obtiene coordenadas (`GetCoordinates`) y, si existe `ConformerData`, vectores de fuerza; copia fuerzas y posiciones a los átomos y solicita actualización del modelo/widget (427–465).
- `AutoOptCommand` guarda una copia inicial, combina ejecuciones sucesivas, y undo desactiva la herramienta y restaura la molécula (569–601).
- El backend de AutoOpt usa `OpenBabel::OBForceField`; el constructor elige UFF inicialmente (70–84), aunque el combo se llena con los force fields registrados.

**Inferencia estática importante, no medida en runtime:** `AutoOptTool::timerEvent()` llama directamente `m_thread->update()` (405–424). `update()` es una llamada ordinaria que ejecuta `Setup`, `SteepestDescent`, `ConjugateGradients` o MD (510–562); `AutoOptThread::run()` sólo llama a `exec()` (505–508). Por la semántica de una llamada directa al objeto `QThread`, esta ruta no demuestra que la parte costosa ocurra en el hilo de trabajo y puede ejecutarse en el hilo que entrega el timer. No se midió bloqueo porque se prohibió ejecutar el binario. No copiar ese patrón.

### ForceFieldExtension: worker incremental distinto

**Demostrado por fuente** (`libavogadro/src/extensions/forcefieldextension.cpp`):

- La extensión estándar inicia con MMFF94, expone optimización, energía, búsqueda de conformeros, restricciones, fijar átomos e ignorar selección (60–132).
- `ForceFieldThread::run()` configura el campo y ejecuta el algoritmo en un `QThread`; en steepest descent toma bloques de 5 pasos, llama `UpdateCoordinates`, copia fuerzas/posiciones bajo un lock y emite `stepsTaken` después de cada bloque (395–507).
- `ForceFieldCommand` posee un snapshot de la molécula; `redo()` inicia el thread, cancelar llama `stop()`, y `undo()` detiene/espera y restaura el snapshot (596–705).

Esto constituye una referencia útil para worker acotado, feedback incremental, restricciones y transacción undoable. ChemUSON conserva su aislamiento de proceso RDKit, no adopta el modelo Open Babel de Avogadro.

### ManipulateTool, selección y cámara

**Demostrado por fuente:**

- `ManipulateTool` usa picking de átomo: left-drag traslada el átomo elegido; con una selección mueve el centroide/selección. Los otros botones permiten inclinación/zoom centrados o rotación en torno al átomo/centro seleccionado. Press/release devuelven comandos de movimiento (104–131, 289–439).
- `SelectRotateTool` usa hits; el clic selecciona y el doble clic recorre el componente conectado mediante DFS de Open Babel y selecciona átomos/enlaces de ese componente (112–175, 411–504).
- `GLWidget::hits`, `computeClickedPrimitive` y `computeClickedAtom` implementan picking; `Camera::project/unProject` convierte entre coordenadas de ventana y espacio, y el drag fija el punto de referencia al plano a profundidad constante (`glwidget.cpp:2224–2361`, `camera.cpp:277–299`).
- El rectángulo de selección usa `gluUnProject`, matriz OpenGL inmediata y dibujo legado (`selectrotatetool.cpp:511–555`). Es evidencia del diseño antiguo, no una recomendación de implementación.

### BondCentricTool, distancias, ángulos y torsiones

**Demostrado por fuente** (`bondcentrictool.cpp`):

- Primero se selecciona un enlace. Un clic sobre uno de sus extremos construye un `SkeletonTree` para operar sobre una parte de la molécula (183–260).
- Arrastre izquierdo del extremo seleccionado produce rotación/ajuste angular del esqueleto; hay visualización de ángulos y ajuste a referencias angulares (`mouseMoveEvent`, 426–622; dibujo de sectores/dihedros, 855–1247).
- Arrastre derecho sobre el extremo altera la longitud proyectada del enlace; arrastrar un átomo vecino al enlace central manipula un dihedro y rota el esqueleto alrededor del eje central (624–718).
- `SkeletonTree::populate()` hace dos recorridos desde extremos opuestos excluyendo el enlace seleccionado y usa `containsAtom()` para no visitar un átomo repetido (`skeletontree.cpp:142–204`). `skeletonRotate/Translate` transforma posiciones de los nodos (`208–263`).

**Límite del manejo de ciclos:** el código evita recursión/duplicación mediante árboles visitados, pero no contiene una regla explícita de “este enlace no es puente, no existe un lado independiente que pueda rotarse”. El recorrido puede asignar átomos de un ciclo a uno de los árboles según el orden de visita; esto no equivale a resolver un anillo como restricción conformacional. ChemUSON debe detectar puentes en su propio grafo y rechazar torsiones independientes en enlaces cíclicos, salvo que incorpore en otro cambio un solver de cierre de anillo.

## Evidencia REA/Ghidra y límites

- Para `manipulatetool.so`, Ghidra importó el ELF de SHA-256 indicado y produjo una vista pseudocódigo de `0x108bc0` (`FUN_00108bc0`, 5,116 unidades UTF-16). Evidencia de análisis padre `ev_1816dc1621dcb495d684959acbefbe9a9815c2eeb8c8e72e7ff29a24fce381ae`; vista `ev_11db7064af4c8ef1af3f3419bb18c58d02551728953a5998ee701a46a4e03dc3`.
- Para `bondcentrictool.so`, se inspeccionó la función candidata `0x107d50` (`FUN_00107d50`, 24,051 unidades UTF-16; salida parcial/truncable). Evidencia padre `ev_1d55c5e54102264e344c847db86a2ce8077943344b6187b0a36c2f61574fa6d6`; vista `ev_0be76fcf7c98b279ef5de4c6f4fb57ebe9db3fbc494bf5be58a4dcc5d1a58355`.
- Para `autoopttool.so`, Ghidra indexó importaciones `SteepestDescent`, `ConjugateGradients`, `SetFixAtom`, `UnsetFixAtom`, `setPos`, `setForceVector`, `setEnergy`, `startTimer` y strings de algoritmos/UI. El script `DumpAutoOptEvidence.java` produjo cuatro funciones ligadas a strings; no recuperó una decompilación semántica verificable de `AutoOptThread::update()`/despacho del campo de fuerza. Resultado local: `/home/ccachyavgp/.hermes/cache/scratch/avogadro-autoopt-ghidra-decomp-v2.txt`.
- Los nombres `FUN_...` y correlaciones por dirección no prueban qué método C++ representa cada función. P-code/pseudocódigo son derivados, truncables y no comportamiento en ejecución. Por ello REA sirve aquí como comprobación estructural/call-surface, no como autoridad sobre el algoritmo.

## Mapa de equivalencias Avogadro → ChemUSON

| Necesidad | Avogadro observado | Equivalente propio recomendado |
|---|---|---|
| Pick y proyección | `GLWidget` hits + `Camera.project/unProject` | picking CPU y cámara existente; plano de drag derivado del rayo y profundidad inicial |
| Mover átomo/componente | `ManipulateTool`, selección local/centroide | modo edición explícito y `CoordinateSet3D` por canvas; un comando de geometría por gesto |
| Componente conectado | doble clic DFS en `SelectRotateTool` | DFS/BFS sobre `MolGraph`, sin depender de Open Babel |
| Distancia/ángulo/dihedro | `BondCentricTool` + `SkeletonTree` | transforms geométricos puros; controles Å/grados; bridge test para torsión |
| Live MM | AutoOpt Open Babel y extensión de force field | RDKit aislado en flujo NDJSON, snapshots reales durante ejecución, revisiones y cancelación |
| Mantener átomo bajo puntero | `SetFixAtom` temporal | `ForceField.AddFixedPoint()` sobre IDs mapeados al índice RDKit durante drag |
| Undo de relajación | snapshot de `AutoOptCommand`/`ForceFieldCommand` | previews fuera del stack; aceptar/parar guarda una sola transacción de coordenadas |
| Render de fuerzas | vectores `ConformerData` | dejar vector overlay fuera de MVP; energía/iteración bastan para probar relajación real |

## Recomendaciones concretas

1. Mantener `MolGraph` como única conectividad y mover la coordenada 3D de estado global de ventana a cada `ChemusonCanvas`; persistir sólo coordenadas de átomos documentales.
2. Separar edición geométrica de cámara y conservar read-only como modo inicial; `Fit` no debe tocar la conformación.
3. Para medir distancias/ángulos/dihedros, usar funciones puras testeadas; al editar una unión, transformar el lado completo del puente. No aplicar torsión independiente sobre un enlace cuyo corte no desconecta el grafo.
4. Crear en M04 un worker aislado que emita cada bloque y haga flush. El backend actual almacena frames hasta `subprocess.run()` final, por lo que `optimize_iter()` existente no cumple el requisito de vivo.
5. Congelar con `AddFixedPoint` los IDs bajo arrastre; coalesce en ~100 ms, mantener un único proceso activo y rechazar frames que no coincidan con token de documento y revisiones de grafo/coordenadas.
6. Mantener MMFF94 preferido; permitir MMFF94s/UFF visible. No reportar éxito cuando faltan parámetros. RDKit devuelve `CalcEnergy()` en kcal/mol y expone `AddFixedPoint()`/`CalcGrad()` en la documentación oficial y en la instalación local.
7. No portar `GL_SELECT`, `gluUnProject`, dibujo inmediato, `SkeletonTree`, QThread de Qt4 ni código fuente GPL. No se reprodujo visualmente Avogadro; toda conclusión es estática.
