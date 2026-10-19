# Tasks — editor molecular 3D y mecánica molecular en vivo

## 1. Fase 3: investigación y OpenSpec

- [x] 1.1 Confirmar `feature/chemuson-3d`, HEAD, worktree y baseline exacto antes de editar.
- [x] 1.2 Inspeccionar modelos, comandos, persistencia, viewport, workers, sincronización de pestañas y catálogo arquitectónico existentes.
- [x] 1.3 Analizar fuente Avogadro 1.2.0 y binarios extraídos con REA/Ghidra sin ejecutar; redactar `avogadro-1.2.0-analysis.md` separando fuente/inferencia/decompilación.
- [x] 1.4 Crear propuesta, diseño, tareas y specs con contratos de coordenadas, interacción, validación, persistencia, Live MM y undo/redo.
- [x] 1.5 Registrar salidas exactas de compileall, colección, suite y Ruff en `baseline/`.
- [x] 1.6 Validar OpenSpec con `openspec validate add-editable-3d-molecular-editor --strict`; resolver errores de esquema antes de producción.

## 2. Estado 3D por documento y persistencia

- [x] 2.1 Añadir estado de coordenadas/revisión y token de documento a cada `ChemusonCanvas`, independiente de cámara y 2D.
- [x] 2.2 Implementar persistencia opcional de coordenadas 3D; validar finitud, IDs, round-trip y archivos antiguos.
- [x] 2.3 Comparar firma química/topológica; ocultar y conservar geometría obsoleta, restaurarla si undo vuelve a firma compatible y preservar cambios sólo de layout 2D.
- [x] 2.4 Probar aislamiento por pestaña, compatibilidad `.cmsn`, estados inválidos y sincronización de revisión.

## 3. Geometría editable y transacciones

- [x] 3.1 Añadir helpers puros de distancia, ángulo, torsión, traslación rígida, plano de drag y rotación de Rodrigues.
- [x] 3.2 Detectar puentes por eliminación temporal de arista; impedir torsión independiente en ciclos, enlaces no simples/aromáticos/query.
- [x] 3.3 Añadir comandos 3D de coordenadas que guardan snapshots y un solo paso undo por gesto.
- [x] 3.4 Añadir controles de ratón y valores numéricos (Å/grados), manteniendo cámara independiente.
- [x] 3.5 Probar manipulación de átomo/fragmento, distancias, ángulos, torsiones, ciclos y undo/redo.

## 4. Constructor y coherencia química

- [x] 4.1 Añadir selección de elemento y construcción 3D de átomo libre/enlazado.
- [x] 4.2 Reusar `MolGraph`, validación y comandos de `gui/commands`; agrupar topología + coordenadas como una sola acción.
- [x] 4.3 Reconciliar topología 2D con firma; no mostrar ni optimizar snapshots incompatibles, preservarlos para undo y pedir regeneración; H se reconstruyen en el siguiente trabajo RDKit.
- [x] 4.4 Completar las pruebas negativas del constructor (valencia, endpoint inexistente y duplicado). Se verificó que los rechazos conservan grafo, geometría 3D e historial undo; la cobertura previa incluye construcción atómica de átomo+enlace, orden doble, H derivados no persistentes, sincronización 2D/3D y undo/redo de topología.

## 5. Live Molecular Mechanics real

- [x] 5.1 Añadir protocolo NDJSON con frames emitidos y flush antes del resultado terminal; conservar aislamiento del proceso RDKit.
- [x] 5.2 Añadir minimización incremental MMFF94/MMFF94s y fallback UFF explícito; reportar kcal/mol y falta de parámetros.
- [x] 5.3 Añadir restricciones temporales `AddFixedPoint` para IDs bajo arrastre.
- [x] 5.4 Añadir debounce/coalescencia, un worker activo, tokens de documento/revisiones y descarte de frames obsoletos.
- [x] 5.5 Implementar cancel-restore desde el toggle Live MM, terminación acotada del subproceso y cancelación al cerrar/cambiar documento.
- [x] 5.6 Probar frames antes de `done`, energía/coordenadas progresivas, cancelación, expiración, error, átomo fijado y descarte de resultado obsoleto.

## 6. Integración, regresión y demostración

- [x] 6.1 Integrar controles en `MolecularViewportDock` y ciclo de vida de `ChemusonWindow` sin crear otro modelo ni paquete raíz.
- [x] 6.2 Ejecutar tests focalizados de geometría, worker, persistencia, comandos, GUI y arquitectura; Ruff scoped y compileall. Evidencia exacta en `validation.md`.
- [ ] 6.3 Prueba gráfica local de butano en `DISPLAY=:0`: el visor OpenGL mostró 14 átomos/13 enlaces; la herramienta fijó C–C–C–C a 60° y RDKit/MMFF94 entregó 30 frames progresivos distintos (energía 281.439546 → −4.293763 kcal/mol). La captura y el informe JSON se guardaron en scratch. La ejecución terminó con SIGSEGV durante liberación de wrappers SIP/QApplication; la causa no está determinada. Etano/benceno no se verificaron gráficamente en esta sesión. Repetición adicional quedó bloqueada; no reintentar desde el agente.
- [x] 6.4 Ejecutar `openspec validate add-editable-3d-molecular-editor --strict`, `git diff --check`, pruebas focalizadas, compileall y Ruff scoped. No repetir `pytest -q`: el baseline registrado ya termina con SIGSEGV; no se atribuye ese fallo global al cambio.
- [x] 6.5 Registrar instrucciones mínimas para CachyOS, archivos/decisiones/limitaciones en `validation.md`; no push/merge/PR.
