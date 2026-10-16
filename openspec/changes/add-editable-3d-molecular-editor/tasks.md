# Tasks — editor molecular 3D y mecánica molecular en vivo

## 1. Fase 3: investigación y OpenSpec

- [x] 1.1 Confirmar `feature/chemuson-3d`, HEAD, worktree y baseline exacto antes de editar.
- [x] 1.2 Inspeccionar modelos, comandos, persistencia, viewport, workers, sincronización de pestañas y catálogo arquitectónico existentes.
- [x] 1.3 Analizar fuente Avogadro 1.2.0 y binarios extraídos con REA/Ghidra sin ejecutar; redactar `avogadro-1.2.0-analysis.md` separando fuente/inferencia/decompilación.
- [x] 1.4 Crear propuesta, diseño, tareas y specs con contratos de coordenadas, interacción, validación, persistencia, Live MM y undo/redo.
- [x] 1.5 Registrar salidas exactas de compileall, colección, suite y Ruff en `baseline/`.
- [x] 1.6 Validar OpenSpec con `openspec validate add-editable-3d-molecular-editor --strict`; resolver errores de esquema antes de producción.

## 2. Estado 3D por documento y persistencia

- [ ] 2.1 Añadir estado de coordenadas/revisión a cada `ChemusonCanvas`, independiente de cámara y 2D.
- [ ] 2.2 Implementar `canvas.geometry3d` opcional versión 1; validar finitud, IDs, round-trip y archivos antiguos.
- [ ] 2.3 Reconciliar ediciones 2D/undo/redo preservando coordenadas de IDs supervivientes y metadata efímera de H.
- [ ] 2.4 Escribir tests RED para aislamiento por pestaña, compatibilidad `.cmsn`, estados inválidos y sincronización.

## 3. Geometría editable y transacciones

- [ ] 3.1 Añadir helpers puros de distancia, ángulo, torsión, traslación rígida, plano de drag y rotación de Rodrigues.
- [ ] 3.2 Detectar puentes por eliminación temporal de arista; impedir torsión independiente en ciclos, enlaces no simples/aromáticos/query.
- [ ] 3.3 Añadir comandos 3D de coordenadas que guardan snapshots y un solo paso undo por gesto.
- [ ] 3.4 Añadir controles de ratón y valores numéricos (Å/grados), manteniendo cámara independiente.
- [ ] 3.5 Añadir tests RED/GREEN de átomo, componente, distancias, ángulos, butano torsional, ciclo, cancelación y undo/redo.

## 4. Constructor y coherencia química

- [ ] 4.1 Añadir selección de elemento/orden de enlace y construcción 3D de átomo libre/enlazado.
- [ ] 4.2 Reusar `MolGraph`, validación y comandos de `gui/commands`; agrupar topología + coordenadas como una sola acción.
- [ ] 4.3 Actualizar escena al cambiar el grafo desde 2D; conservar IDs/coordenadas supervivientes y recalcular sólo H derivados.
- [ ] 4.4 Probar valencia/endpoints/duplicados, H implícitos, coherencia 2D/3D y undo/redo de topología.

## 5. Live Molecular Mechanics real

- [ ] 5.1 Añadir protocolo NDJSON con frames emitidos y flush antes del resultado terminal; conservar aislamiento del proceso RDKit.
- [ ] 5.2 Añadir minimización incremental MMFF94/MMFF94s y fallback UFF explícito; reportar kcal/mol y falta de parámetros.
- [ ] 5.3 Añadir restricciones temporales `AddFixedPoint` para IDs bajo arrastre.
- [ ] 5.4 Añadir debounce/coalescencia, un worker activo, tokens de documento/revisiones y descarte de frames obsoletos.
- [ ] 5.5 Implementar stop-accept/cancel-restore, terminación acotada del proceso y limpieza al cerrar/cambiar documento.
- [ ] 5.6 Probar frames que llegan antes de `done`, energía/coords progresivas, cancelación, expiración, error, pinned atom y resultado tardío.

## 6. Integración, regresión y demostración

- [ ] 6.1 Integrar controles en `MolecularViewportDock` y ciclo de vida de `ChemusonWindow` sin crear otro modelo ni paquete raíz.
- [ ] 6.2 Ejecutar tests focalizados de geometría, worker, persistencia, comandos, GUI y arquitectura; Ruff scoped y compileall.
- [ ] 6.3 Probar gráficamente etano, butano y benceno; para butano registrar varios frames reales de coordenadas/energía y captura de torsión/relajación.
- [ ] 6.4 Ejecutar `openspec validate ... --strict`, `git diff --check` y `pytest -q`; documentar el SIGSEGV si repite y no atribuirlo sin evidencia.
- [ ] 6.5 Escribir instrucciones mínimas para CachyOS, lista de archivos/decisiones/limitaciones; no push/merge/PR.
