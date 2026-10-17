# Tasks — fix-3d-dock-tools-and-conformers

## Hito A — Dock flotante y estado de edición

- [x] A1. QSS/paleta clara/oscura específica del dock; controles, disabled/checked/hover legibles; reaccionar a tema y `topLevelChanged` sin tocar KDE ni otros docks.
- [x] A2. Reflejar modo/herramienta coherentemente desde inicio; elegir herramienta activa edición y navegación queda claramente diferenciada.
- [x] A3. Probar estado acoplado/flotante, retema y preservación de tool/escena/cámara/selección/contexto a nivel de wiring offscreen.
- [x] A4. Commit local Hito A y verificación de regresiones focalizadas.

## Hito B — Interacción y handles

- [x] B1. Highlight visible de átomos documentales/derivados y enlaces; selección múltiple ordenada y selección de fragmento.
- [x] B2. Incorporar cambio de elemento seleccionado mediante comando undoable y validar valencia/estructura antes de confirmar; cambio de orden de enlace mediante gesto entre extremos.
- [x] B3. Desplazamiento en profundidad con Alt+arrastre o botón central; transformar H derivados al mover átomos/fragmentos; navegación con modificadores documentados.
- [x] B4. Handles undoables de longitud, ángulo y dihedro; rechazar ciclos/geometrías degeneradas y restaurar preview con Escape.
- [x] B5. Tests de interacción/señales, valencia, cancelación y undo/redo; etano/etanol/butano; commit local Hito B.

## Hito C — Generación robusta y diagnósticos

- [ ] C1. Diagnósticos estructurados por etapa/conteos/estereoquímica/tiempos/código RDKit y errores diferenciados en viewport + panel CompChem.
- [ ] C2. Intentos ETKDG finitos/diferenciados; opción macro ciclo condicionada; validación geométrica/estereoquímica; resultado no optimizado claro ante fallo MMFF/UFF.
- [ ] C3. Progreso por etapas real, timeout y cancelación del subproceso aislado; descartar resultados obsoletos sin ocultar estado del documento actual.
- [ ] C4. Regresión auténtica de tetrandrina CID 73078 desde la referencia documentada, además de etano, butano, benceno y otro policíclico; preservar pruebas estéreo/aromaticidad.
- [ ] C5. Verificar escena 3D preparada y solicitud de upload/render; registrar como no visual si no se inicia GUI.
- [ ] C6. Ejecuciones focalizadas, arquitectura, Ruff scoped, compileall, OpenSpec strict, `git diff --check`; commit local Hito C y estado Git limpio sin publicar.

## Aceptación manual separada

- [ ] En CachyOS, comprobar con compositor KDE los temas claro/oscuro, flotar/mover/redimensionar/acoplar el dock y confirmar contexto OpenGL, colores y escena; no es sustituible por pruebas headless.
- [ ] Probar el archivo `.cmsn` original de tetrandrina del usuario cuando esté disponible; la referencia SMILES no se considera el documento original.
