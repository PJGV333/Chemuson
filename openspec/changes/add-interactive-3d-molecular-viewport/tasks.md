# Tasks — visor molecular 3D interactivo

## 1. OpenSpec y baseline

- [x] 1.1 Verificar `feature/chemuson-3d`, worktree independiente, HEAD y árbol limpio antes de modificar.
- [x] 1.2 Leer los artefactos OpenSpec activos, `geometry3d`, `MolGraph`, GUI y catálogo arquitectónico.
- [x] 1.3 Registrar baseline de compileall, colección, suite completa y Ruff en `baseline.md`.
- [x] 1.4 Validar este cambio con `openspec validate add-interactive-3d-molecular-viewport --strict` antes de implementar producción.

## 2. Núcleo puro del visor

- [x] 2.1 Añadir estado de cámara acotado con perspectiva/ortográfica, orbit, pan, zoom y fit estable.
- [x] 2.2 Añadir geometría de rayos, intersecciones esfera/cápsula y picking determinista.
- [x] 2.3 Añadir generación de mallas propia para esferas y varillas, con descarte de degenerados.
- [x] 2.4 Añadir pruebas unitarias de invariantes, valores inválidos, moléculas vacías y ausencia de mutaciones.

## 3. Adaptación OpenGL/PyQt6

- [x] 3.1 Implementar renderer shader/VBO/VAO sobre `QOpenGLWidget` y `QtOpenGL`, sin PyOpenGL, `GL_SELECT` ni renderizado inmediato.
- [x] 3.2 Implementar interacción de ratón/rueda, selección local, perspectiva/ortográfica y encuadre.
- [x] 3.3 Implementar `MolecularViewportDock` con degradación visible si el contexto OpenGL falla.
- [x] 3.4 Añadir pruebas GUI que no requieran un contexto real y pruebas de wiring del dock.

## 4. Integración ChemUSON

- [x] 4.1 Montar el dock independiente oculto por defecto y añadir `Ver → Visor molecular 3D`.
- [x] 4.2 Reutilizar `SceneMolecule3D` y el `CoordinateSet3D` del controlador CompChem; actualizar al terminar generación.
- [x] 4.3 Limpiar la escena al cambiar de documento y demostrar que el canvas 2D, `MolGraph`, selección documental y undo no cambian.
- [x] 4.4 Actualizar `architecture/modules.yml` y contratos de tests M08 sin crear una dependencia nueva desde M05.

## 5. Validación y entrega

- [x] 5.1 Ejecutar tests focalizados, compileall, Ruff scoped, validación OpenSpec estricta y `git diff --check`.
- [x] 5.2 Ejecutar una prueba visual reproducible en el display disponible; verificar la captura de esferas y varilla renderizadas.
- [x] 5.3 Ejecutar un subconjunto de regresión GUI existente; no rerun de la suite monolítica si reproduce el segfault baseline.
- [x] 5.4 Registrar archivos, limitaciones y commits locales separados; no publicar ni modificar otras ramas.

## 6. Estabilización de corrección

- [x] 6.1 Añadir baseline pre-cambio en `baseline.md` para el HEAD `e2a5a3d`; conservar salida, conteo, fallos visibles y segfault como no atribuidos, sin inventar nombres de tests ni causas.
- [ ] 6.2 Cubrir en tests RED la subida GPU diferida, destrucción/recreación de contexto y estado de inicialización fallida.
- [ ] 6.3 Cubrir en tests RED el fit con aspect ratio/radios, auto-fit inicial/resize, preservación tras navegación manual y clipping de escena amplia.
- [ ] 6.4 Cubrir en tests RED cápsulas exactas: primer contacto lateral, hemisferios, origen interior, tangente, paralela y segmento corto.
- [ ] 6.5 Incluir hidrógenos añadidos por RDKit y sus enlaces en `SceneMolecule3D` a partir de metadata efímera; probar no mutación ni persistencia `.cmsn`.
- [ ] 6.6 Ejecutar pruebas focalizadas y smoke visual con etano/benzene; revisar diff y conservar cambios sólo en este worktree/branch.
- [ ] 6.7 Re-ejecutar checks globales; registrar claramente fallos residuales sin atribuirlos al visor salvo evidencia diferencial.
