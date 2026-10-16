# Design — visor molecular 3D interactivo

## Contexto actual

`geometry3d.model` ya define `CoordinateSet3D`, `SceneAtom3D`, `SceneBond3D` y `SceneMolecule3D`; `scene_molecule_from_graph` construye una escena inmutable derivada de `MolGraph`. La generación RDKit y el controlador `CompChem3DController` ya existen. La GUI usa `ChemusonWindow`, `SidePanel` y `CompChemDock`; el nuevo visor no debe duplicar ninguno de esos modelos ni convertir el SidePanel en un octavo tab.

El baseline previo a esta especificación queda registrado en `baseline.md`.

## Decisiones

1. **Separación por capas dentro de M08.**
   - `camera.py`: estado de cámara y operaciones puras de orbit/pan/zoom/fit.
   - `geometry.py`: vectores, matrices/proyección, rayos e intersecciones numéricamente defensivas.
   - `picking.py`: selección de la primitiva más cercana de `SceneMolecule3D`.
   - `mesh.py`: generación propia de esferas y varillas en buffers intercalados; no contiene Qt.
   - `renderer.py`: adaptación OpenGL de PyQt6, shaders, VBO/VAO, profundidad y dibujo.
   - `widget.py`: `QOpenGLWidget`, eventos de ratón/rueda, estado visual de selección y señales.
   - `dock.py`: `QDockWidget` discreto con el widget, mensaje de estado y contrato de actualización de escena.

2. **Renderer Qt6 moderno sin PyOpenGL.** `QOpenGLWidget` fija un formato 4.1 core solicitado, usa shader vertex/fragment, profundidad y buffers. El renderer sólo recibe la escena y la selección; no conoce `MolGraph`, documentos ni undo/redo. Si el contexto o shaders no están disponibles, el widget emite un estado controlado y no intenta renderizado inmediato legado.

3. **Cámara estable y portable.** La cámara se representa con triples inmutables y valores acotados. La proyección perspectiva usa FOV vertical; la ortográfica usa escala vertical. El renderer convierte el estado a `QMatrix4x4` sólo en la frontera Qt. El picking reutiliza exactamente los parámetros de cámara y coordenadas de viewport, evitando una segunda convención de ejes.

4. **Picking CPU.** Un clic construye un rayo en coordenadas de escena. Los átomos se prueban como esferas con el radio visual de ball-and-stick; los enlaces se prueban como cápsulas alrededor del segmento entre sus átomos. Se devuelve la menor distancia positiva y, en empate, se prioriza átomo sobre enlace. No se usa `GL_SELECT`, color picking ni lectura de framebuffer.

5. **Selección local.** `MolecularViewportWidget` mantiene `selected_atom_id` y `selected_bond_id`, emite `primitiveSelected(kind, id)` y recolorea la malla seleccionada. No escribe `MolGraph.selected_atoms`, `MolGraph.selected_bonds`, `ChemState` ni la selección de `ChemusonCanvas`.

6. **Integración incremental.** `MolecularViewportDock` se monta como dock clásico independiente y oculto por defecto. `Ver → Visor molecular 3D` usa su `toggleViewAction`. Al mostrarse, la ventana reutiliza el `CoordinateSet3D` generado por CompChem; si no existe, inicia la generación RDKit mediante el controlador existente. Al terminar un trabajo visible, la ventana actualiza la escena derivada. Cambiar de documento descarta la escena/coords efímeras como ya hace `_reset_compchem_state`.

7. **Datos y degradación.** Moléculas vacías muestran un estado informativo. Escenas sin coordenadas suficientes omiten primitivas incompletas sin mutar el modelo. Un fallo del backend de generación o de OpenGL deja el editor 2D y el documento operativos y muestra el motivo estable en el dock; no se presenta una captura o éxito visual ficticio.

## Endurecimiento de corrección

1. `set_scene` y cambios de selección sólo invalidan el buffer GPU y solicitan repintado. La subida y liberación de recursos se ejecutan con el contexto propietario actual; `aboutToBeDestroyed` libera los recursos y un contexto recreado vuelve a compilar y subir la escena.
2. El fit usa el aspect ratio real y los radios visuales para contener la esfera de bounds completa bajo perspectiva u ortográfica. La escena inicial se encuadra al conocerse el tamaño del viewport y durante resize mientras el usuario no haya navegado; orbit/pan/zoom desactivan ese auto-fit. El fit explícito siempre vuelve a encuadrar. Los planos near/far se derivan del radio y distancia para no truncar moléculas grandes.
3. El picking de enlace busca la primera intersección positiva con la superficie de una cápsula finita: cilindro lateral y hemisferios externos, incluidos rayos que empiezan dentro, tangencias y segmentos cortos. Se compara esa distancia de superficie con las esferas de los átomos.
4. El worker RDKit devuelve coordenadas y enlaces de los hidrógenos añadidos como metadata transitoria; `scene_molecule_from_graph` los agrega a `SceneMolecule3D`. No se asignan IDs al `MolGraph`, no se cambian selecciones/documentos ni se serializa metadata generada en `.cmsn`.

## Archivos afectados

- `src/chemuson/gui/viewport3d/*.py` — nueva implementación M08 y estabilización del renderer/cámara/picking.
- `src/chemuson/geometry3d/model.py`, `src/chemuson/chemio/rdkit_safe.py` y `src/chemuson/chemio/_rdkit_worker.py` — metadata transitoria para hidrógenos y enlaces añadidos por RDKit.
- `src/chemuson/gui/shell/assembly.py` — composición del dock y actualización con trabajos CompChem.
- `src/chemuson/gui/main_window.py` — estado/refresh del visor y limpieza al cambiar de documento.
- `src/chemuson/gui/main_window_ui_builder.py` — entrada `Ver → Visor molecular 3D`.
- `architecture/modules.yml` — rutas, internals y tests de M08.
- `tests/test_viewport3d_*.py` y una regresión GUI focalizada.
- `openspec/changes/add-interactive-3d-molecular-viewport/*` — contrato y evidencia.

## Riesgos y mitigaciones

- **Contexto OpenGL no disponible en CI/offscreen.** Las pruebas de cámara/picking/malla son puras; la prueba GUI no llama `show()` ni exige inicialización GPU. La prueba visual se reportará separadamente.
- **Diferencias de versión Qt/OpenGL.** Se comprueba el contexto y la inicialización de funciones; el dock conserva estado de error y no rompe el shell.
- **Bloqueo al generar el primer conformero.** Se reutiliza `CompChem3DController` y su worker existente, no se ejecuta RDKit desde el evento de pintura.
- **Coordenadas y radios degenerados.** Las funciones puras rechazan NaN/inf, segmentos de longitud casi cero y viewport no positivo de forma determinista.
- **Regresión de documentos.** La escena es una proyección inmutable; los tests comparan IDs, coordenadas químicas, selección 2D y undo antes/después de abrir/cerrar el visor.

## Pruebas necesarias

- Unitarias: vectores, cámara, matrices/proyección, orbit/pan/zoom/fit con aspect ratio y resize, rayos, esfera, cápsula (interior/tangente/paralela/segmento corto), selección más cercana y malla degenerada.
- Integración GUI offscreen: dock oculto por defecto, acción `Ver`, actualización con escena, cierre sin destruir el canvas y ausencia de mutaciones químicas.
- Renderer: ciclo de contexto, subida diferida durante `paintGL`, teardown y reconstrucción de escena/buffers al recrear contexto.
- Geometría 3D: hidrógenos y enlaces generados aparecen en escena sin mutar grafo ni persistir en `.cmsn`.
- Arquitectura: catálogo M08, imports de `geometry3d` sin Qt inverso y ausencia de Avogadro/GL_SELECT/PyOpenGL en el nuevo subpaquete.
- Verificación: compileall, tests focalizados, Ruff selectivo, OpenSpec estricto y `git diff --check`. La suite monolítica seguirá documentando cualquier fallo baseline no causado por este cambio.
