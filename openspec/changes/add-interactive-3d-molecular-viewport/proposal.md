# Proposal — visor molecular 3D interactivo

## Por qué

ChemUSON ya puede generar conformeros 3D y exponerlos como `SceneMolecule3D`, pero la GUI sólo ofrece controles de generación/proyección y no una superficie nativa para inspeccionar la geometría espacial. El primer visor debe aportar navegación, representación y picking sin convertir la fase inicial en un editor químico 3D ni acoplar `geometry3d` a Qt.

La referencia funcional de Avogadro 1.2.0 se usa únicamente para identificar capacidades esperadas. No se porta código, API Qt4, `GL_SELECT`, display lists ni código GPL.

## Alcance

- Crear una especificación y una implementación aislada bajo `src/chemuson/gui/viewport3d/`.
- Consumir `SceneMolecule3D` y `CoordinateSet3D` sin cambiar `MolGraph` ni `.cmsn`; completar sólo metadatos 3D efímeros de hidrógenos añadidos por RDKit para poder construir la escena completa.
- Añadir cámara orbital con perspectiva/ortográfica, pan, zoom y encuadre automático.
- Renderizar átomos como esferas y enlaces como varillas mediante `QOpenGLWidget` y shaders/VBO propios.
- Implementar picking CPU por rayos sobre esferas y cápsulas de enlace, con selección local del visor.
- Endurecer ciclo de vida OpenGL, encuadre respecto al aspect ratio, intersección finita de cápsulas y clipping para escalas moleculares amplias.
- Exponer el visor como un `QDockWidget` independiente, accesible desde `Ver`, sin alterar el `SidePanel` de siete docks.
- Añadir pruebas puras de cámara, transformaciones, picking, invariantes y pruebas GUI acotadas que no requieran un contexto OpenGL real.

## Fuera de alcance

- Edición 3D de átomos, enlaces, longitudes, ángulos, torsiones o fragmentos.
- Modificación de `MolGraph` desde el visor o sincronización de selección 2D.
- Persistencia de cámara/selección y cambios al formato `.cmsn`.
- Portar código de Avogadro o introducir Qt4, Python 2, PyOpenGL u otra dependencia runtime.
- Cambiar Clean2D, ChemName, el canvas 2D estabilizado o la infraestructura de documentos.
- Prometer validación visual acelerada por GPU si el host no puede crear un contexto OpenGL.

## Módulos afectados

- M08 GUI: nuevo subpaquete `gui/viewport3d/`, composición de `ChemusonWindow`, menú `Ver` y tests GUI.
- M05 geometry3d y M04 chemio: completar metadata transitoria de átomos/enlaces hidrógeno del conformero RDKit para `SceneMolecule3D`; no modificar `MolGraph` ni persistencia.
- Catálogo arquitectónico: registrar la propiedad del nuevo subpaquete y sus tests dentro de M08.

## Dependencias y decisión inicial

`QOpenGLWidget` y `PyQt6.QtOpenGL` forman parte de PyQt6. El renderer usará `QOpenGLShaderProgram`, VBO/VAO y `QOpenGLFunctions_4_1_Core` cuando el contexto esté disponible. No se añade PyOpenGL: el cálculo de cámara, mallas y picking permanece en Python puro y el backend Qt se degrada con un estado visible si no puede inicializar OpenGL.
