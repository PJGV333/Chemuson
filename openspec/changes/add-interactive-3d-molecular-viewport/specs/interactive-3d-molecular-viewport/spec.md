## ADDED Requirements

### Requirement: El visor consume una escena 3D derivada sin mutar el documento

El visor SHALL aceptar `SceneMolecule3D`/`CoordinateSet3D` existentes y SHALL mantener la escena, la cámara y la selección visual como estado propio. La apertura, navegación, renderizado y selección SHALL NOT modificar `MolGraph`, `ChemState`, la selección del canvas 2D, la pila undo/redo ni el formato `.cmsn`.

#### Scenario: Abrir y cerrar el visor conserva el documento
- **GIVEN** un documento con átomos, enlaces, coordenadas 2D, selección y estado undo observables
- **WHEN** se abre, navega, selecciona y cierra el visor 3D
- **THEN** los átomos, enlaces, coordenadas, selección documental, índice undo y estado dirty son iguales antes y después

#### Scenario: Escena incompleta
- **GIVEN** una escena con coordenadas ausentes para parte de los átomos
- **WHEN** el visor recibe la escena
- **THEN** omite sólo primitivas no representables, conserva los IDs disponibles y no muta el modelo químico

### Requirement: La cámara soporta navegación y proyecciones deterministas

La cámara SHALL exponer orientación orbital, centro de rotación, distancia, perspectiva y ortográfica. Orbit, pan, zoom y fit SHALL estar acotados, ser finitos y mantener el centro de rotación explícito. Perspective y orthographic SHALL producir rayos de picking compatibles con la misma convención de viewport.

#### Scenario: Orbit/pan/zoom
- **GIVEN** una cámara válida y un viewport positivo
- **WHEN** se aplican orbit, pan y zoom con deltas finitos
- **THEN** la cámara conserva parámetros finitos, respeta los límites de distancia/escala y no altera la escena

#### Scenario: Encuadre automático
- **GIVEN** una escena con una o más primitivas válidas
- **WHEN** se solicita fit
- **THEN** el centro se aproxima al centro de sus bounds, la distancia/escala es positiva y la escena queda dentro del volumen de vista con padding

#### Scenario: Viewport degenerado
- **GIVEN** un ancho o alto cero, negativo o no finito
- **WHEN** se solicita un rayo o una transformación de viewport
- **THEN** el contrato devuelve un resultado controlado o una excepción de validación estable, nunca NaN silencioso

### Requirement: El renderizado usa OpenGL moderno de PyQt6

El renderer SHALL usar `QOpenGLWidget`, shaders y buffers Qt (`QOpenGLShaderProgram`, VBO/VAO y funciones OpenGL core cuando estén disponibles). SHALL renderizar con prueba de profundidad, átomos como esferas y enlaces como varillas. SHALL NOT usar `GL_SELECT`, display lists, renderizado inmediato legado, PyOpenGL, Qt4 ni código ejecutable de Avogadro.

#### Scenario: Renderizado de escena válida
- **GIVEN** un contexto OpenGL inicializado y una escena con átomos/enlaces
- **WHEN** se pinta el viewport
- **THEN** se emiten primitivas con colores por elemento, profundidad habilitada y geometría derivada de la escena

#### Scenario: Contexto no disponible
- **GIVEN** un host sin contexto OpenGL utilizable o con shader incompatible
- **WHEN** el viewport intenta inicializarse
- **THEN** el dock muestra un estado de indisponibilidad controlado
- **AND** el editor 2D, sus documentos y sus acciones permanecen utilizables

### Requirement: Picking CPU identifica la primitiva visible más cercana

El picking SHALL construir un rayo desde la cámara y probar átomos como esferas y enlaces como cápsulas finitas. SHALL devolver como máximo una selección con tipo (`atom`/`bond`), ID y distancia positiva; SHALL ignorar geometría no finita o degenerada y SHALL escoger la menor distancia válida.

#### Scenario: Seleccionar átomo
- **GIVEN** un átomo visible y un punto de cursor que atraviesa su esfera
- **WHEN** se solicita picking
- **THEN** el resultado identifica `atom` y el ID exacto del átomo

#### Scenario: Seleccionar enlace
- **GIVEN** un enlace visible entre dos átomos y un cursor que atraviesa su cápsula
- **WHEN** se solicita picking sin impactar una esfera más cercana
- **THEN** el resultado identifica `bond` y el ID exacto del enlace

#### Scenario: Sin impacto o empate
- **GIVEN** un cursor que no interseca primitivas o primitivas a igual distancia
- **WHEN** se solicita picking
- **THEN** se devuelve `None` para ausencia y se aplica la prioridad determinista átomo antes que enlace para empate

### Requirement: La selección del visor es local y observable

La selección visual SHALL vivir dentro del widget 3D, SHALL resaltarse sin cambiar el modelo y SHALL emitir una señal pública con el tipo e ID seleccionado. Un clic vacío SHALL limpiar la selección local.

#### Scenario: Selección local
- **GIVEN** una escena renderizada y un átomo o enlace elegible
- **WHEN** el usuario hace clic sobre la primitiva
- **THEN** el visor resalta esa primitiva y emite su tipo e ID
- **AND** `MolGraph.selected_atoms` y `MolGraph.selected_bonds` permanecen intactos

#### Scenario: Clic vacío
- **GIVEN** una selección local existente
- **WHEN** el usuario hace clic en espacio sin geometría
- **THEN** la selección local queda vacía y el documento 2D no cambia

### Requirement: Integración GUI discreta y reversible

Chemuson SHALL exponer un `QDockWidget` de visor 3D oculto por defecto y una única entrada `Ver → Visor molecular 3D` basada en su acción de visibilidad. La ventana SHALL reutilizar el `CoordinateSet3D` del flujo CompChem y SHALL actualizar la escena al terminar una generación válida. Cambiar de pestaña SHALL limpiar la escena efímera del documento anterior.

#### Scenario: Abrir desde Ver
- **GIVEN** una `ChemusonWindow` con documentos y el menú Ver construido
- **WHEN** se activa `Visor molecular 3D`
- **THEN** se muestra el dock existente sin crear otro canvas o modelo químico
- **AND** el viewer recibe la escena derivada o un estado informativo si aún no hay coordenadas

#### Scenario: Generación 3D termina con el visor abierto
- **GIVEN** el dock visible y un trabajo CompChem que produce coordenadas
- **WHEN** el controlador emite el resultado válido
- **THEN** el visor reconstruye `SceneMolecule3D` para el canvas activo sin mutar el documento

#### Scenario: Cambio de documento
- **GIVEN** el visor visible con una escena del documento A
- **WHEN** se activa el documento B
- **THEN** se descarta la escena/coordset efímera de A mediante la ruta existente de reset
- **AND** no se reutiliza la selección 3D de A en B
