## ADDED Requirements

### Requirement: La geometría 3D es estado editable propio de cada documento

Cada `ChemusonCanvas` SHALL poseer una geometría 3D identificada por IDs estables de `MolGraph`. La geometría SHALL ser independiente de la cámara, de la selección y de `Atom.x/y`. El viewport activo SHALL mostrar la geometría del canvas activo y SHALL NOT transferirla ni descartarla al cambiar pestañas.

#### Scenario: Dos documentos conservan conformaciones distintas
- **GIVEN** dos documentos con el mismo grafo químico pero coordenadas 3D diferentes
- **WHEN** el usuario alterna entre sus pestañas
- **THEN** cada escena muestra su propio `CoordinateSet3D`
- **AND** las coordenadas 2D y cámaras permanecen independientes

#### Scenario: Abrir documento sin geometría 3D
- **GIVEN** un `.cmsn` anterior que no contiene `canvas.geometry3d`
- **WHEN** se carga el documento
- **THEN** la química 2D se carga normalmente y el dock muestra que falta generar geometría
- **AND** no se inventa una conformación ni se informa éxito de optimización

### Requirement: La geometría 3D tiene persistencia opcional compatible

Al guardar un documento que tenga posiciones 3D explícitas, ChemUSON SHALL serializar `canvas.geometry3d` con `schema_version: 1`, metadatos de método/campo de fuerza/energía y coordenadas finitas en Å indexadas por `atom_id`. SHALL guardar sólo IDs existentes en `MolGraph`. SHALL aceptar archivos sin el campo y SHALL recuperar de manera controlada frente a entradas 3D inválidas, sin perder el modelo 2D. Cámara, selección, trabajo activo, frames parciales e hidrógenos generados SHALL NOT persistirse.

#### Scenario: Round-trip de coordenadas
- **GIVEN** una geometría 3D válida para átomos documentales
- **WHEN** se guarda y vuelve a cargar el `.cmsn`
- **THEN** cada `atom_id` conserva su tripleta XYZ, fuente, método y energía
- **AND** `Atom.x/y`, conectividad, estereoquímica y datos de canvas 2D permanecen iguales

#### Scenario: Entrada 3D malformada
- **GIVEN** coordenadas con NaN/inf, IDs duplicados/desconocidos, dimensión incorrecta o versión no soportada
- **WHEN** se carga el `.cmsn`
- **THEN** la geometría inválida se descarta con diagnóstico estable
- **AND** el documento químico 2D y sus anotaciones se cargan

#### Scenario: Guardar tras borrar un átomo
- **GIVEN** un modelo con coordenadas conservadas temporalmente para undo de un átomo eliminado
- **WHEN** el documento se guarda
- **THEN** `canvas.geometry3d.positions` sólo contiene IDs presentes en el `MolGraph` actual

### Requirement: La edición 3D requiere activación explícita y no cambia la cámara

El dock SHALL ofrecer un modo de edición explícito, desactivado por defecto. Con el modo desactivado SHALL mantener navegación y picking actuales. En modo activo, los gestos de geometría SHALL modificar coordenadas/topología según herramienta sin cambiar la orientación, el centro ni la escala de cámara. `Fit` y cambio de proyección SHALL modificar sólo la cámara.

#### Scenario: Modo read-only inicial
- **GIVEN** un viewport recién abierto
- **WHEN** el usuario arrastra con los gestos de navegación existentes
- **THEN** sólo cambia la cámara y la geometría/documento permanecen intactos

#### Scenario: Arrastre de átomo
- **GIVEN** una geometría válida, modo edición activo y un átomo seleccionado
- **WHEN** el usuario arrastra el átomo
- **THEN** cambia la posición 3D del ID elegido sobre un plano paralelo a cámara a la profundidad inicial
- **AND** no cambian su ID, `Atom.x/y`, conectividad, cámara ni el resto de coordenadas

#### Scenario: Traslación de componente conectado
- **GIVEN** un componente seleccionado mediante doble clic
- **WHEN** el usuario lo arrastra como fragmento
- **THEN** se aplica el mismo vector 3D a todos sus átomos
- **AND** las distancias internas y la topología del componente se conservan

#### Scenario: Fit no es edición química
- **GIVEN** una geometría editada
- **WHEN** el usuario activa `Fit` o cambia perspectiva/ortográfica
- **THEN** sólo cambian parámetros de cámara y el estado undo/redo no cambia

### Requirement: Distancias, ángulos y dihedros se editan con ratón y valores numéricos

El editor SHALL ofrecer controles numéricos con unidades explícitas —distancia en Å y ángulos en grados— para el enlace o conjunto de átomos seleccionado. SHALL ofrecer gestos de ratón equivalentes para estirar un enlace, variar un ángulo o cambiar un dihedro válido. Las transformaciones SHALL conservar IDs/conectividad, rechazar valores no finitos y mostrar errores de geometría degenerada sin mutar el estado.

#### Scenario: Fijar una distancia de enlace
- **GIVEN** un enlace covalente editable seleccionado y un valor de distancia positivo y finito
- **WHEN** el usuario confirma el valor numérico o el gesto de extensión
- **THEN** la distancia geométrica entre sus extremos coincide con el valor dentro de la tolerancia definida
- **AND** el enlace sigue usando los mismos IDs y propiedades químicas

#### Scenario: Fijar un ángulo
- **GIVEN** un átomo central y dos vecinos con geometría no degenerada
- **WHEN** el usuario cambia el valor numérico o arrastra el ángulo
- **THEN** el ángulo geométrico resultante coincide con el objetivo
- **AND** sólo rota el lado seleccionado como unidad alrededor del centro

#### Scenario: Fijar un dihedro acíclico
- **GIVEN** cuatro átomos consecutivos cuyo enlace central es simple, covalente, no aromático y puente del grafo
- **WHEN** se aplica un objetivo de dihedro
- **THEN** se rota un lado completo alrededor del eje central sin cambiar longitudes internas
- **AND** se conserva conectividad y los átomos fuera del lado móvil

#### Scenario: Dihedro sobre enlace cíclico
- **GIVEN** un enlace central cuyos extremos siguen conectados al retirar temporalmente ese enlace
- **WHEN** se intenta arrastrar o fijar un dihedro independiente
- **THEN** la acción se rechaza como no rotatable por cierre de anillo
- **AND** no se cambian coordenadas ni topología

### Requirement: El constructor 3D modifica el único `MolGraph`

La herramienta SHALL permitir seleccionar un elemento compatible y crear átomos libres o enlazados, crear/modificar un enlace compatible y conservar los IDs estables. Las operaciones SHALL usar `MolGraph`, validación existente y los comandos existentes/compuestos de canvas; SHALL rechazar endpoint inexistente, enlace duplicado, autoenlace y orden/valencia no admitidos. H generados por RDKit SHALL permanecer derivados salvo que se creen explícitamente.

#### Scenario: Añadir átomo enlazado
- **GIVEN** un átomo ancla válido, elemento seleccionable y geometría disponible
- **WHEN** el usuario crea un átomo enlazado
- **THEN** se crea un único átomo documental y un único enlace en el `MolGraph`
- **AND** se asignan posiciones 2D/3D coherentes y una sola acción de undo
- **AND** el renderer muestra la conectividad nueva con los H visuales recalculados

#### Scenario: Cambiar orden de enlace
- **GIVEN** un enlace no-query editable
- **WHEN** se elige un orden permitido que supera la validación química existente
- **THEN** el `MolGraph` se actualiza y la escena 3D refleja el cambio
- **AND** no se crean enlaces paralelos

#### Scenario: Cambio estructural iniciado desde 2D
- **GIVEN** coordenadas 3D guardadas y una operación 2D de añadir, borrar o cambiar conectividad
- **WHEN** la operación se confirma, deshace o rehace
- **THEN** la escena se reconcilia con el `MolGraph` actualizado
- **AND** posiciones de IDs supervivientes se conservan
- **AND** H generados se invalidan y recalculan para la nueva conectividad

### Requirement: Cada gesto 3D es una transacción undoable

Una edición geométrica confirmada SHALL producir exactamente una operación del `QUndoStack` del canvas activo, con snapshots antes/después. Una edición topológica más su geometría SHALL ser un comando compuesto/macro que se deshace como una acción lógica. Undo/redo SHALL restaurar una combinación consistente de topología, coordenadas, selección documental y geometría renderizada. Movimiento de cámara, selección, frames de preview y eventos intermedios del puntero SHALL NOT crear comandos.

#### Scenario: Undo/redo de arrastre de átomo
- **GIVEN** un átomo movido y gesto confirmado
- **WHEN** se ejecuta undo y luego redo
- **THEN** las coordenadas 3D vuelven exactamente al snapshot anterior y al final
- **AND** la coordenada 2D y el grafo no varían

#### Scenario: Undo/redo de construcción
- **GIVEN** un átomo y enlace creados por una acción 3D
- **WHEN** se ejecuta undo y redo
- **THEN** undo elimina átomos/enlaces y posiciones creadas en una sola acción
- **AND** redo restaura IDs y geometría original

#### Scenario: Cancelar un gesto
- **GIVEN** un arrastre 3D no confirmado
- **WHEN** el usuario cancela con Escape
- **THEN** el snapshot inicial se restaura y no se añade entrada al historial
