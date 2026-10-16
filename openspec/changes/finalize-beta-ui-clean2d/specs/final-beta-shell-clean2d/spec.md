# Spec Delta — cabecera compacta y Clean2D estereoquímicamente seguro

## MODIFIED Requirements

### Requirement: La cabecera beta tiene dos filas compactas

La superficie superior SHALL mostrar por defecto la barra de menús tradicional con Archivo, Editar, Ver, Estructura, Reacción y Ayuda, seguida por una única franja compacta de documentos/acciones. La AppBar SHALL NOT mostrar el nombre grande ChemUSON, un matraz ni la etiqueta de versión; esos datos SHALL permanecer en metadatos, título principal y Acerca de donde corresponda. La franja compacta SHALL exponer Nuevo, Abrir, Guardar y Pegar mediante las QAction existentes; Pegar SHALL aparecer una sola vez en el menú Editar y ningún comando/handler/shortcut SHALL duplicarse. Deshacer/Rehacer SHALL seguir sólo a la derecha de la AppBar. La barra superior adicional de acciones SHALL NOT aparecer en defaults nuevos ni en la migración de preferencias anteriores.

#### Scenario: Inicio con preferencias anteriores o nuevas
- **WHEN** la aplicación abre una configuración nueva o antigua
- **THEN** se ven primero los seis menús y debajo una sola franja compacta de documentos/acciones
- **AND** no aparecen la marca grande, matraz, versión ni una tercera fila de toolbar
- **AND** el título, Acerca de y metadatos conservan ChemUSON y la versión oficial.

#### Scenario: Acciones compartidas
- **WHEN** se pulsa Nuevo, Abrir, Guardar o Pegar desde la franja compacta o desde su menú correspondiente
- **THEN** se activa la misma QAction/controlador, sin acción duplicada
- **AND** Ctrl+S guarda el documento activo y Undo/Redo no aparecen dos veces en la cabecera.

#### Scenario: Cambios de tema, documento y panel
- **WHEN** se alterna tema claro/oscuro, se modifica una pestaña o se redimensiona la ventana
- **THEN** los cuatro controles permanecen alineados/accesibles, el punto dirty no desplaza la X, y Save/Discard/Cancel conserva el contenido
- **AND** el panel lateral inicia oculto y sigue redimensionable de forma independiente.

### Requirement: Clean2D cuenta centros estereoquímicos distintos

La clasificación de complejidad SHALL contar cada centro estereoquímico estructural una sola vez por ID de átomo. Un enlace dirigido wedge/hash SHALL identificar el centro en su extremo de origen, no crear un centro separado en ambos extremos. Varios enlaces wedge/hash que compartan origen SHALL contribuir a un solo centro. La regla existente `multiple_stereo_centers → preserve_only` SHALL seguir aplicándose a dos o más centros distintos. Una cuña geométrica aislada SHALL NOT bastar para afirmar una asignación CIP química. Para un centro wedge/hash con asignación R/S independiente, el controlador GUI SHALL verificar la configuración del dibujo antes/después en el worker RDKit aislado, omitiendo CIP cacheado del grafo; SHALL rechazar la aplicación si la asignación cambia o deja de estar disponible tras un baseline asignado.

#### Scenario: Un centro real con una sola cuña
- **GIVEN** `C[C@H](F)CO` con un centro CIP independiente y una cuña/hash importada
- **WHEN** se clasifica Clean2D sobre la molécula completa
- **THEN** se informa un centro estereoquímico, no `multiple_stereo_centers`
- **AND** candidatos reales se evalúan con los filtros existentes; si ninguno es seguro se preserva el dibujo.

#### Scenario: Dos cuñas sobre un solo centro tetraédrico
- **GIVEN** dos enlaces wedge/hash con el mismo átomo de origen
- **WHEN** se clasifica la molécula
- **THEN** ambos se asocian a un único centro
- **AND** el estilo, configuración y enlaces permanecen intactos al limpiar.

#### Scenario: Dos centros quirales distintos
- **GIVEN** una molécula con dos átomos estereoquímicos diferentes
- **WHEN** se calcula el perfil
- **THEN** se cuentan dos centros y la política preexistente puede preservar la geometría global.

#### Scenario: Una cuña no asigna quiralidad por sí sola
- **GIVEN** una estructura dibujada con cuña cuyos sustituyentes son equivalentes o cuya configuración química está sin especificar
- **WHEN** se asigna estereoquímica química independiente
- **THEN** no se exige un valor R/S inventado, aunque la firma geométrica siga protegiendo estilo y orden de la cuña.

#### Scenario: Selección parcial y multicomponente
- **GIVEN** varias moléculas desconectadas y una selección limitada a una molécula, o una selección parcial que incluye un centro
- **WHEN** se ejecuta Quick, publicación o propuesta
- **THEN** sólo pueden cambiar los átomos seleccionados; los vecinos de frontera participan en la protección geométrica
- **AND** no se altera conectividad, identidad, cuñas ni configuración; si no hay candidato seguro se conserva la geometría con motivo comprensible.

#### Scenario: Una propuesta GUI no puede voltear R/S por geometría implícita
- **GIVEN** un centro con cuña/hash y una asignación R/S detectable de forma independiente
- **WHEN** Clean2D Quick/Publicación intenta un fallback automático de propuesta o se aplica una propuesta
- **THEN** el controlador obtiene la asignación desde la cuña/hash sin confiar en CIP cacheado
- **AND** rechaza cualquier candidato que cambie o pierda esa asignación.

#### Scenario: Round-trip de enantiómeros
- **GIVEN** los enantiómeros de referencia `C[C@H](F)CO` y `C[C@@H](F)CO`
- **WHEN** se calcula Clean2D en Quick, Publicación y propuesta, y se exporta/reimporta por el camino de prueba soportado
- **THEN** la asignación independiente R/S de cada uno se conserva y ambos permanecen distintos.
