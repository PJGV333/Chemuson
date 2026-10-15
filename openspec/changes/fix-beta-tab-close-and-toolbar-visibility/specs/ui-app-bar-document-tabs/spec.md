# Spec Delta — shell de la beta: app bar, menús, barra de acciones y pestañas

## MODIFIED Requirements

### Requirement: Barra de aplicación integrada en la ventana real

La aplicación SHALL exponer una `AppBar` (`src/chemuson/gui/app_bar.py`) de altura fija 54 px en la parte superior del área central (sobre el `QTabWidget` de documentos), con marca, `DocumentTabBar`, botón `+`, píldora de búsqueda y botones undo/redo/tema/ajustes existentes. El `QMenuBar` histórico SHALL permanecer presente y funcional, con los seis menús tradicionales (`Archivo`, `Editar`, `Ver`, `Estructura`, `Reacción`, `Ayuda`) **visibles por defecto** y en una fila compacta (altura ≤ 36 px), sin duplicar `QMenu`, `QAction`, atajos ni controladores. La barra clásica superior de acciones (`main_toolbar`) SHALL conservar su `QToolBar` y sus QActions originales, ser visible por defecto, contener exclusivamente las acciones de Nuevo, Abrir y Guardar, y quedar controlada por la preferencia persistente **Mostrar barra de herramientas**. La barra izquierda de dibujo SHALL continuar reemplazada visualmente por el rail y el panel derecho SHALL mantener su visibilidad independiente y oculta al iniciar.

#### Scenario: La ventana real monta las superficies
- **GIVEN** un `ChemusonWindow` instanciado y mostrado con preferencias recién creadas
- **WHEN** se inspecciona la jerarquía de widgets
- **THEN** existe `window.app_bar` visible de altura 54, y la tab bar nativa de `tabs` permanece oculta
- **AND** `main_toolbar` es visible por defecto y sus únicas acciones son las QActions originales de nuevo, abrir y guardar
- **AND** `menuBar()` está visible, contiene los seis menús originales y su altura no supera 36 px ni se solapa con la app bar
- **AND** la barra lateral de dibujo sigue en el rail y el panel derecho comienza oculto e independiente de ambas barras.

#### Scenario: El modo hamburguesa conserva los menús originales
- **GIVEN** la ventana con `QMenuBar` oculto desde Ver
- **WHEN** se pulsa la hamburguesa/Alt o se abre Ver
- **THEN** las mismas seis QMenu históricas y sus QActions funcionan y las opciones Mostrar barra de menús y Mostrar barra de herramientas aparecen agrupadas al principio de Ver.

## ADDED Requirements

### Requirement: La barra de menús tradicional es visible por defecto y compacta

La barra de menús tradicional SHALL mostrarse por defecto con los seis menús existentes y SHALL poder ocultarse desde `Ver → Mostrar barra de menús`, cuya preferencia persistida SHALL seguir gobernando la visibilidad. La fila de menús SHALL ser compacta (altura ≤ 36 px) y legible en temas claro y oscuro, ventanas pequeñas y DPI elevado. Al ocultarla, el acceso hamburguesa/Alt SHALL seguir disponible sobre los mismos `QMenu`. El default nuevo SHALL migrar una sola vez sobre preferencias escritas antes del cambio, mediante un marcador de versión de defaults, y toda elección posterior del usuario SHALL respetarse.

#### Scenario: Arranque con preferencias nuevas
- **WHEN** se construye la ventana sin preferencias de shell persistidas
- **THEN** la barra de menús es visible, sus seis ítems tienen geometría válida y su altura es ≤ 36 px
- **AND** la app bar de 54 px comienza debajo de la barra de menús sin solapamiento.

#### Scenario: Ocultar la barra conservando el acceso a los menús
- **GIVEN** la barra de menús visible
- **WHEN** el usuario desactiva Mostrar barra de menús
- **THEN** la barra se oculta y el botón hamburguesa de la app bar sigue visible y abre el agregador con exactamente los mismos seis `QMenu`
- **AND** la preferencia queda persistida y se reaplica al reiniciar.

#### Scenario: Preferencias previas al cambio de default
- **GIVEN** un archivo de configuración con `ui/show_menu_bar=false` escrito antes de este cambio y sin marcador de defaults
- **WHEN** se abre la ventana
- **THEN** se aplican los defaults nuevos de la beta (barra de menús y barra de acciones visibles)
- **AND** tras la primera elección explícita del usuario esa elección persiste y no vuelve a migrarse.

### Requirement: La barra de herramientas superior contiene sólo Nuevo, Abrir y Guardar

`main_toolbar` SHALL añadir únicamente las QActions existentes de Nuevo, Abrir y Guardar, con iconos theme-aware de 18 px y botones de tamaño uniforme en una fila de altura ≤ 40 px. Las acciones de Deshacer, Rehacer, rotaciones, volteos, Clean2D y SMILES NO SHALLN estar en esa barra, y su única representación visible de Undo/Redo SHALLN ser los botones de la app bar; cada acción retirada SHALL seguir alcanzable por sus rutas actuales (menús, atajos, rail y paleta de comandos) sin crear QActions, señales ni atajos duplicados.

#### Scenario: Contenido y uniformidad de la barra
- **WHEN** se inspecciona `main_toolbar.actions()` en tema claro u oscuro y en ventanas de 1440×900, 980×600 o 900×560
- **THEN** las acciones no separador son exactamente Nuevo, Abrir y Guardar (las mismas instancias que en Archivo)
- **AND** los tres botones tienen el mismo tamaño, icono no vacío y caben dentro de la barra
- **AND** undo, redo, rotaciones, volteos, Clean2D y SMILES no están en la barra pero sí en sus menús/atajos.

#### Scenario: La opción de Ver controla sólo esa barra
- **GIVEN** la ventana con barra de herramientas visible
- **WHEN** se desactiva Mostrar barra de herramientas
- **THEN** sólo se oculta `main_toolbar`; la barra de menús, la app bar, el rail y el panel lateral no cambian
- **AND** los atajos de Guardar, Undo y Redo siguen alcanzando sus QActions originales.

### Requirement: Los iconos y el área de cierre de pestaña permanecen visibles

Cada pestaña de `DocumentTabBar` SHALL inicializar al construirla un icono X no vacío de `IconProvider`, visible sin hover, con contraste del token de icono del tema, área de interacción de al menos 24×24 px, tooltip y nombre accesible “Cerrar pestaña”. El icono de cierre SHALL mantenerse separado del punto dirty de 7 px, actualizarse al cambiar de tema y seguir unido a su pestaña lógica al añadir/reordenar/seleccionar documentos. La activación SHALL emitir el índice actual y delegar al flujo de confirmación/cierre existente; no SHALL modificar `CanvasTabManager` ni crear un segundo flujo.

#### Scenario: Pestaña nueva y tema
- **GIVEN** una ventana en tema claro u oscuro
- **WHEN** se crea o sincroniza una pestaña limpia o dirty
- **THEN** su X se renderiza visiblemente antes de hover con el tinte correcto
- **AND** el dirty dot refleja sólo el estado clean/dirty y no oculta ni sustituye la X.

#### Scenario: Reordenar y cerrar la pestaña objetivo
- **GIVEN** varias pestañas, algunas reordenadas
- **WHEN** se hace clic en la X de una pestaña
- **THEN** se cierra sólo el documento debajo de ese botón y el flujo Save/Discard/Cancel existente conserva su comportamiento.

### Requirement: La geometría de la pestaña no cambia entre limpio y modificado

El bloque lateral de cada pestaña SHALL tener tamaño fijo con el espacio del indicador dirty reservado de forma permanente; pasar de documento limpio a modificado (y viceversa) SHALL alterar únicamente la pintura del indicador. El botón de cierre SHALL conservar su tamaño, su posición dentro del bloque y su holgura respecto al borde derecho de la pestaña, y su rect completo SHALL quedar dentro del `tabRect`, sin solaparse con el botón `+` ni con otras pestañas.

#### Scenario: Mismo título, sólo cambia el estado dirty
- **GIVEN** una pestaña con título fijo y documento limpio
- **WHEN** el documento pasa a modificado y luego vuelve a limpio
- **THEN** el `tabRect`, el tamaño del bloque lateral, el tamaño de la X, su posición dentro del bloque y su holgura al borde derecho son idénticos en los tres estados
- **AND** el punto dirty se pinta con el color de acento sólo mientras el documento está modificado.

#### Scenario: Temas, tamaños y varias pestañas
- **WHEN** se cambia de tema, se redimensiona la ventana (1440×900, 980×600, 900×560), se añaden pestañas con títulos largos o se reordenan
- **THEN** la geometría del bloque lateral y de la X permanece invariante y dentro de su pestaña
- **AND** el clic en la X cierra el documento correspondiente con las confirmaciones existentes.

### Requirement: La visibilidad de barras persiste independientemente

Ver SHALL presentar al inicio del menú las acciones checkable distintas **Mostrar barra de menús** y **Mostrar barra de herramientas**. La segunda SHALL controlar la visibilidad del `QToolBar main_toolbar` existente, cuya preferencia `ui/show_main_toolbar` SHALL ser booleana y por defecto true. `ui/show_menu_bar` SHALL conservar su significado con default true en la beta y migrarse una sola vez desde preferencias anteriores mediante `ui/shell_defaults_version`. Ambas preferencias SHALL persistir independientemente de `ui/side_panel/visible`; la opción auxiliar copiar/pegar/zoom SHALL seguir controlando sólo la presencia de esas acciones auxiliares dentro de la barra principal. Cambiar tema o una preferencia no SHALL resetear las otras visibilidades.

#### Scenario: Primer arranque y persistencia
- **WHEN** no hay valores `ui/show_main_toolbar` ni `ui/show_menu_bar` y se construye la ventana
- **THEN** ambas acciones de Ver están marcadas y `main_toolbar` y la barra de menús están visibles
- **WHEN** el usuario alterna Mostrar barra de herramientas y reinicia
- **THEN** la misma barra existente reabre con la visibilidad elegida.

#### Scenario: Independencia de panel lateral y menú
- **GIVEN** el panel lateral derecho está oculto
- **WHEN** se muestra u oculta la barra de herramientas
- **THEN** sólo cambia `main_toolbar`, no el rail de dibujo, el panel lateral ni la visibilidad del menú
- **AND** alternar la barra de menús no cambia la preferencia de toolbar.

#### Scenario: Ventana pequeña y tema
- **WHEN** se usa la ventana a 980×600 y 900×560 y se cambia entre tema claro/oscuro
- **THEN** la barra de acciones cabe o usa el overflow normal de QToolBar, el AppBar/rail/canvas conservan geometría válida y el panel lateral no se abre automáticamente.

### Requirement: Atajos existentes siguen alcanzables desde las barras

Las QActions compartidas de Guardar, Undo y Redo SHALL conservar sus atajos/controller existentes al mostrarse/ocultarse la barra. Ninguna acción, conexión o shortcut se SHALL duplicar; Ctrl+S SHALL seguir llamando al flujo de guardado CMSN ya aceptado.

#### Scenario: Guardar, deshacer y rehacer por teclado
- **WHEN** se invocan los atajos Ctrl+S y los atajos estándar de Undo/Redo en una ventana con toolbar visible u oculta
- **THEN** cada atajo llega a su QAction/controller original y el estado del documento refleja Undo/Redo correctamente.
