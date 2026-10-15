## MODIFIED Requirements

### Requirement: Barra de aplicación integrada en la ventana real

La aplicación SHALL exponer una `AppBar` (`src/chemuson/gui/app_bar.py`) de altura fija 54 px en la parte superior del área central (sobre el `QTabWidget` de documentos), con marca, `DocumentTabBar`, botón `+`, píldora de búsqueda y botones undo/redo/tema/ajustes existentes. El `QMenuBar` histórico SHALL permanecer presente y funcional, con visibilidad opcional. La barra clásica superior de acciones (`main_toolbar`) SHALL conservar su `QToolBar` y QActions originales, ser visible por defecto y quedar controlada por la preferencia persistente **Mostrar barra de herramientas**. La barra izquierda de dibujo SHALL continuar reemplazada visualmente por el rail y el panel derecho SHALL mantener su visibilidad independiente y oculta al iniciar.

#### Scenario: La ventana real monta las superficies
- **GIVEN** un `ChemusonWindow` instanciado y mostrado
- **WHEN** se inspecciona la jerarquía de widgets
- **THEN** existe `window.app_bar` visible de altura 54, y la tab bar nativa de `tabs` permanece oculta
- **AND** `main_toolbar` es visible por defecto y contiene las QActions originales de nuevo, abrir, guardar, deshacer y rehacer
- **AND** la barra de menús puede seguir compacta/oculta y la barra lateral de dibujo sigue en el rail
- **AND** el panel derecho comienza oculto e independiente de ambas barras.

#### Scenario: El modo hamburguesa conserva los menús originales
- **GIVEN** la ventana con `QMenuBar` oculto
- **WHEN** se pulsa la hamburguesa/Alt o se abre Ver
- **THEN** las mismas seis QMenu históricas y QActions funcionan y las opciones Mostrar barra de menús y Mostrar barra de herramientas aparecen agrupadas al principio de Ver.

### Requirement: Los iconos y el área de cierre de pestaña permanecen visibles

Cada pestaña de `DocumentTabBar` SHALL inicializar al construirse un icono X no vacío de `IconProvider`, visible sin hover, con contraste del token de icono del tema, área de interacción de al menos 24×24 px, tooltip y nombre accesible “Cerrar pestaña”. El icono de cierre SHALL mantenerse separado del punto dirty de 7 px, actualizarse al cambiar de tema y seguir unido a su pestaña lógica al añadir/reordenar/seleccionar documentos. La activación SHALL emitir el índice actual y delegar al flujo de confirmación/cierre existente; no SHALL modificar `CanvasTabManager` ni crear un segundo flujo.

#### Scenario: Pestaña nueva y tema
- **GIVEN** una ventana en tema claro u oscuro
- **WHEN** se crea o sincroniza una pestaña limpia o dirty
- **THEN** su X se renderiza visiblemente antes de hover con el tinte correcto
- **AND** el dirty dot refleja sólo el estado clean/dirty y no oculta ni sustituye la X.

#### Scenario: Reordenar y cerrar la pestaña objetivo
- **GIVEN** varias pestañas, algunas reordenadas
- **WHEN** se hace clic en la X de una pestaña
- **THEN** se cierra sólo el documento debajo de ese botón y el flujo Save/Discard/Cancel existente conserva su comportamiento.

### Requirement: La visibilidad de barras persiste independientemente

Ver SHALL presentar al inicio del menú las acciones checkable distintas **Mostrar barra de menús** y **Mostrar barra de herramientas**. La segunda SHALL controlar la visibilidad del `QToolBar main_toolbar` existente, cuya preferencia `ui/show_main_toolbar` SHALL ser booleana y por defecto true. La preferencia `ui/show_menu_bar` conserva su significado y su default compacto. Ambas preferencias SHALL persistir independientemente de `ui/side_panel/visible`; la opción auxiliar copiar/pegar/zoom SHALL seguir controlando sólo la presencia de esas acciones auxiliares dentro de la barra principal. Cambiar tema o una preferencia no SHALL resetear las otras visibilidades.

#### Scenario: Primer arranque y persistencia
- **WHEN** no hay valor `ui/show_main_toolbar` y se construye la ventana
- **THEN** el QAction de Ver está marcado y `main_toolbar` está visible
- **WHEN** el usuario alterna Mostrar barra de herramientas y reinicia
- **THEN** la misma barra existente reabre con la visibilidad elegida.

#### Scenario: Independencia de panel lateral y menú
- **GIVEN** el panel lateral derecho está oculto y el menú clásico está compacto
- **WHEN** se muestra u oculta la barra de herramientas
- **THEN** sólo cambia `main_toolbar`, no el rail de dibujo, el panel lateral ni la visibilidad del menú
- **AND** alternar la barra de menús no cambia la preferencia de toolbar.

#### Scenario: Ventana pequeña y tema
- **WHEN** se usa la ventana a 980×600 y 900×560 y se cambia entre tema claro/oscuro
- **THEN** la barra de acciones cabe o usa el overflow normal de QToolBar, el AppBar/rail/canvas conservan geometría válida y el panel lateral no se abre automáticamente.

## ADDED Requirements

### Requirement: Atajos existentes siguen alcanzables desde las barras

Las QActions compartidas de Guardar, Undo y Redo SHALL conservar sus atajos/controller existentes al mostrarse/ocultarse la barra. Ninguna acción, conexión o shortcut se SHALL duplicar; Ctrl+S SHALL seguir llamando al flujo de guardado CMSN ya aceptado.

#### Scenario: Guardar, deshacer y rehacer por teclado
- **WHEN** se invocan los atajos Ctrl+S y los atajos estándar de Undo/Redo en una ventana con toolbar visible u oculta
- **THEN** cada atajo llega a su QAction/controller original y el estado del documento refleja Undo/Redo correctamente.
