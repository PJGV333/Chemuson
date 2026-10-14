# ui-side-panel-status-bar Specification

## Purpose
TBD - created by archiving change 2026-09-27-modernize-ui-side-panel-status-bar. Update Purpose after archive.
## Requirements
### Requirement: Una única superficie lateral reutiliza los docks funcionales existentes

La aplicación SHALL exponer `SidePanel` en
`src/chemuson/gui/side_panel.py`, con ancho inicial de 340 px y mínimo legible
300 px, integrado al lado del lienzo mediante un `QSplitter(Qt.Horizontal)` y
sin chrome de docking clásico. El usuario SHALL poder ajustar el divisor; el
lienzo ocupa el espacio restante y el rail no participa del divisor.
`SidePanel` SHALL presentar Inspector, Validación, Propiedades, Plantillas y
Apariencia como tabs principales, con Inspector activo por defecto. SHALL
mantener Espectroscopía y CompChem accesibles mediante un control `…` dentro
de la misma región lateral.

Las siete instancias QDockWidget históricas SHALL conservarse como única
fuente de contenido funcional, estado, señales, métodos y conexiones. El
SidePanel SHALL alojar esas mismas instancias como páginas de un
`QStackedWidget`; SHALL NOT duplicar sus widgets internos, datos o handlers;
ninguna página SHALL mostrarse como QDockWidget clásico flotante o acoplado.

#### Scenario: La ventana monta cinco tabs y reutiliza siete instancias
- **GIVEN** un `ChemusonWindow` creado en Qt offscreen sin preferencias del panel
- **WHEN** se inspecciona `window.side_panel`
- **THEN** el panel existe con ancho dentro del rango de 336–340 px
- **AND** Inspector es la página inicial
- **AND** los cinco tabs principales aparecen en el orden especificado
- **AND** las páginas del stack son exactamente las instancias de los siete atributos dock históricos
- **AND** el widget de contenido de cada dock sigue siendo el mismo objeto
- **AND** no se crea ni registra un dock funcional duplicado.

### Requirement: Los cinco tabs principales permanecen legibles y sin clipping

`SideTabRow` SHALL conservar los cinco textos completos como contenido de cada
botón, exponerlos como tooltip y nombre accesible aunque la vista se recorte, y
mantener fuente de al menos 10 px, padding horizontal real de al menos 3 px por
tab y separación de al menos 3 px entre tabs cuando el ancho inicial es 340 px.
SHALL conservar el acento/underline del
tab activo y mantener Espectroscopía y CompChem dentro de `…`. La aceptación
visual SHALL basarse en texto íntegro, separación y contención en la fila; no
SHALL exigir `horizontalScrollBar().maximum() == 0` como objetivo independiente.

#### Scenario: Tabs legibles en temas y tamaños solicitados
- **GIVEN** SidePanel visible con ancho entre 336 y 340 px
- **WHEN** la ventana se captura a 1440×900 y 980×600 en temas light y dark
- **THEN** los cinco labels principales están completos, visibles, separados y sin clipping
- **AND** cada tab conserva padding horizontal y el activo conserva el acento/underline
- **AND** Espectroscopía y CompChem siguen accesibles desde `…`.

#### Scenario: Tooltip de una pestaña recortada
- **GIVEN** SidePanel reducido por el usuario hasta un ancho permitido
- **WHEN** una etiqueta se recorta visualmente
- **THEN** el tooltip y el nombre accesible conservan la etiqueta completa
- **AND** la navegación y el overflow existentes no cambian.

#### Scenario: La página activa viene del overflow
- **GIVEN** SidePanel visible
- **WHEN** se elige Espectroscopía o CompChem desde `…`
- **THEN** se activa la página de la instancia histórica correspondiente dentro del mismo stack
- **AND** el control `…` representa visualmente que está activo un tab de overflow
- **AND** ningún QDockWidget clásico aparece flotante ni en otra región.

### Requirement: Las acciones Ver navegan al SidePanel en vez de alternar docks

El menú Ver SHALL exponer acciones Inspector, Validación, Propiedades
químicas, Plantillas, Apariencia, Espectroscopía y CompChem. El SidePanel SHALL
iniciar oculto en toda ventana nueva independientemente de `visible` guardado.
Cada acción SHALL delegar a `side_panel.show_page(key)`, SHALL hacer visible el panel si
estaba oculto y SHALL seleccionar la página pedida. Las acciones SHALL NOT
usar `toggleViewAction()` de los docks ni mostrar un dock clásico.
Una acción separada de visibilidad puede ocultar/restaurar el SidePanel.

#### Scenario: Ver selecciona cada una de las siete páginas
- **GIVEN** el menú Ver de la ventana
- **WHEN** se activa cualquiera de sus siete acciones de panel
- **THEN** SidePanel queda visible con la página solicitada activa
- **AND** el widget de página activa es el QDockWidget histórico correspondiente
- **AND** los siete docks conservan `isFloating() == False` y área de docking nula.

#### Scenario: Inicio oculto con la preferencia anterior abierta
- **GIVEN** settings con `visible=true` y `active_tab=compchem`
- **WHEN** se crea una ventana nueva
- **THEN** el panel permanece oculto, CompChem sigue seleccionada y el valor guardado no fuerza su apertura.

#### Scenario: Ocultar y restaurar el panel
- **GIVEN** SidePanel visible con una página seleccionada
- **WHEN** el usuario desactiva y vuelve a activar la acción Panel lateral
- **THEN** el panel se oculta y restaura sin recrear páginas ni perder la selección
- **AND** la preferencia `visible` refleja el estado elegido.

#### Scenario: La ruta de Apariencia del AppBar conserva los contratos existentes
- **GIVEN** un `ChemusonWindow` con el AppBar
- **WHEN** se usa la hamburguesa del AppBar y se elige Ver → Apariencia
- **THEN** SidePanel se muestra con Apariencia activa
- **AND** la QAction de preferencias del botón dedicado del AppBar sigue abriendo su PreferencesDialog histórico.

### Requirement: Estado del SidePanel persistido y normalizado

La aplicación SHALL persistir `ui/side_panel/active_tab` y
`ui/side_panel/visible` mediante `chemuson.platform.settings`, sin usar la
visibilidad persistida para abrir una ventana nueva. La pestaña activa SHALL
restaurarse, pero cada ventana inicia con el panel oculto. Claves de tab
desconocidas o valores de visibilidad inválidos SHALL degradar a `inspector` y
`true` respectivamente como valores normalizados. Cambiar tab o visibilidad SHALL guardar
el estado sin duplicar el estado funcional de los docks.

#### Scenario: Restaurar página sin restaurar visibilidad
- **GIVEN** settings con `active_tab=compchem` y cualquier valor de `visible`
- **WHEN** se crea SidePanel
- **THEN** CompChem es la página activa y el panel comienza oculto.

#### Scenario: Defaults y valores inválidos
- **GIVEN** settings ausentes, una clave de tab desconocida o una visibilidad inválida
- **WHEN** se cargan las preferencias
- **THEN** cada valor ausente/inválido toma su default seguro (`inspector`/`true`)
- **AND** la ventana inicia oculta aunque el valor normalizado de `visible` sea `true`.

#### Scenario: Persistir navegación
- **GIVEN** SidePanel creado con un store de settings
- **WHEN** se muestra Apariencia y luego se oculta el panel
- **THEN** el store contiene `active_tab=appearance` y `visible=false`.

### Requirement: Los contratos funcionales de los siete docks permanecen activos

Inspector SHALL seguir recibiendo la selección real; Validación SHALL seguir
mostrando issues y conservando selección/corrección cruzada; Propiedades
SHALL seguir actualizándose; Plantillas SHALL conservar señales y acciones;
Espectroscopía SHALL conservar selección cruzada; CompChem SHALL conservar sus
acciones; Apariencia SHALL seguir modificando el estilo existente. El montaje
visual SHALL NOT reemplazar esas APIs por implementaciones paralelas.

#### Scenario: Interacción existente conserva los widgets reales
- **GIVEN** cualquiera de las siete páginas activada en SidePanel
- **WHEN** se usan su tabla, botón o señal funcional ya existente
- **THEN** se ejecuta el handler histórico y se actualiza el canvas/estado correspondiente
- **AND** no interviene una copia del contenido del dock.

### Requirement: Los controles de Validación caben en el panel embebido

El layout de presentación de Validación SHALL conservar las mismas instancias
históricas de botones y combo, señales y handlers, y SHALL mantener todos los
controles dentro de los límites visibles del panel lateral de hasta 340 px incluso
cuando se muestran navegación, exportación y corrección. El QSS compacto SHALL
aplicarse solo a widgets de docks embebidos.

#### Scenario: Controles de Validación no se recortan
- **GIVEN** Validación activa dentro del SidePanel
- **WHEN** se muestran los controles de navegación, reporte y corrección
- **THEN** cada control existente queda dentro de los límites del panel
- **AND** sus instancias, señales y handlers históricos se mantienen.

### Requirement: Barra de estado refinada sin datos inventados

La aplicación SHALL conservar una única instancia `QStatusBar` de 34 px, su
API `showMessage()`, `_update_status()`, el indicador IUPAC y el indicador de
carga. Su presentación SHALL usar jerarquía semántica y tokens light/dark
existentes. SHALL NOT mostrar fórmula, posición de cursor, estado de autosave
u otro valor que no tenga ya contrato barato y reutilizable; SHALL NOT añadir
polling ni cálculos químicos.

#### Scenario: Actualización de indicadores existentes y mensaje temporal
- **GIVEN** la ventana con SidePanel en modo light u oscuro
- **WHEN** cambia la herramienta/documento y se invoca `statusBar().showMessage()`
- **THEN** herramienta, IUPAC y carga siguen actualizándose por sus rutas actuales
- **AND** el mensaje temporal sigue accesible por `showMessage()`
- **AND** no se muestran datos de fórmula/cursor/autosave inventados.

### Requirement: Tema, splitter y tamaño mínimo utilizable

Los tabs, overflow y páginas embebidas SHALL seguir los tokens del tema
existente al cambiar entre light/dark. El SidePanel SHALL tener mínimo de
300 px; el lienzo conserva un mínimo de 400 px y el rail permanece intacto.
Ocultar el panel SHALL devolver el ancho disponible al lienzo y volver a
mostrarlo SHALL recuperar un ancho razonable. Con la ventana a 980×600, ningún
control SHALL quedar perdido ni el lienzo SHALL tener geometría negativa.

#### Scenario: Tema y ventana compacta
- **GIVEN** la ventana con SidePanel visible
- **WHEN** se alterna light→dark→light y se redimensiona a 980×600
- **THEN** el SidePanel permanece funcional y temáticamente coherente
- **AND** lienzo, rail y panel permanecen visibles sin geometría negativa ni superposición.
