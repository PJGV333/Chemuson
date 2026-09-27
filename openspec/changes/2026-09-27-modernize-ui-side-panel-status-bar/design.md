# Design: Panel lateral moderno y barra de estado (Fase 5)

## Contexto y auditoría previa

Ver `baseline.md` para el estado Git, entorno y salidas de validación previas.
El spike aprobado define `SideTabRow + QStackedWidget`, ancho `sideW=324` y
una fila de tabs desplazable; los cinco labels principales de esta fase son
Inspector, Validación, Propiedades, Plantillas y Apariencia. El contenido de
demo del spike no se migra.

`docks.py` define las siete clases funcionales históricas. `shell/assembly.py`
las instancia y conecta sus señales; `main_window.py` actualiza sus datos y
atiende acciones. `main_window_ui_builder.py` actualmente incorpora sus
`toggleViewAction()` al menú Ver. Las rutinas de validación también llaman
`validation_dock.show()` directamente y deberán usar la ruta pública del panel.

La barra de estado ya es un `QStatusBar` de 34 px: `_update_status()` llama a
`showMessage()`, mientras `_iupac_name_label` y `_total_charge_label` son
widgets permanentes actualizados por las rutinas existentes. No se añadirá
ninguna fuente de datos nueva.

## Decisiones

### D1 — Una sola superficie lateral; cinco páginas principales y overflow

`SidePanel` (`gui/side_panel.py`) se monta como hermano del `QTabWidget` de
documentos en `_central_body_layout`, después del ToolRail. El ancho objetivo
es `METRICS["sideW"] = 324` px, dentro del rango 300–340 px. La fila superior
reutiliza el patrón del spike: tabs checkables en una tira con scroll
horizontal sin scrollbar visible, más un botón `…` fijo. El menú `…` contiene
Espectroscopía y CompChem. Seleccionar una página de overflow mantiene `…`
visualmente activo y selecciona el widget correspondiente en el mismo stack.

El stack inicial es Inspector. `SidePanel.show_page(key)` es la única ruta
pública para navegación: valida la clave, sincroniza tab/stack, hace visible el
panel y persiste página y visibilidad. `set_panel_visible(bool)` centraliza
ocultar/mostrar y persistir. Las acciones del menú Ver y las rutas internas que
antes mostraban Validación delegan en estos métodos.

### D2 — Reutilizar cada QDockWidget como página, no extraer/copiar su contenido

Las siete instancias existentes se crean una sola vez y se pasan al
`QStackedWidget` como páginas. `QDockWidget` es un `QWidget`; el experimento
Qt offscreen previo al cambio confirmó que el mismo objeto puede añadirse a un
`QStackedWidget` conservando tanto `dock.widget()` como el widget de contenido.
El objeto no se registra en `QMainWindow.addDockWidget()`, recibe
`NoDockWidgetFeatures` y un title bar vacío de altura cero, para eliminar el
chrome clásico. No se crea un segundo inspector, validador, tabla ni handler.

Los atributos de compatibilidad en la ventana conservan exactamente las
instancias (`window.inspector_dock`, etc.). Sus métodos, señales, tablas,
botones y conexiones actuales permanecen sin cambios. El panel sólo coordina
página, contención, visibilidad y tema. Los tests de contrato comprueban
identidad de objeto, identidad del contenido y ausencia de páginas duplicadas.

### D3 — Páginas y acciones del menú Ver

Claves estables: `inspector`, `validation`, `properties`, `templates`,
`appearance`, `spectroscopy`, `compchem`. Se corresponden 1:1 con las siete
instancias existentes. Ver contiene acciones normales (no toggle de dock) con
los labels solicitados: Inspector, Validación, Propiedades químicas,
Plantillas, Apariencia, Espectroscopía y CompChem. Todas llaman a
`side_panel.show_page(key)`; ninguna llama `show()`, `raise_()` o
`toggleViewAction()` sobre un dock clásico.

Se conserva la QAction de la toolbar de símbolos en Ver. La acción de
preferencias del AppBar permanece conectada a `PreferencesDialog` porque la
auditoría del código confirma que ése es su contrato actual; el acceso del
AppBar a Apariencia es la ruta hamburguesa → Ver → Apariencia. El menú también
ofrece una acción de visibilidad del panel para que el usuario pueda
ocultarlo/restaurarlo sin convertir las siete acciones de página en toggles.

### D4 — Preferencias del panel

`SidePanelPreferences(active_tab="inspector", visible=True)` y las funciones
`load_side_panel_preferences()` / `save_side_panel_preferences()` viven en
`chemuson.platform.settings`, junto a `ui/theme`. Las claves son
`ui/side_panel/active_tab` y `ui/side_panel/visible`. Tabs desconocidos y
valores booleanos inválidos se normalizan al default seguro: Inspector y visible.
No se implementa orden movible.

### D5 — Barra de estado: pulido sin sustituir datos ni API

Se conserva la instancia de `QStatusBar`, altura de 34 px, `showMessage()`,
`_update_status()`, `_iupac_name_label` y `_total_charge_label`. Se asignan
object names semánticos y reglas QSS de tokens para distinguir herramienta a
la izquierda, IUPAC con espacio flexible y carga a la derecha, manteniendo los
widgets existentes y su orden. No se cambia el cálculo ni el ciclo de
actualización de esos indicadores.

Fórmula molecular, cursor y autosave quedan diferidos: la implementación
actual no ofrece para ellos un contrato de UI directo, barato y transversal
que cumpla estas restricciones sin ampliar la fase.

### D6 — Arquitectura y límites

`side_panel.py` sólo importa widgets Qt, `docks.py`, tema/tokens e
infraestructura de preferencias; no importa canvas internals ni Clean2D, no
calcula química y no crea docks funcionales. Se añade el path y la nota de
responsabilidad a M08 (`gui`) en `architecture/modules.yml`; no se crea un
paquete nuevo ni una dependencia externa.

No se modifica AppBar, la lógica ni los contratos de `docks.py`, la jerarquía
de controllers o los contratos químicos. La fila de controles de Validación
solo se adapta en presentación para el ancho del panel: se conservan los mismos
botones, combo, señales y handlers, refluídos en filas compactas; los cambios
QSS se limitan a widgets de docks embebidos.

### D7 — Reflujo de controles de Validación con widgets históricos

La captura de 324 px confirmó que el QSS global (mínimo de 80 px más padding
horizontal de 22 px por botón) hacía que la fila histórica de Validación
recortara `Siguiente`; al mostrar reportes y correcciones el desbordamiento era
mayor. Se mantiene cada instancia y conexión existente, pero el layout de
presentación agrupa navegación, reportes y corrección en filas separadas. QSS
específico de docks embebidos reduce padding y mínimo horizontal de botones y
combo. Esto no cambia handlers, datos ni señales y evita controles fuera del
límite visible del SidePanel.

## Verificación visual

Las capturas offscreen cubren light/dark 1440×900 en Inspector, light
Validación, dark Propiedades, overflow abierto, Apariencia y 980×600. Se
adjunta un script reproducible bajo la carpeta de capturas. Sirven para
inspección de regresiones visuales; la aceptación final exige prueba manual
por el usuario en KDE/Wayland real y no se declara por el agente.
