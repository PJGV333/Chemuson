# Spec Delta — visibilidad de la fila superior de menús

## MODIFIED Requirements

### Requirement: La opción superior controla únicamente la fila de menús

La acción checkable de `Ver` SHALL mostrarse con el texto «Mostrar barra de herramientas» y SHALL controlar exclusivamente la visibilidad de la `QMenuBar` tradicional (Archivo, Editar, Ver, Estructura, Reacción y Ayuda). SHALL reutilizar la preferencia persistida `ui/show_menu_bar`; estará marcada por defecto para configuración nueva y la elección SHALL sobrevivir el reinicio. Desmarcar SHALL ocultar sólo la fila de menús, sin ocultar/mover la AppBar, perder acciones ni mostrar la QToolBar histórica. Al marcarla de nuevo SHALL reaparecer en su posición original.

La misma QAction dentro del mismo menú `Ver` SHALL permanecer accesible mediante la hamburguesa cuando la QMenuBar esté oculta, para que el usuario siempre pueda restaurarla. SHALL NOT existir una segunda opción de Ver con la misma finalidad que cambie la visibilidad de una QToolBar. La toolbar antigua de Nuevo/Abrir/Guardar SHALL permanecer oculta incluso si `ui/show_main_toolbar` está obsoleto y vale true.

#### Scenario: Nueva configuración
- **WHEN** una instalación nueva inicia con preferencias vacías
- **THEN** la QMenuBar está visible, la opción «Mostrar barra de herramientas» está marcada y la QToolBar antigua está oculta
- **AND** la AppBar compacta conserva sus acciones y pestañas.

#### Scenario: Ocultar la fila superior
- **WHEN** el usuario desmarca la opción
- **THEN** sólo desaparece la QMenuBar y se persiste `ui/show_menu_bar=false`
- **AND** la AppBar, sus acciones, pestañas y el panel lateral conservan su estado.

#### Scenario: Restaurar desde la hamburguesa
- **GIVEN** la QMenuBar está oculta
- **WHEN** el usuario abre la hamburguesa, accede a `Ver` y marca la misma QAction
- **THEN** la fila reaparece en su posición original y la opción queda marcada
- **AND** la QToolBar antigua sigue oculta.

#### Scenario: Persistencia
- **WHEN** la aplicación se reinicia tras ocultar o mostrar la fila
- **THEN** la QMenuBar y el estado marcado/no marcado reflejan la elección persistida
- **AND** un antiguo `ui/show_main_toolbar=true` nunca hace visible la toolbar de iconos.

#### Scenario: Acciones existentes y shell
- **WHEN** se alterna la visibilidad con tema claro/oscuro y tamaños habituales
- **THEN** Ctrl+S, Nuevo/Abrir/Guardar/Pegar, Undo/Redo, tabs, indicador dirty, X y panel lateral siguen disponibles y no duplicados
- **AND** ningún atajo ni callback se reconecta a una nueva ruta.

### Requirement: El popup Ver y sus submenús caben en pantalla

El popup de `Ver` SHALL colocar «Mostrar barra de herramientas» como primera acción visible, checkable y marcada cuando la fila de menús esté visible. El menú superior SHALL agrupar las opciones existentes en submenús concisos para que su geometría completa y las filas primera/última queden dentro de `QScreen.availableGeometry()`. Los submenús SHALL conservar las QActions originales y sus conexiones. La hamburguesa SHALL exponer esos mismos QMenus y la misma acción de recuperación.

La opción obsoleta «Mostrar copiar/zoom en barra de herramientas» SHALL estar ausente del menú y de la paleta de comandos; SHALL NOT eliminar ni recrear las acciones funcionales Copiar o Zoom, ni cambiar sus atajos/rutas. La QToolBar histórica SHALL permanecer oculta e incapaz de reaparecer por su toggle nativo o por preferencias antiguas.

#### Scenario: Menú compacto y completamente visible
- **WHEN** se abre `Ver` con interacción de mouse en tema claro u oscuro y con las ventanas configuradas a 1600×900, 1280×720 o 1024×768
- **THEN** la primera fila es la opción marcada «Mostrar barra de herramientas»
- **AND** el popup y las filas primera y última están completamente dentro de `QScreen.availableGeometry()`.

#### Scenario: Submenús accesibles
- **WHEN** se despliega cada submenú de Ver con el puntero
- **THEN** cada popup y sus opciones primera/última caben dentro del área disponible
- **AND** las QActions de visualización, zoom, numeración, lienzo, paletas y paneles son las existentes, sin duplicados.

#### Scenario: Restauración por mouse y reinicio
- **WHEN** el usuario oculta la fila haciendo clic en la primera acción y reinicia la aplicación
- **THEN** la fila continúa oculta, pero la hamburguesa permite abrir `Ver` y hacer clic en la misma acción para restaurarla
- **AND** el estado visible se persiste tras otro reinicio, con AppBar visible y QToolBar histórica oculta.

#### Scenario: Sin opción auxiliar obsoleta
- **WHEN** se abre `Ver` o se busca la opción en la paleta de comandos
- **THEN** «Mostrar copiar/zoom en barra de herramientas» no aparece
- **AND** Copiar permanece en Editar y Zoom en Ver con sus QActions y atajos existentes.
