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
