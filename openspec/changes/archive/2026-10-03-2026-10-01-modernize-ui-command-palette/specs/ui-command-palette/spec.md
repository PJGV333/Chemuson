# UI Command Palette Specification

## Purpose

Define la paleta de comandos (Ctrl+P) de la Fase 6 como una superficie de
presentación, filtro y ejecución de las `QAction` existentes de la aplicación,
sin duplicar lógica química ni handlers. La paleta se abre con `Ctrl+P` (única
`QAction` global que lo posee, `action_command_palette`) y con el clic en la
píldora de búsqueda del AppBar. `Ctrl+K` es el atajo histórico de
`action_clean_2d_full` (limpia 2D, 1 paso) y NO abre la paleta.

## ADDED Requirements

### Requirement: La paleta presenta, filtra y ejecuta QActions existentes

La aplicación SHALL exponer `CommandPalette` en
`src/chemuson/gui/command_palette.py` que presente, filtre y ejecute `QAction`
ya existentes en la ventana. Al ejecutar una entrada, la paleta SHALL disparar
`triggered` de la `QAction` objetivo (o `activate` equivalente) y SHALL NO crear
una segunda `QAction` para representar una `QAction` existente ni duplicar
handlers. La paleta SHALL conservar en cada fila la identidad, el estado
`enabled`, el estado `checkable`/`checked`, el shortcut mostrado y el icono de
la `QAction` correspondiente.

#### Scenario: Ejecutar una entrada dispara su QAction histórica
- **GIVEN** una `ChemusonWindow` con una `CommandPalette` montada
- **WHEN** el usuario selecciona la entrada "Limpiar 2D (1 paso)" y pulsa Enter
- **THEN** se dispara `triggered` de `action_clean_2d_full` exactamente una vez
- **AND** la paleta se cierra
- **AND** no existe ninguna `QAction` duplicada para ese comando.

### Requirement: El registro deduplica por identidad de QAction

`CommandRegistry` SHALL indexar entradas por identidad de `QAction`
(`id(action)`), no por texto. Registrar la misma `QAction` más de una vez SHALL
no crear entradas duplicadas. El registro SHALL contener al menos 60 comandos
únicos y realmente ejecutables/buscables.

#### Scenario: No duplica la misma QAction
- **GIVEN** el registro construido por la ventana
- **WHEN** se enumeran las `QAction` representadas
- **THEN** ninguna `QAction` aparece en dos entradas distintas
- **AND** el número de entradas únicas es ≥ 60.

#### Scenario: Reabrir la paleta no crea registros duplicados
- **GIVEN** una `ChemusonWindow`
- **WHEN** se abre y cierra la paleta varias veces
- **THEN** el número de entradas del registro permanece constante.

### Requirement: Las siete páginas del SidePanel son buscables

La paleta SHALL exponer las siete páginas del SidePanel (Inspector, Validación,
Propiedades químicas, Plantillas, Apariencia, Espectroscopía y CompChem)
reutilizando las `QAction` creadas en el menú Ver
(`window.side_panel_actions`), sin llamar a `side_panel.show_page()` desde una
copia paralela.

#### Scenario: Cada página del panel es alcanzable desde la paleta
- **GIVEN** la ventana con `side_panel_actions` pobladas
- **WHEN** se busca "Inspector" (y análogamente cada una de las 7 etiquetas)
- **THEN** aparece la entrada correspondiente
- **AND** ejecutarla activa la misma `QAction` del menú Ver (mismo objeto).

### Requirement: Las exportaciones mínimas son buscables

La paleta SHALL exponer como entradas las exportaciones existentes representadas
por `QAction` reales: PNG, SVG, PDF, CML y SMILES, además de las demás
exportaciones ya presentes como `QAction`.

#### Scenario: Buscar "export" lista las exportaciones
- **GIVEN** la ventana montada
- **WHEN** se busca "export"
- **THEN** aparecen las entradas PNG, SVG, PDF, CML y SMILES (entre otras)
- **AND** cada una apunta a su `QAction` de exportación correspondiente.

### Requirement: Ctrl+P abre la paleta; Ctrl+K ejecuta Clean2D quick

La ventana SHALL exponer una única `QAction` de apertura
(`action_command_palette`) con `QKeySequence("Ctrl+P")`, contexto
`WindowShortcut`, registrada mediante `window.addAction(...)`, cuyo `triggered`
abre la `CommandPalette`. `Ctrl+P` SHALL poseerla **exclusivamente** esa
`QAction` global (sin conflictos). `action_clean_2d_full` SHALL poseer el
atajo histórico `QKeySequence("Ctrl+K")` (contexto `WindowShortcut`,
`window.addAction(...)`), cuyo `triggered` ejecuta *Limpiar 2D (1 paso)*;
`Ctrl+K` SHALL NO abrir la paleta. `Ctrl+Shift+K` y `Ctrl+Alt+K` SHALL
permanecer intactos.

#### Scenario: Ctrl+P abre la paleta
- **GIVEN** una `ChemusonWindow` visible
- **WHEN** se pulsa Ctrl+P
- **THEN** la `CommandPalette` se muestra (foco en el input)
- **AND** `action_clean_2d_full.triggered` NO se dispara.

#### Scenario: Ctrl+K ejecuta Clean2D quick y no abre la paleta
- **GIVEN** una `ChemusonWindow` visible
- **WHEN** se pulsa Ctrl+K
- **THEN** `action_clean_2d_full.triggered` se dispara exactamente una vez
- **AND** la `CommandPalette` NO se muestra.

#### Scenario: Clean2D quick sigue accesible por menú y por paleta
- **GIVEN** la ventana
- **WHEN** se busca "Limpiar 2D (1 paso)" en la paleta y se ejecuta
- **THEN** `action_clean_2d_full.triggered` se dispara
- **AND** la acción sigue presente en *Estructura → Limpiar 2D (1 paso)*.

#### Scenario: Solo una QAction global posee Ctrl+P
- **GIVEN** la ventana montada
- **WHEN** se enumeran las `QAction` de la ventana con `QKeySequence("Ctrl+P")`
- **THEN** aparece exactamente una: `action_command_palette`.

#### Scenario: Ctrl+Shift+K y Ctrl+Alt+K conservan su comportamiento
- **GIVEN** la ventana
- **WHEN** se pulsa Ctrl+Shift+K y luego Ctrl+Alt+K
- **THEN** `action_clean_2d_publication` y `action_clean_2d_propose` se disparan
  respectivamente.

### Requirement: La píldora del AppBar abre la misma ruta de Ctrl+P

La `SearchPill` del AppBar SHALL convertirse en entrada real a la paleta: su
clic y Ctrl+P SHALL abrir exactamente la misma instancia/ruta (una misma
`QAction` de apertura). El badge visual `Ctrl P` SHALL mostrarse y el tooltip
SHALL dejar de indicar "(próximamente)". La `SearchPill` SHALL NO convertirse
en un editor permanente; el campo editable pertenece a `CommandPalette`.

#### Scenario: Clic en la píldora abre la misma ruta que Ctrl+P
- **GIVEN** la ventana
- **WHEN** se hace clic en `app_bar.search_pill`
- **THEN** se dispara la misma `QAction` de apertura (`action_command_palette`)
  que Ctrl+P
- **AND** la paleta se muestra.

#### Scenario: La píldora muestra el hint Ctrl P y no es un editor
- **GIVEN** el AppBar montado
- **WHEN** se inspecciona la `SearchPill`
- **THEN** el badge `Ctrl P` es visible
- **AND** el tooltip no contiene "(próximamente)"
- **AND** la `SearchPill` no contiene un `QLineEdit`.

### Requirement: Filtro substring con prioridad por prefijo

La paleta SHALL filtrar case-insensitive por substring sobre el título y los
keywords (y sección). Una coincidencia por **prefijo** SHALL rankear antes que
una coincidencia por substring. El ranking SHALL ser simple, determinista y
testeable. No se usa fuzzy-search externo.

#### Scenario: El prefijo rankea antes que el substring
- **GIVEN** entradas "Exportar SMILES" y "Importar SMILES"
- **WHEN** se busca "export"
- **THEN** "Exportar SMILES" (prefijo) aparece antes que cualquier entrada que
  solo contenga "export" como substring en otra posición.

#### Scenario: El matching por keywords funciona
- **GIVEN** una entrada cuyo título no contiene el query pero cuyo keyword sí
- **WHEN** se busca el keyword
- **THEN** la entrada aparece en los resultados.

### Requirement: Navegación y ejecución por teclado

Con la paleta abierta, `↑`/`↓` SHALL mover la selección; `Enter` SHALL
ejecutar la fila seleccionada exactamente una vez y cerrar la paleta; `Esc`
SHALL cerrar sin ejecutar; el clic en una fila SHALL ejecutar esa fila. Las
`QAction` disabled SHALL no ejecutarse y SHALL reflejarse visualmente.

#### Scenario: ↑/↓ cambia la selección
- **GIVEN** la paleta abierta con ≥ 2 resultados
- **WHEN** se pulsa ↓ y luego ↑
- **THEN** la fila seleccionada cambia correspondientemente.

#### Scenario: Enter dispara exactamente una vez la QAction elegida
- **GIVEN** la paleta abierta con una fila seleccionada habilitada
- **WHEN** se pulsa Enter
- **THEN** la `QAction` seleccionada se dispara exactamente una vez
- **AND** la paleta se cierra.

#### Scenario: Esc cierra sin disparar
- **GIVEN** la paleta abierta con una fila seleccionada
- **WHEN** se pulsa Esc
- **THEN** la paleta se cierra y ninguna `QAction` se dispara.

#### Scenario: Una QAction disabled no se ejecuta
- **GIVEN** una fila que representa una `QAction` disabled
- **WHEN** se selecciona y se pulsa Enter (o clic)
- **THEN** la `QAction` no se dispara.

#### Scenario: Una QAction checkable conserva su semántica
- **GIVEN** una fila checkable (p. ej. "Mostrar carbonos")
- **WHEN** se ejecuta desde la paleta
- **THEN** el estado `checked` de la `QAction` se alterna por delegación a su
  `triggered`.

### Requirement: Presentación visual dentro de la ventana

La paleta SHALL abrirse centrada sobre la ventana, con ancho objetivo de 560 px
limitado a la ventana, superficie/borde por tokens, input superior y lista de
resultados con icono, nombre, sección y shortcut cuando existan. La selección
SHALL usarse con `accent`/`accentSoft`. El estilo SHALL ser coherente en light y
dark usando solo tokens del sistema (sin colores hardcodeados fuera de tokens).
SHALL funcionar correctamente al menos en 1440×900 y 980×600, y en ventanas
estrechas no SHALL salirse de pantalla.

#### Scenario: 980×600 mantiene la paleta dentro de la ventana
- **GIVEN** la ventana a 980×600
- **WHEN** se abre la paleta
- **THEN** la tarjeta queda íntegramente dentro del rectángulo de la ventana.

#### Scenario: Light → dark → light no rompe la paleta
- **GIVEN** la ventana
- **WHEN** se cicla light→dark→light y se abre la paleta
- **THEN** la paleta se re-estiliza correctamente y sigue funcional.

#### Scenario: Los atajos de herramienta no interfieren al escribir
- **GIVEN** la paleta abierta con foco en su `QLineEdit`
- **WHEN** se teclean letras (p. ej. "b", "r", "c")
- **THEN** se escriben en el input y no se activa ninguna herramienta del rail.

### Requirement: Fronteras arquitectónicas de la paleta

`command_palette.py` SHALL importar únicamente `QAction`/metadatos, `theme`
(tokens/QSS/IconProvider) y el shell. SHALL NOT importar `chemuson.clean2d`,
`chemuson.chemname`, `chemuson.chemio.persistence`, internals de `gui.canvas` ni
controllers químicos concretos.

#### Scenario: Contrato de imports de la paleta
- **GIVEN** `src/chemuson/gui/command_palette.py`
- **WHEN** se analiza su conjunto de imports
- **THEN** no contiene import alguno de `chemuson.clean2d`,
  `chemuson.chemname`, `chemuson.chemio.persistence`, `chemuson.gui.canvas` ni de
  un controller químico concreto.
