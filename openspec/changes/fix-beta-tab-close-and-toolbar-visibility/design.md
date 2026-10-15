# Diseño

1. **Icono de cierre ausente (diagnóstico Qt).** El primer documento obtiene icono
   X cuando el tema de la ventana se aplica al final del ensamblaje, pero
   `_add_tab()` creaba después los botones laterales sin `setIcon()`. Cambiar de
   tema los refresca, lo que explica que aparecieran sólo tras hover/cambio de
   tema. `_TabSideButton` inicializa su icono con el `IconProvider` compartido y
   el tinte vigente; el mismo `refresh_icons()` sigue actualizando el tinte ante
   cambios de tema.
2. **A3 — causa raíz del desplazamiento de la X.** El bloque lateral
   (`_TabSideButton`) medía 28×24 px mientras su layout exigía 40 px
   (márgenes 2+2 + punto 7 + separación 5 + botón 24). La QSS
   `QToolButton[tabClose="true"]` fijaba `min-width/min-height: 24px`, que en el
   modelo de caja QSS se suma al borde (1 px) y al `setFixedSize(24, 24)`,
   inflando el botón a 26×26. Con la pestaña limpia el botón ocupaba todo el hueco
   (`x=2`); al aparecer el punto dirty, Qt le restaba 12 px y lo empujaba a
   `x=14`, desbordando el rect de la pestaña (14+26 > 28). Corrección: el bloque
   lateral pasa a tamaño fijo `_SIDE_WIDTH × _CLOSE_BTN_SIZE` (40×24), el punto
   dirty **ocupa siempre su celda** y sólo cambia su pintura mediante la propiedad
   dinámica `dirty` (`#dirtyDot[dirty="false"]` transparente, con
   `unpolish/polish`), y la QSS deja de imponer `min-width/min-height` al botón de
   cierre. Con eso `tabRect`, tamaño del bloque, tamaño de la X y su posición
   dentro del bloque son idénticos en limpio y sucio, y la holgura respecto al
   borde derecho de la pestaña es constante.
3. **A1 — barra de menús tradicional.** Se conserva el único `QMenuBar` de
   `QMainWindow` con los seis `QMenu` originales; el agregador hamburguesa sigue
   apuntando a esos mismos objetos (`_menu_aggregator.addMenu(top_menu)`), así que
   no hay una segunda copia de menús, acciones ni atajos. El default pasa a
   **visible** (`UiPreferences.show_menu_bar = True`) y la QSS compacta la fila:
   `QMenuBar { padding: 1px 6px; min-height: 22px }` y
   `QMenuBar::item { padding: 3px 10px; margin: 0px 1px }` → 30 px de alto
   (baseline: 43 px), por debajo de los 54 px de la app bar y sin solapamiento.
   La hamburguesa permanece visible en ambos modos (es también el acceso por Alt).
4. **Migración de defaults.** Un INI escrito antes de este cambio puede contener
   `ui/show_menu_bar=false` generado por el default anterior, indistinguible de
   una elección explícita. Se añade el marcador `ui/shell_defaults_version`
   (`SHELL_DEFAULTS_VERSION = 1`): mientras no esté escrito, `load_ui_preferences`
   devuelve los defaults nuevos de la beta (ambas barras visibles); el primer
   guardado de preferencias (toggle, tema) sella el marcador y a partir de ahí se
   respeta lo persistido. La migración es idempotente y no toca el tema.
5. **A2 — barra de tres botones.** `build_main_toolbar` añade únicamente
   `action_new`, `action_open`, `action_save` (mismas QActions de Archivo, mismos
   atajos y mismos controladores). Se retiran de esa barra undo/redo, rotaciones,
   volteos, Clean2D y SMILES: siguen en Editar/Estructura, sus atajos, el rail y la
   paleta de comandos; Undo/Redo quedan representados en la app bar
   (`setDefaultAction`). La barra se compacta con `objectName("mainToolBar")` +
   QSS dedicada (`padding: 3px`, `spacing: 4px`, botones `padding: 4px`, icono 18
   px) → 37 px de alto y botones uniformes de 31×31 (baseline: 65 px). La QSS
   general de `QToolBar`/`QToolButton` (paletas verticales) no se toca.
   La opción heredada `Mostrar copiar/pegar/zoom en barra superior` se conserva
   como opt-in explícito: es la única que puede ampliar esa barra y no cambia su
   visibilidad global.
6. **Preferencias y despacho.** `Ver` mantiene al inicio las dos acciones
   checkable distintas (`Mostrar barra de menús`, `Mostrar barra de herramientas`
   = `main_toolbar.toggleViewAction()`), cada una persistida por separado
   (`ui/show_menu_bar`, `ui/show_main_toolbar`) e independiente de
   `ui/side_panel/visible`. No se crean QActions ni señales paralelas.
7. **Riesgos de regresión y pruebas.** La QSS de `QMenuBar` afecta a todos los
   temas: se comprueba altura ≤ 36 px y geometría de los seis ítems en claro/oscuro
   y a 1440×900, 980×600 y 900×560. Retirar acciones de `main_toolbar` podría
   romper contratos de tests previos: se actualizan explícitamente los que
   codificaban la barra de nueve botones y el modo compacto por defecto. Los tests
   de defaults usan el fixture `isolated_shell_settings` (INI temporal) porque
   varios tests construían `ChemusonWindow` con el INI real del equipo. Suite
   amplia en shards con timeout externo ≤ 600 s por invocación.
