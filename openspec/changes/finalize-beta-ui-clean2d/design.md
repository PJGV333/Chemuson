# Diseño — consolidación UI + Clean2D

## 1. Integración trazable

La rama conjunta parte de `9bc8192db3835262191ddb5247a8ff87f4d9169f` e incorpora con merges no fast-forward las ramas remotas originales, sin reescribirlas. La unión del manifiesto CI debe coincidir exactamente con la colecta del árbol combinado; los node IDs antiguos de tests UI sustituidos por la fase A no se conservan como entradas obsoletas, y todas las pruebas aún colectables sí.

## 2. Cabecera

La fila tradicional mantiene `Archivo | Editar | Ver | Estructura | Reacción | Ayuda`, visible por defecto. La fila de documentos/acciones conserva las pestañas, controles de accesibilidad y acciones actuales. Se retiran de la AppBar el QLabel de marca, el pixmap de matraz y el QLabel de versión, sin trasladarlos a otro lugar. El título principal, metadatos y Acerca de siguen con nombre/versión oficiales.

Las acciones Nuevo, Abrir, Guardar y Pegar son las mismas QAction ya conectadas a sus menús/controladores. La acción Pegar aparece una sola vez en Editar y su botón compacto de la AppBar usa esa misma QAction; no se crea handler ni shortcut alternativo. El `+` de pestañas se oculta en la ventana real para que Nuevo tenga un único botón dedicado en el grupo de acciones; el widget aislado conserva su API de prueba. Undo/Redo mantienen solamente sus botones existentes a la derecha. La QToolBar superior migra al default oculto (versión de defaults 2) y permanece opcional desde Ver; sus acciones auxiliares conservan Copiar/Zoom, pero no repiten Pegar. La fila de menú y la AppBar son las dos franjas por defecto.

La prueba visual cuantifica geometrías, botones y listas de QActions en temas claro/oscuro y tamaños compactos. La fila del menú y la AppBar juntas son el header de dos filas; no se altera el panel ni la jerarquía general de GUI.

## 3. Diagnóstico Clean2D reproducible

El primer smoke GUI sobre `C[C@H](F)CO` importado como grafo real mostró una cuña `HASHED` y 5 átomos; `classify_clean2d_complexity()` informó `stereo_center_count=2`, `multiple_stereo_centers`, `preserve_only=true`. `run_clean2d_engine(quick)` retornó sólo `current`, sin candidatos geométricos ni movimiento; el filtro geométrico, backend y comparación de firma no fueron la causa porque el retorno temprano sucede antes de generarlos. `Ctrl+K` llega a la acción de limpieza completa (modo engine `quick`) y muestra el mismo mensaje. La causa concreta está en `core.layers`: cada extremo del mismo bond `UP/DOWN` generaba un bloque `STEREO_CENTER`; la política existente interpretaba dos extremos como dos centros.

Se corregirá la construcción de bloques para agregar cada centro una vez por ID. Los `stereo_cip` explícitos cuentan por átomo y los enlaces estereo dirigidos cuentan en el extremo de origen (`a1_id`), acumulando varios enlaces que comparten centro en el mismo bloque. Un centro no asignado químicamente no recibe R/S por el mero hecho de contener una cuña; la firma de layout continúa protegiendo la representación visual. La regla `>=2 centros distintos → preserve_only` permanece igual.

Se registra por separado: una molécula real desconectada de otras y selección de sólo una; el centro con wedge+hash; selección total; selección parcial incluyendo centro y enlaces frontera; quick/Ctrl+K, publicación/Ctrl+Shift+K y propuesta/Ctrl+Alt+K. El diagnóstico conserva todos los candidatos con `source`, `rejected`, `rejection_reason` y motivo estable; no se aplicará un candidato si no pasa los gates existentes.

## 4. Verificación estereoquímica independiente

Los contratos visuales (orden cíclico y estilos de enlace) no sustituyen la configuración química. La prueba de la GUI detectó que un fallback de propuesta en Publicación puede mantener el orden cíclico, pero cambiar la asignación R/S al serializar la geometría alrededor de un H implícito. Por ello, antes de aplicar cualquier candidato GUI sobre cuñas/hash, el controlador pide al worker RDKit aislado las asignaciones independientes anterior y posterior, quitando CIP cacheado de las copias de verificación. Si un R/S asignado cambia o deja de ser asignable, el candidato se rechaza. Si el baseline no permite asignar R/S (como una cuña sobre grupos equivalentes), se conservan los gates geométricos sin inventar quiralidad.

Las pruebas exportan/reimportan los enantiómeros `C[C@H](F)CO` y `C[C@@H](F)CO` en Quick, Publicación y propuesta; comparan CIP independiente, SMILES isomérico, firma visual, enlaces y coordenadas fuera de selección. No se modifica el worker ni la conversión ChemIO de producción.

## 5. Restricciones y reversión

Solo cambian AppBar/shell/settings necesarios, `core.layers` para la deduplicación estructural demostrada, el gate R/S aislado del controlador GUI y las pruebas/OpenSpec/manifiesto. No se cambia `complex_policy` ni el motor de `local_graph_cleaner`; tampoco se relajan filtros geométricos. Si los candidatos no pasan el gate R/S o de orden dibujado, se conserva la estructura y se informa la razón.
