# Diseño — depiction segura y ergonomía de dibujo

## Estado inicial

La rama es `feat/molecular-depiction-excellence`, HEAD `c38d954b9b7fb88a1aba04c06fc9376d5494ecc4`, sincronizada con su upstream. La Fase previa centra las rayas hash en `src/chemuson/gui/bond_junctions.py`; no modificó Clean2D ni límites del canvas. La galería Qt de uniones y sus pruebas son la base reutilizable. El inventario de baseline está en `baseline.md`.

## Primera reproducción Clean2D (sin parche)

En un ciclo cerrado de seis átomos construido con coordenadas manuales, variantes de enlaces simples/Kekulé/aromáticos produjeron candidatos seguros en quick/publication; la ruta real de ventana seleccionó `rdkit_isolated` (simple no aromático) o `simple_aromatic_template` (aromático), conservó todos los endpoints/órdenes/estilos y mejoró la razón de área del anillo de `0.753` a `0.907`, sin cruces. La API del motor tampoco muta el grafo por sí misma; la aplicación del controlador sí fue comprobada en el flujo GUI.

Se reprodujo una situación visualmente ambigua distinta: seis átomos en una forma casi cerrada pero con **sólo cinco enlaces** y terminales a `5.83` unidades. El grafo es una cadena, no un ciclo. Ctrl+K eligió `rdkit_isolated` y reubicó la cadena como zigzag; los seis átomos y los cinco endpoints permanecieron intactos. Esto explica cómo una figura que parecía un anillo puede «abrirse» al limpiar, pero no prueba que sea el mismo documento/caso del propietario. El test Qt de herramienta de enlace también confirma cierre cuando el usuario crea los seis enlaces y alcanza el átomo inicial.

Por tanto, todavía no hay evidencia de que Clean2D elimine o cambie conectividad de un ciclo de seis enlaces. No se añadirá un enlace automáticamente ni se cambiará la política de limpieza basándose en una figura de conectividad ausente. Queda pendiente obtener el archivo/caso exacto o fijar un repro manual reproducible de seis enlaces que falle.

## Resultado de la fase de etiquetas

La galería y la geometría mostraron que la separación se derivaba del margen fijo aplicado sobre bounding-ellipses y de sumar casi todo el ancho de cuña en su extremo ancho; el hash añadía además un intervalo vacío antes del primer mark. Se cambió el cálculo para proyectar el contorno real visible de etiqueta/carga sobre cada rayo y aplicar clearance según el tipo/ancho de trazo. La cuña sólida deja un margen nominal de `1.5` unidades; cada extremo hash se posiciona considerando el radio del trazo y mantiene al menos `1.0` unidad visible en los casos medidos. El conteo y espaciado uniforme interno de marks no cambian. Una matriz automatizada cubre F/Cl/Br/N⁺/CH₂OH, cuatro ángulos, punta/base y enlaces simple/doble/triple/bold/wedge/hash.

La brecha se cerró sin intersecciones en la matriz (0.5–5.0 unidades); la lámina before/after mantiene juntas las 20 uniones originales y añade 20 orientaciones de etiqueta, en claro y oscuro. Esto no concluye la revisión visual del propietario. Tampoco se presupone que toda forma de bifurcación sea un defecto: sólo se ajustó el trimming que explicaba la métrica observada.

## Secuencia y límites técnicos

1. **Clean2D primero**: construir repros con ciclo de seis enlaces dibujado manualmente (incluyendo enlaces simples no aromáticos y aromaticidad explícita) y plantilla equivalente. Registrar clasificación, política, selección, fuentes/rechazos, estado final, grafo, coordenadas, longitudes, cruces y degeneración del anillo en quick/publication/propose según proceda. Si un candidato aceptado visualmente abre/colapsa el ciclo, rastrear su fuente y compuertas; no relajar invariantes para forzar un resultado.
2. **Etiquetas y bifurcaciones**: extender la galería existente con F/Cl/Br y etiquetas anchas, cuña sólida/hash, enlaces simple/múltiple/bold y dos ramas salientes. Comparar el borde visible del glyph/label con el trazo y el contorno de cuña en la misma escala. Cambiar sólo la operación de trimming/unión que explique la métrica observada; conservar el sentido del enlace estéreo.
3. **Borde del canvas**: seguir el flujo de eventos de herramientas y el cálculo de `sceneRect`. Probar commits cerca del viewport con ítems parcialmente fuera de la zona visible, y verificar que el viewport puede navegar al resultado. No limitar posiciones ni reordenar despacho de eventos.
4. **Anillos/Haworth**: probar la ruta actual de plantillas, snapping y rotación/colocación con Undo. Un ajuste se limita a geometría/local snapping explícitamente reproducible; preservar topología y semántica de sustituyentes.

## Evidencia y arquitectura

Reusar `BondItem`, `wedge_geometry.py`, `canvas_structure.py`, `canvas_view.py`, `canvas_tools_rings_chains.py` y los módulos Clean2D existentes; no extraer paquetes ni sumar dependencias sin OpenSpec y catálogo. Mantener render Qt de referencia. Conservar capturas exploratorias y comparaciones intermedias bajo `/tmp`; el repositorio tendrá el generador reproducible y una única lámina comparativa.

## Riesgos y gates

- Una coordenada candidata puede mantener los enlaces en el grafo pero colapsar el polígono del anillo; validar ambos aspectos de forma independiente.
- Trimming excesivo puede invadir/ocultar glyphs, cortar hash strokes o alterar la lectura estereoquímica.
- La geometría de cuña en una bifurcación puede ocultar enlaces vecinos; no cambiarla por simetría estética sin pruebas.
- Expandir `sceneRect` puede mover el centro o saltar el scroll; preservar la posición visible y probar navegación/Undo.
- No declarar una fase finalizada por una captura aislada; la comparación debe incluir tema claro/oscuro cuando aplique.
