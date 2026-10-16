# Diseño — depiction segura y ergonomía de dibujo

## Estado inicial

La rama es `feat/molecular-depiction-excellence`, HEAD `c38d954b9b7fb88a1aba04c06fc9376d5494ecc4`, sincronizada con su upstream. La Fase previa centra las rayas hash en `src/chemuson/gui/bond_junctions.py`; no modificó Clean2D ni límites del canvas. La galería Qt de uniones y sus pruebas son la base reutilizable. El inventario de baseline está en `baseline.md`.

La observación de que Clean2D «abre» un anillo aún no identifica si el defecto está en el grafo, un candidato geométrico, una validación incompleta, una selección parcial o sólo en la visualización. No se presume conectividad alterada hasta comparar snapshots del grafo antes/después y coordenadas elegidas. Tampoco se presupone que toda separación cuña-etiqueta o forma de bifurcación sea un defecto: primero se mide y compara con etiquetas y ángulos variados.

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
