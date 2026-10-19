# Propuesta — geometría 2D inicial de calidad editorial

## Por qué

El dibujo esqueletal inicial de cadenas puede heredar el ángulo tetraédrico de 109.5° y no la convención 2D de zigzag de 120°. Además, se reprodujeron dos fallos en la asistencia de dibujo: el snapping absoluto puede alejar un sustituyente de la geometría local al rotar un anillo, y un clic corto en un átomo de anillo puede resolver un enlace vecino en vez de crear el nuevo enlace. Esto impide obtener una estructura ordenada desde el primer trazado.

La auditoría separó estos defectos de la geometría del benceno y del renderer: el hexágono producido por la herramienta de anillo ya tiene lados de 40 unidades y ángulos internos de 120° dentro de error de punto flotante; las capturas del renderer actual no muestran una irregularidad reproducible del anillo aislado. No se justifica cambiar los offsets de dobles/aromáticos ni otros avances visuales.

## Alcance

- Definir una política determinista de ángulos esqueletales 2D, separada de la referencia tetraédrica 109.5° y de cualquier geometría 3D.
- Usar 120° interiores para el zigzag de cadenas acíclicas y las sugerencias automáticas sp3; mantener el ángulo trigonal sp2 de 120°, sp lineal de 180° y el polígono regular propio de cada tamaño de anillo.
- Hacer que una sugerencia por clic corto respete el contexto angular local y coincida entre preview y commit; el snapping global se mantiene para la orientación explícita mediante arrastre y para los controles existentes.
- Corregir el hit-test de liberación que bloquea un clic corto iniciado sobre un átomo con enlaces vecinos.
- Añadir pruebas numéricas, gestuales, de undo/redo y de render a escala 1×/1.25×, y conservar un registro antes/después reproducible.

## Fuera de alcance

No cambiar Clean2D ni sus presets, QSS, fuentes, espesores, junctions, cuñas, offsets/representación de enlaces dobles o aromáticos, proyecciones, química/modelos, formato `.cmsn`, selección/importación antigua, versión, empaquetado, OpenGL 3D, Live MM ni la política de publicación. No añadir optimizador global, dependencias o paquetes nuevos. No alterar la intención de arrastre ni modificar coordenadas ya existentes por insertar un átomo.

## Módulos probables

`src/chemuson/gui/geom.py`, `canvas/canvas_bond_state.py`, `canvas/canvas_bond_geometry.py`, `canvas/canvas_bond_drag.py`, `canvas/canvas_tools_rings_chains.py` y pruebas Qt/numéricas. Los comandos existentes de enlace/anillo/cadena y `BondItem`/`AtomItem` se auditan para preservar el contrato; no se prevé modificar renderizado ni arquitectura de paquetes.
