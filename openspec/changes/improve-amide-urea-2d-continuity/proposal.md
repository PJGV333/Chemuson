# Propuesta — continuidad 2D de amidas y ureas

## Por qué

La geometría inicial 2D aprobada conserva longitudes de enlace y ángulos razonables al construir amidas/ureas, pero los enlaces a `NH` quedan visualmente truncados. La primera reproducción con eventos reales del editor dio 40 unidades entre centros y ángulos trigonal-planos, mientras dos enlaces NH–carbonilo mostraron sólo ~24.3 unidades de trazo visible. La causa sospechada es medir la separación texto-enlace por la proyección de todos los puntos de un glifo, aunque estén lejos de la trayectoria real del enlace.

La asistencia local también debe distinguir N de amida/urea por conectividad y orden de enlace, sin basarse en los textos visibles `NH`/`NH2` ni imponer la misma geometría a aminas, N aromáticos o centros de otra coordinación.

## Alcance

- Reproducir con QTest el fragmento `Ar–NH–C(=O)–NH–CH3`, variantes de urea/amida y orientaciones de etiquetas; separar coordenadas, dirección, glifos y trazo visible.
- Clasificar los entornos de nitrógeno relevantes usando el grafo químico y seleccionar direcciones locales planas para amidas/ureas cuando corresponda, conservando centros existentes y longitudes de referencia.
- Sustituir sólo el recorte excesivo confirmado por una intersección rayo–texto/forma ocupada más precisa y un margen óptico acotado.
- Verificar preview/commit, estilos actuales, persistencia y temas/escalas sin cambiar la silueta aprobada de beta.

## Fuera de alcance

No invocar Clean2D, modificar la identidad química, valencias, hidrógenos asignados, formato `.cmsn`, interfaz 3D, tipografía, colores, anchajes aprobados, estilo de enlaces múltiples, cuñas, hashes, junctions o presets globales. No mover átomos existentes ni acortar enlaces centro-a-centro para resolver colisiones de etiqueta. No aplicar una regla amida a aminas, nitrógenos aromáticos o centros coordinados sin evidencia y clasificación propia.

## Módulos probables

`src/chemuson/gui/canvas/canvas_bond_geometry.py`, `canvas/canvas_text.py`, `canvas/canvas_structure.py`, `src/chemuson/gui/items.py` y pruebas Qt/geométricas de canvas, depiction y persistencia. No se prevén paquetes nuevos ni dependencias.
