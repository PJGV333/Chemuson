# Propuesta — disposición visual de H en etiquetas N–H

## Problema

En `Ar–NH–C(=O)–NH–R`, el elemento N conserva su anclaje óptico pero la etiqueta completa hereda un offset local y el H se dibuja siempre inline a la derecha (o a la izquierda por reflujo textual). En las capturas Qt de referencia ambos H quedan en una línea que compite visualmente con los enlaces vecinos; la corrección previa del recorte y sus 40 unidades centro-a-centro están presentes y deben conservarse.

## Alcance

- Determinar una posición visual para el H (lado, arriba o abajo de N) con base en enlaces vecinos y espacio gráfico real.
- Separar el glifo H del ancla N sin girar letras; centrar el glifo N sobre el centro del átomo.
- Cubrir `NH`, `HN`, `NH2` y `H2N` sin cambiar sus conteos ni la química de heteroátomos.
- Evitar enlaces, carbonilos, cargas y etiquetas vecinas usando las formas de dibujo existentes.
- Añadir reproducción QTest con gestos reales para urea N,N'-sustituida, acetanilida y aminas primaria/secundaria; capturas comparativas claro/oscuro y aserciones geométricas.

## Invariantes/fuera de alcance

No cambiar coordenadas de átomos, longitud o ángulo de enlaces, conectividad, valencias, Clean2D, estilos del renderer, persistencia CMSN ni undo/redo. No cambiar el recorte de etiquetas implementado por `improve-amide-urea-2d-continuity`; la forma ocupada del label sí debe incluir el H recolocado para conservar el despeje correcto. Otros heteroátomos con H implícitos mantienen su tratamiento actual.

## Módulos probables

`src/chemuson/gui/canvas/canvas_text.py`, `src/chemuson/gui/items.py`, pruebas de etiquetas/continuidad y evidencia en este cambio OpenSpec. Sin dependencias ni paquetes nuevos.
