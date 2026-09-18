## Context

`CanvasSelectionMixin` implementa siete reglas deterministas dispersas entre
la rotación de fragmentos y el escalado de selecciones:

- `_signed_angle_delta_deg`
- `_rotate_scene_point`
- `_optional_float_equal`
- `_point_equal`
- `_scale_point_from_anchor`
- `_normalize_label_scale`
- `_normalize_custom_stroke`

Las cinco primeras ya son estáticas. `_normalize_label_scale` no consulta
`self`. `_normalize_custom_stroke` sólo lee
`self.drawing_style.stroke_px`; la regla puede recibir ese valor como
argumento explícito. Otros mixins consumen algunos nombres a través del MRO de
`ChemusonCanvas`, por lo que sus puntos de llamada privados deben conservarse.

## Goals / Non-Goals

**Goals:**

- Dar ownership explícito a reglas numéricas sin estado.
- Preservar exactamente firmas privadas, tolerancias, límites y resultados.
- Mantener `QPointF` como tipo de valor para evitar conversiones.
- Hacer comprobable la lógica sin importar QtGui, QtWidgets ni el canvas.

**Non-Goals:**

- Mover handlers, drag, hit testing o overlays.
- Cambiar comandos undo/redo o transformación de objetos gráficos.
- Reorganizar rotación de fragmentos, ramas o trackball 3D.
- Cambiar escalado visual, química, renderizado o persistencia.
- Crear una API pública de M09.

## Decisions

### One internal value-geometry module

Las reglas vivirán en `gui/canvas/selection_geometry.py`. El módulo podrá
importar `math` y `PyQt6.QtCore.QPointF`, pero no QtGui, QtWidgets, el canvas,
modelos químicos, comandos ni items gráficos. Sus funciones serán
deterministas y no producirán efectos secundarios.

### Preserve private MRO compatibility

`CanvasSelectionMixin` conservará seis nombres como aliases estáticos:
delta angular, rotación, comparaciones, escalado y normalización de etiqueta.
Esto preserva llamadas desde `canvas_selection.py`, `canvas_text.py` y
`canvas_structure.py`.

`_normalize_custom_stroke` conservará su firma de un argumento como wrapper
de instancia porque debe pasar el grosor heredado actual a la función pura.

### Literal numeric extraction

No se cambiarán las fórmulas ni las tolerancias:

- delta angular en `[-180, 180)`;
- tolerancia opcional `0.05`;
- tolerancia de puntos `1e-4`;
- escala mínima de etiqueta `0.2` y herencia dentro de `0.02`;
- grosor mínimo `0.6` y herencia dentro de `0.05`.

## Risks / Trade-offs

- [Cambio numérico] → traslado literal más regresiones de límites, rotación,
  escalado y tolerancias.
- [Binding incorrecto] → contratos AST verifican aliases `staticmethod` y el
  wrapper de grosor.
- [Importación GUI accidental] → el módulo sólo admite QtCore como tipo de
  valor y prohíbe dependencias del canvas.
- [Ruptura entre mixins] → se preservan todos los nombres privados existentes.

## Migration Plan

1. Capturar baseline y consumidores.
2. Añadir pruebas numéricas y contratos AST inicialmente rojos.
3. Crear el módulo y trasladar literalmente las reglas.
4. Sustituir implementaciones por aliases/wrapper compatibles.
5. Actualizar M09 y ejecutar validaciones.

Rollback: devolver las siete reglas a `CanvasSelectionMixin` y retirar el
módulo interno; no existe migración de datos.
