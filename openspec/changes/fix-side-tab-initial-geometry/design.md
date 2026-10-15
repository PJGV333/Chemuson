# Diseño

1. Inspeccionar el orden real: `SideTabRow` fija anchos con `QFontMetrics` en el constructor, ajusta `_strip` y su mínimo una sola vez; el tema/QSS principal se aplica después durante el ensamblaje y el panel permanece oculto hasta la acción Ver. Medir el estado inicial antes de clicks y tras `StyleChange`, `FontChange`, show y resize.
2. Implementar una única rutina idempotente de sincronización de geometría basada en la fuente tipográfica efectiva y el estilo Qt aplicado. Actualizar anchos, layout del strip y su tamaño mínimo juntos; no depender de `set_active()` como side effect para corregir el primer paint.
3. Invocar la sincronización desde los eventos necesarios y diferirla al siguiente ciclo Qt sólo si las mediciones muestran que QSS/layout todavía no está asentado. Proteger contra ciclos/reentrancia y conservar los botones existentes, su orden, el QScrollArea, overflow, tooltips y navegación.
4. La prueba de regresión abre desde Ver, procesa el primer paint y mide etiquetas/anchos sin seleccionar ninguna pestaña; repite hide/show, resize y light/dark. También mide el rect del canvas antes/después para asegurar que el refresco no modifica el splitter.

## Riesgo

QFontMetrics calcula advance del texto pero no incluye necesariamente padding/bordes del QSS. `adjustSize()` antes de aplicar el stylesheet puede dejar `_strip` y su `minimumWidth` desactualizados. La prueba debe identificar si el fallo es el ancho del botón, el viewport/strip, o ambos antes de elegir qué geometría refrescar.
