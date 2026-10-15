# Propuesta — geometría inicial de pestañas del SidePanel

## Problema

En la rama `fix/beta-properties-panel-polish`, al abrir SidePanel desde Ver las cinco pestañas primarias pueden aparecer truncadas hasta que un clic en una pestaña provoca una actualización de estilo/geometría. No debe requerirse interacción sobre las pestañas para que se lean.

## Alcance

- Determinar por qué el ancho medido y el tamaño del strip no coinciden con el texto pintado después de aplicar el tema y revelar el panel.
- Actualizar la geometría/presentación determinísticamente al revelar y tras cambios de estilo, fuente o tamaño que afecten el viewport.
- Con espacio suficiente, mostrar completos los cinco labels antes de cambiar la página; con espacio insuficiente, mantener la navegación, el desplazamiento actual y tooltips completos.
- Añadir regresión sobre apertura desde Ver antes de cualquier clic de pestaña, incluyendo hide/show, resize y light/dark.

## Fuera de alcance

No cambiar splitter, preferencias, inicio oculto, lógica de navegación/overflow, química, cálculos, serialización/exportación ni otros subsistemas.
