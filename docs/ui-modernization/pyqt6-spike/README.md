# Referencia visual archivada — PyQt6 UI spike

El spike se usó en septiembre de 2026 para decidir si PyQt6/QtWidgets podía expresar el nuevo lenguaje visual antes de modificar producción. La modernización Fases 0–8 ya está cerrada; el demo ejecutable y sus widgets/iconos de prueba se retiraron porque no son código de producto ni tienen consumidores. No se documenta como herramienta soportada ni se puede ejecutar desde este directorio.

`theme.py` se conserva como fuente exacta de tokens/métricas de la propuesta: la especificación global `openspec/specs/ui-theme-foundation/spec.md` lo sigue identificando como referencia normativa. No importar este archivo desde la aplicación.

## Decisiones validadas

- Mantener PyQt6 y `QGraphicsView` para el editor; evitar migración de toolkit y doble modelo QML/web.
- Tokens claro/oscuro centralizados, SVG y métricas coherentes; prototipo sujeto a comparación real, no a pixel-perfect.
- Shell final: app bar de 54 px, rail de 58 px, botones de 42 px, estado de 34 px; la app real conserva canvas, acciones y handlers en lugar de simularlos.
- `CommandPalette` terminó usando `Ctrl+P`, no el `Ctrl+K` inicial, que sigue asignado a Clean2D.
- Las sombras, fuentes, transparencia y viewport/QSS requieren verificación en Qt real; el modo offscreen no sustituye una revisión KDE/Wayland.

La narración compacta, los fallos experimentales y el resultado final están en [`../../history/CAMPAIGNS.md`](../../history/CAMPAIGNS.md). Evidencia final: [`../after/README.md`](../after/README.md). El spike no debe ampliarse ni regenerarse como una fase implícita.
