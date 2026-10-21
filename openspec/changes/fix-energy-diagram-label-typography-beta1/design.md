# Diseño — etiquetas de presentación separadas de presets

## Contrato actual

`diagram_presets.py` define IDs estables (`name`), nombres ASCII (`display_name`) y datos de construcción separados (`title`, parámetros científicos). Los submenús de `SymbolPaletteToolbar` usan actualmente `display_name` en `QAction`, pero capturan `preset.name` al emitir `electronic_diagram_preset_requested`. Los presets de herramienta (`energy_diagrams.py`) separan `kind`, `display_name` y `label`; la UI consume el nombre en tooltips y el `label` también alimenta la creación/persistencia de un diagrama.

## Decisiones

1. Mantener inalterados los objetos, nombres, IDs, títulos, parámetros y labels de todos los presets. No sustituir `display_name` ni etiquetas de diagramas persistentes.
2. Introducir funciones específicas de presentación de menú, basadas en IDs/kinds conocidos y no en una sustitución indiscriminada de dígitos. Aplicarlas a `QAction.setText()` de los presets semánticos y a tooltips/glifos de la paleta de cajas. Mantener `energy_diagram_display_name()` sin cambios porque también alimenta detalles de selección y estado fuera de este menú.
3. Conectar cada `QAction` con el mismo ID ASCII original. La representación Unicode (`₂`, `⁻`, `²`, `⁴`–`⁸`) es solo UI; HTML literal no se pasa a `QAction`.
4. No formatear `title`, `DiagramLevel.id/label`, texto editable del objeto, serialización, payloads, orden orbital ni valores calculados. `2s` y `2p` permanecen en línea base.
5. Mantener las fuentes, colores, estilos, tamaños y dimensiones existentes. Verificar métricas/layout de acciones y glifos mediante Qt; capturas offscreen documentan evidencia, sin presentarlas como aceptación nativa de KDE.

## Evidencia visual

Capturar los submenús MO y campo de ligandos y la paleta híbrida antes y después en tema claro/oscuro. Registrar fuente Qt efectiva, escala/DPR y verificar bounds sin recorte y glifos no vacíos; si el escritorio KDE real no está disponible, declarar esa limitación.

## Riesgos

- Modificar `display_name` o `title` directamente podría cambiar contratos ajenos a la UI; por eso la traducción visual vive en helpers de presentación.
- Un reemplazo general de números tipográficos dañaría `2s`/`2p`; solo se admiten mapeos explícitos para los IDs y tokens aprobados.
- El superíndice en la paleta de iconos puede ser más ancho; no se cambia tamaño ni fuente, y una regresión de recorte debe detectarse en pruebas/capturas.
