# Diseño

1. Usar caracteres Unicode científicos en el `QTableWidgetItem` sólo para la fila/campo de fórmula, etiqueta NMR y aducto MS reconocidos. Conservar `_rows` intacto como fuente de copia TSV; los datos de magnitudes y desplazamientos no se transforman. Ajustar la columna Propiedad a contenido y dar el espacio restante a Valor para evitar clipping del texto científico según la fuente/plataforma. La tabla conserva selección, texto seleccionable y alineación estándar Qt.
2. Sustituir el `QHBoxLayout` canvas/panel por `QSplitter(Qt.Horizontal)` después del rail. Mantener el ancho inicial de 340 px, mínimo de 300 px y el canvas expandible; el splitter conserva su tamaño entre ocultar/mostrar dentro de la misma sesión. No persistir un nuevo ancho.
3. Cargar la pestaña activa guardada, pero iniciar `SidePanel` explícitamente oculto sin reescribir `ui/side_panel/visible`; las acciones de página de Ver siguen mostrando panel y página solicitada. No cambiar los contratos de `QSettings` existentes.
4. El tercer paso del onboarding conserva su posición y texto adaptado, pero dirige el foco al botón hamburguesa visible del AppBar para que el usuario abra `Ver`; no se muestra automáticamente el panel ni se usa su geometría oculta.
5. Asignar a cada botón primario de `SideTabRow` tooltip y nombre accesible idénticos a la etiqueta completa, sin alterar la navegación ni el control de overflow.

## Riesgos y límites

Qt puede devolver anchuras diferentes por handle/márgenes; las pruebas deben medir el splitter tras mostrarlo, en ventanas pequeñas/grandes, y verificar ocultar/mostrar sin pérdida de ancho. El formato Unicode es sólo de presentación; `copy_properties_tsv` sigue serializando `_rows` sin glifos de presentación. La aceptación visual final sigue siendo manual.
