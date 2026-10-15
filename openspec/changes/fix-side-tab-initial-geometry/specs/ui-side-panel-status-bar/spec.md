# Delta de especificación — geometría inicial de SideTabRow

## MODIFIED Requirements

### Requirement: Los cinco tabs principales permanecen legibles y sin clipping

`SideTabRow` SHALL calcular y aplicar geometría compatible con el texto pintado por el estilo/fuente Qt vigentes después de aplicar tema y al revelar SidePanel. El primer paint después de abrir desde Ver SHALL mostrar inmediatamente los cinco labels completos si caben en el viewport, sin requerir clic/selección de pestaña. Si no caben, SHALL conservar el mecanismo existente de navegación/scroll y tooltips/nombres accesibles completos, sin pintar fuera del viewport. Cambios de tamaño, fuente o tema SHALL volver a sincronizar strip, botones y viewport de forma determinista.

#### Scenario: Primera presentación desde Ver antes de cualquier clic
- **GIVEN** SidePanel inicia oculto y la ventana ya aplicó el tema
- **WHEN** el usuario abre el panel desde Ver
- **THEN** antes de que seleccione una pestaña, cada label que quepa está completo y sus geometrías están dentro del viewport
- **AND** la actualización no cambia la geometría del canvas ni del rail.

#### Scenario: Recalcular tras cambios de presentación
- **GIVEN** SidePanel abierto o vuelto a abrir
- **WHEN** cambia la visibilidad, tamaño de ventana, fuente/stylesheet o tema light/dark
- **THEN** la geometría converge sin clics de pestaña
- **AND** overflow, tooltips y navegación continúan funcionales incluso cuando el viewport es más angosto que el strip.
