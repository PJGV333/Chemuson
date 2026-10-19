## ADDED Requirements

### Requirement: El dock 3D permanece legible acoplado y flotante

El viewport SHALL aplicar estilos claros/oscuros scoped al `MolecularViewportDock`, sus controles y estados disabled/hover/checked, y SHALL reaplicarlos al cambiar tema o `topLevelChanged`. No SHALL modificar paleta global ni estilos de otros docks. Flotar, mover, redimensionar o acoplar SHALL preservar herramienta, escena, cámara y selección; la destrucción/recreación del contexto OpenGL SHALL liberar recursos GPU y permitir reinitializar renderer y subir de nuevo la escena retenida.

#### Scenario: El tema se mantiene al desacoplar
- **GIVEN** el dock acoplado con una herramienta seleccionada y una escena
- **WHEN** se desacopla, se redimensiona y se vuelve a acoplar en tema claro u oscuro
- **THEN** textos, checks, botones, combos y estados disabled conservan contraste
- **AND** la herramienta y escena siguen seleccionadas y visibles tras reconstruir el contexto
- **AND** los otros docks conservan sus estilos

#### Scenario: Cambio de tema con dock flotante
- **GIVEN** el viewport como ventana flotante
- **WHEN** el usuario cambia entre tema claro y oscuro
- **THEN** el QSS scoped/paleta del dock se actualizan sin modificar configuración del sistema
