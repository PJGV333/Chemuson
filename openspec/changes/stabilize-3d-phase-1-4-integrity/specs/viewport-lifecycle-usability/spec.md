## ADDED Requirements

### Requirement: La investigación de cierres nativos distingue causas y conserva evidencia

Los SIGSEGV de la suite de pruebas y del cierre Qt/SIP SHALL tratarse como incidentes independientes hasta disponer de evidencia diferencial. El diagnóstico SHALL incluir stack/core disponibles, estado de workers y renderer, alcance de la reproducción y limitaciones. SHALL NOT ocultar crashes, forzar salida del proceso, instalar componentes, elevar privilegios ni eludir una aprobación de GUI. Pruebas offscreen SHALL NOT ser criterio suficiente para declarar estable el cierre en display real.

#### Scenario: Reproducción gráfica bloqueada
- **GIVEN** que la herramienta solicita una aprobación adicional para abrir GUI y ésta no se concede
- **WHEN** continúa la campaña
- **THEN** no se intenta otra vía gráfica
- **AND** se documenta la guía manual y se deja la aceptación visual abierta

#### Scenario: Suite completa conocida
- **GIVEN** un baseline previo con SIGSEGV en `concurrent.futures.thread._worker`
- **WHEN** se ejecutan pruebas de esta estabilización
- **THEN** se priorizan shards focalizados y se informa la suite monolítica como pendiente/conocida sin reruns innecesarios

### Requirement: El dock 3D prioriza el viewport y mantiene controles principales

El dock SHALL reservar la mayor parte del espacio para el viewport, agrupar controles secundarios en secciones compactas/contextuales y mantener edición y Live MM accesibles. La UI SHALL evitar que el layout o el resize oculten botones principales o degraden el canvas 2D. Contraste, color por elemento, selección, clipping y enlaces múltiples SHALL seguir siendo distinguibles.

#### Scenario: Dock estrecho
- **GIVEN** un dock 3D estrecho y redimensionado junto al canvas 2D
- **WHEN** se abre la vista
- **THEN** el viewport conserva el área prioritaria y los controles no quedan cortados
- **AND** Live MM permanece visible y operable

#### Scenario: Render de contraste
- **GIVEN** una escena con carbonos, heteroátomos, enlaces múltiples y selección
- **WHEN** se presenta sobre el tema oscuro actual
- **THEN** átomos, enlaces y resaltado pueden distinguirse
- **AND** la aceptación de iluminación/layout se confirma sólo en display real