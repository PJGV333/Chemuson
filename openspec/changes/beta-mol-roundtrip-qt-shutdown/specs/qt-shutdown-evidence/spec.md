# Spec Delta — evidencia de cierre y recuperación Qt

## ADDED Requirements

### Requirement: Cierre normal tras flujo de archivo
La ventana y sus recursos Qt SHALL poder cerrarse tras ventana limpia, exportación MOL, documento modificado con confirmación aprobada y apertura de snapshot autosave recuperado. Las pruebas con posibilidad de bloqueo SHALL ejecutarse en subprocess con timeout externo.

#### Scenario: Cierre en varios estados de documento
- **WHEN** se cierra una ventana limpia, exportada, modificada con descarte aprobado o recuperada
- **THEN** el proceso termina limpiamente dentro del timeout

### Requirement: Destrucción segura de canvas y señales async
La destrucción de canvas/ventana SHALL esperar los QThreads descendientes según el coordinador existente y no SHALL permitir callbacks tardíos sobre widgets destruidos. Las pruebas cubrirán desconexión/retención de señales y salida de proceso, no sólo `window.close()`.

#### Scenario: Canvas y resultados tardíos
- **WHEN** el canvas se destruye o un worker controlado termina durante el cierre
- **THEN** no se invoca ningún callback sobre widgets destruidos y ningún QThread sobrevive a su propietario

### Requirement: Investigación de SIGSEGV no reproducido
Un SIGSEGV manual reportado SHALL permanecer riesgo abierto si no se reproduce ni se obtiene causa respaldada. La ausencia de reproducción no basta para declararlo corregido. Una traza concreta SHALL guiar cualquier cambio productivo.

#### Scenario: Crash reproducible
- **WHEN** una reproducción controlada termina con SIGSEGV
- **THEN** el informe identifica código de salida, etapa y traza antes de aceptar una corrección

### Requirement: Recuperación sin duplicar estructura
Reabrir un snapshot autosave SHALL preservar el grafo, sin duplicar nodos/enlaces ni crear etiquetas de carbono redundantes. Si la representación visual observada no es reproducible, el informe SHALL dejarlo explícito.

#### Scenario: Recuperar autosave
- **WHEN** se restaura un snapshot creado por autosave
- **THEN** el grafo y el número/visibilidad de etiquetas corresponden al snapshot sin duplicación

## Invariantes

- No se modifica formato o sistema de persistencia por una observación no reproducida.
- Sin `QThread.terminate()`, supresión de crash, catch indiscriminado ni sleeps arbitrarios.
- Una salida no cero o timeout nunca se cuenta como PASS.
