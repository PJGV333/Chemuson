# Auditoría de resiliencia operativa

Fecha: 2026-09-19

## Ownership explícito

- **M22 — resilience:** crash logging, excepthook, autosave rotativo,
  recuperación de snapshots y aislamiento de fallos de runtime.
- **M08/M10 — GUI y controllers:** containment de tareas Qt, workers,
  diálogos de recuperación y coordinación de lifecycle. No se mueve al núcleo
  de M22 porque requiere widgets, tabs o controllers.
- **M14 — update:** telemetry específica del subsistema de actualizaciones,
  separada de crash logging y observabilidad general.
- **M19 — bootstrap:** instalación inicial del crash hook y composición del
  event loop; no posee la implementación de resiliencia.

## Decisión

La frontera operacional es suficiente con M22 más los owners existentes de GUI,
controllers, bootstrap y update. No new module is created. M23 y M24 siguen
reservados para futuras fronteras justificadas.
