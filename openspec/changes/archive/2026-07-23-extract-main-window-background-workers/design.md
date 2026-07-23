## Context

`main_window.py` posee `_DescriptorWorker` y `_NameToStructureWorker`. Ambos
son `QObject` privados con una señal `finished`, estado de entrada mínimo y un
método `run()` que difiere el import de su backend. La ventana conserva los
`QThread`, conexiones, jobs y presentación de resultados.

## Goals / Non-Goals

**Goals:**

- Dar ownership explícito a las implementaciones de los workers.
- Preservar exactamente constructor, señal, payload, timeout y errores.
- Mantener el lifecycle de threads y toda respuesta UI en la ventana.
- Evitar un ciclo `background_workers` → `main_window`.

**Non-Goals:**

- Rediseñar la infraestructura asíncrona.
- Introducir pools, tareas genéricas, DI o dependencias nuevas.
- Cambiar cancelación, red, RDKit o Name→Structure.
- Mover handlers, docks, canvas o comportamiento químico.

## Decisions

### One internal worker module

Los dos workers vivirán juntos en `gui/background_workers.py`: ambos son
adaptadores Qt privados, pequeños y poseídos por M08. No se crea un paquete ni
una API pública nueva.

### Preserve private identities

Los nombres `_DescriptorWorker` y `_NameToStructureWorker` permanecen privados.
`main_window.py` los importa explícitamente para conservar los puntos de
construcción existentes.

### Keep orchestration in the window

`ChemusonWindow` seguirá creando `QThread`, conectando señales, registrando jobs
y mostrando progreso o resultados. El módulo nuevo sólo ejecuta el backend y
emite el resultado.

### Preserve deferred backend imports

Los imports de `molecular_descriptors_isolated` y
`resolve_name_to_structure` seguirán dentro de `run()`. Importar la ventana no
debe adelantar RDKit ni el resolvedor.

## Risks / Trade-offs

- [Cambio de señal] → prueba AST verifica argumentos y emisión.
- [Ciclo GUI] → el módulo nuevo no importa `main_window`.
- [Cambio de threading] → los métodos `_start_*_job` no se trasladan.
- [API accidental] → el catálogo y `__all__` conservan sólo
  `ChemusonWindow` como API pública de esta costura.

## Migration Plan

1. Capturar baseline y reconocer consumidores.
2. Añadir caracterización arquitectónica.
3. Trasladar literalmente ambos workers.
4. Actualizar M08 y documentación.
5. Ejecutar validaciones focalizadas y globales disponibles.

Rollback: devolver las dos clases a `main_window.py` y retirar el módulo
interno; no existe migración de datos.
