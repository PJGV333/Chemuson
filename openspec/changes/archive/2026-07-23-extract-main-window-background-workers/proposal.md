## Why

Después de extraer el application shell, `gui/main_window.py` todavía define
dos objetos `QObject` dedicados a ejecutar trabajo bloqueante en segundo
plano. Esto mezcla la coordinación de la ventana con las implementaciones de
descriptores RDKit y Name→Structure, y mantiene imports de señales Qt que no
pertenecen a `ChemusonWindow`.

## What Changes

- Crear `gui/background_workers.py` como implementación interna de M08.
- Trasladar literalmente `_DescriptorWorker` y `_NameToStructureWorker`.
- Mantener en `ChemusonWindow` la creación de `QThread`, los diálogos de
  progreso, la cancelación y los callbacks que actualizan la interfaz.
- Caracterizar ownership, señales, imports diferidos y ausencia de ciclos.
- Actualizar el catálogo y la ficha M08.

## Capabilities

### New Capabilities

- `gui-background-workers`: workers Qt internos y acotados para operaciones
  bloqueantes iniciadas por la ventana.

### Modified Capabilities

- `module-catalog`: M08 registra el módulo y su prueba arquitectónica.

## Impact

El cambio reduce `main_window.py` sin alterar señales, timeouts, threading,
red, RDKit, Name→Structure, UI, canvas, química, persistencia ni formatos.
