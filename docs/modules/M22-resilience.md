# M22 — Resiliencia y recuperación

## Responsabilidad

M22 (`resilience`) posee autosave rotativo, recuperación de documentos y
registro de crashes. Sus servicios canónicos no importan módulos GUI; el
notificador de crash sólo usa Qt como dependencia externa.

## API

- `AutosaveManager` y sus protocolos de colaboración.
- `read_autosave_metadata`, `list_autosave_entries` y `archive_autosave` para
  leer, ordenar y archivar snapshots sin depender de la GUI.
- `install` y `write_crash_log` para el excepthook y reportes.

## Compatibilidad

`chemuson.utils.autosave` y `chemuson.utils.crash_reporter` permanecen como
shims de importación. No contienen implementaciones duplicadas.

## Límites

M22 no posee widgets, tabs, PersistenceManager ni composición de aplicación.
`recovery.py` sólo contiene política filesystem de snapshots y no importa
`chemuson.gui` ni `PersistenceManager`. La GUI y el bootstrap pueden consumir
M22; M22 no importa `chemuson.gui`.
