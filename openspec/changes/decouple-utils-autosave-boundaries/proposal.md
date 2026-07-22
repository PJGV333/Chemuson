# Decouple Utils Autosave Boundaries

## Why

M15 (`utils`) declara una arquitectura objetivo sin dependencias ChemUSON, pero
`utils.autosave` importa `chemio.persistence` en runtime y `gui.canvas` bajo
`TYPE_CHECKING`. Ambas dependencias son deuda temporal documentada, impiden que
autosave sea reutilizable e impiden reducir el baseline de excepciones.

## What Changes

- Invertir las dependencias de autosave mediante contratos estructurales y
  colaboradores inyectados desde el punto de composición GUI.
- Conservar comportamiento, rutas, nombres, frecuencia, recuperación y formato
  JSON/CMSN de autosave.
- Eliminar las dos excepciones M15, dejar sus dependencias actuales y objetivo
  vacías, y reducir el baseline congelado de 10 a 8 identidades.
- Añadir pruebas de comportamiento e aislamiento de importación.

## Non-Goals

- Rediseñar la persistencia completa o mover `autosave.py`.
- Reescribir canvas, cambiar formato CMSN, UI, preferencias, recuperación,
  undo/redo o frecuencia de autosave.
- Eliminar las otras ocho excepciones ni resolver los ciclos M01<->M02 o
  M03<->M04.
