## Why

Después de extraer el application shell y los workers de segundo plano,
`gui/main_window.py` todavía implementa cuatro transformaciones geométricas
puras usadas al integrar coordenadas Clean2D y proyecciones 3D. Estas
funciones no dependen de una ventana ni de Qt, pero su ubicación hace que
`ChemusonWindow` posea lógica numérica además de coordinación de interfaz.

## What Changes

- Crear `gui/clean2d_geometry.py` como implementación interna y pura de M08.
- Trasladar literalmente centro, longitud media, reescalado, alineación rígida
  y proyección de hidrógenos faltantes.
- Conservar aliases estáticos privados en `ChemusonWindow` para compatibilidad
  con consumidores y pruebas existentes.
- Caracterizar pureza, ownership, firmas y ausencia de ciclos.
- Actualizar el catálogo y la ficha M08.

## Capabilities

### New Capabilities

- `gui-clean2d-geometry`: transformaciones geométricas puras usadas por la
  integración GUI de Clean2D.

### Modified Capabilities

- `module-catalog`: M08 registra el módulo interno y su caracterización.

## Impact

El cambio reduce `main_window.py` sin alterar coordenadas resultantes, canvas,
renderizado, química, persistencia, UI, acciones ni formatos.
