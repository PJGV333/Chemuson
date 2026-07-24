## Context

`ChemusonWindow` define `_coords_center`, `_average_bond_length`,
`_rescale_coords_to_bond_length`, `_align_coords_to_reference` y
`_project_missing_hydrogen_coords`. Aunque se exponen como helpers privados de
la clase, sólo operan sobre coordenadas, enlaces y símbolos de elemento.
Clean2DController y pruebas existentes acceden a ellas mediante la ventana.

## Goals / Non-Goals

**Goals:**

- Dar ownership explícito a las transformaciones geométricas puras.
- Preservar exactamente entradas, salidas, límites numéricos y reflexión.
- Mantener compatibilidad temporal con los helpers privados de la ventana.
- Evitar dependencias de Qt, `main_window` o subsistemas de dominio.

**Non-Goals:**

- Cambiar algoritmos Clean2D, thresholds o selección de candidatos.
- Unificar helpers homónimos de `clean2d` o `Clean2DController`.
- Cambiar llamadas de controllers, canvas o tests funcionales.
- Crear una API pública de M08.

## Decisions

### One pure internal module

Las transformaciones vivirán en `gui/clean2d_geometry.py`. El módulo sólo
importará `math` y tipos estándar; no importará Qt, la ventana, el canvas ni
modelos químicos.

### Preserve private window compatibility

`ChemusonWindow` conservará los cinco nombres privados como aliases estáticos
de las funciones extraídas. Esto permite que `Clean2DController` y las pruebas
existentes mantengan sus puntos de llamada durante esta fase.

### Literal extraction before deduplication

No se intentará consolidar `_coords_center` ni otros helpers similares de
`Clean2DController` o del motor Clean2D. Esa deduplicación requeriría comparar
contratos semánticos y queda fuera del alcance estructural.

## Risks / Trade-offs

- [Cambio numérico] → las implementaciones se trasladan literalmente y se
  ejecutan las regresiones de escala, alineación e hidrógenos.
- [API accidental] → el módulo se registra como `internal_api`; la API pública
  de M08 sigue siendo únicamente `ChemusonWindow`.
- [Ciclo] → el módulo puro no importa `main_window`.
- [Alias con binding incorrecto] → pruebas AST y de ejecución verifican que los
  aliases sean estáticos y acepten las mismas llamadas.

## Migration Plan

1. Capturar baseline y consumidores.
2. Añadir caracterización arquitectónica.
3. Trasladar literalmente las cinco funciones.
4. Sustituir los métodos por aliases estáticos privados.
5. Actualizar M08 y ejecutar validaciones.

Rollback: devolver las funciones al cuerpo de `ChemusonWindow` y retirar el
módulo interno; no existe migración de datos.
