# M02 - Clean2D

## Responsabilidad
Clean2D es owner de depiction quality, ranking e imported depiction, incluido scaffold y block unwrap.

## Dependencias
- Runtime: `core` (M00), `chemio` (M01) para parsing y workers RDKit aislados.
- La dirección es exclusivamente M02 -> M01; no existe ciclo.

## Métricas
La deuda pasó de seis excepciones a una: cero runtime, una TYPE_CHECKING de persistencia y cero ciclos.
