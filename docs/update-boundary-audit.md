# Auditoría de la frontera de auto-actualización

Fecha: 2026-09-19

## Decisión

M14 (`update`) está **audited / no structural change required**. El paquete ya
separa de forma coherente:

- `types.py`: contratos y enums.
- `semver.py` y `policy.py`: decisión pura de elegibilidad.
- `provider.py`: consulta y cache de releases.
- `security.py`: HTTPS, hashes y firmas.
- `portable.py`, `windows.py` y `rollback.py`: aplicación por plataforma.
- `telemetry.py`: registro de eventos de actualización.
- `core.py`: orquestación del flujo.

## Ownership y dependencias

La ruta canónica es `src/chemuson/update/`, el catálogo conserva M14 y la
superficie pública se mantiene en `chemuson.update`. No hay dependencias GUI en
el paquete ni una excepción arquitectónica pendiente.

## M23

M23 queda reservado para una futura frontera que demuestre una cohesión nueva.
No se crea un módulo vacío ni se fragmenta M14 sólo para aumentar el catálogo.
No new module is created.
