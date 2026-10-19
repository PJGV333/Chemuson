# Validación — geometría inicial 2D

## Registro antes/después

| Caso | Ángulo previsto | Antes | Después | Longitud prevista/obtenida | Cruces/excepción |
|---|---:|---:|---:|---|---|
| Benceno aislado | 120° | 120° (máx. error 1.2e-13°) | Pendiente | 40.0 / 40.0 | Cero; sin excepción |
| Benceno con sustituyente, vértice rotado | 120° local y radial | 137°/103°, radial error 17° | Pendiente | 40.0 / pendiente | Pendiente |
| sp3/zigzag acíclico | 120° | 109.5° | Pendiente | 40.0 / pendiente | Pendiente |
| Centro carbonilo | ~120° | 119.69–120.18° en gestos de rejilla | Pendiente | 40.0 / ~40 | Pendiente |
| 1-fenilbutan-2-ona completa | conectividad C10H12O y geometría 2D coherente | 11 átomos/11 enlaces; gesto dirigido 119.69–120.32° | Pendiente | 40.0 / ~40 | Pendiente |

## Baseline

Ver `baseline.md`, logs completos en `baseline/` y capturas/modelo en `evidence/before/`. El fallo de menú beta, timeout monolítico y F401 Ruff son basales y se mantienen separados de cualquier regresión de esta campaña.

## Evidencia posterior y comandos

Pendiente de implementación. Se registrarán comandos y exit codes exactos, nodos pytest, matriz de escala/tema y hashes de capturas. Ningún passing de QImage/Qt offscreen se presentará como aceptación visual manual del propietario. La suite completa no se ejecutará sin límite; todas las pruebas relevantes usarán timeout externo.
