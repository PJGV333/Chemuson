# CI-WORKER-LIFECYCLE — contratos de regresión del worker asíncrono

## Purpose

Hacer que las regresiones existentes de CompChem y cierre de workers midan el contrato real de señales/vida útil de Qt sin dobles obsoletos, timeouts arbitrarios ni supresión de fallos.

## ADDED Requirements

### Requirement: Worker doubles implement the production signal contract

Un doble de `CompChem3DWorker` que se inyecte en `CompChem3DController.start_job` SHALL exponer las señales conectadas por el controller (`frame_ready`, `progress_ready`, `finished`) con firmas compatibles. Si el doble modela progreso o frames, SHALL emitirlos sólo según el escenario declarado. La prueba MUST fallar antes de editar producción si el doble queda desincronizado del contrato.

#### Scenario: Descendant-worker close test substitutes CompChem worker
- **GIVEN** el controller conecta `progress_ready` aunque el worker no emita progreso en esta prueba
- **WHEN** el worker simulado se instala con monkeypatch
- **THEN** la conexión funciona y todos los workers se mantienen retenidos hasta que la prueba inicia cierre/liberación
- **AND** se verifica la terminación de los seis QThreads y la ausencia de callbacks/UI tardíos.

### Requirement: CompChem async completion is bounded and fully observed

La prueba async de CompChem SHALL conservar un límite finito y comprobar la ruta `worker.finished` → `QThread.finished` → `job_finished`, que no queden jobs activos ni resultados pendientes. Si no llega una señal, el fallo MUST conservar estado, señales recibidas y límite utilizado; MUST NOT convertirse en skip, éxito parcial o un timeout general ampliado.

#### Scenario: Fake backend returns a completed conformer
- **WHEN** el fake backend responde con coordenadas válidas
- **THEN** el worker concluye, el thread termina y el controller entrega un resultado válido dentro del límite existente; `active_jobs` y `pending_results` quedan vacíos.

#### Scenario: A Qt signal remains pending
- **WHEN** la señal worker llega pero la señal del thread/controller no llega al límite
- **THEN** la regresión conserva el estado real del QThread y del controller para diagnóstico, y el test falla hasta que la causa sea corregida o el entorno quede validado con evidencia repetible.

#### Scenario: Worker completes while the GUI event loop is not being pumped
- **WHEN** `worker.finished` se emite desde el QThread y el hilo principal está esperando o temporalmente ocupado
- **THEN** la solicitud thread-safe de salida del QThread se procesa desde el hilo emisor y no queda encolada detrás del event loop GUI
- **AND** el QThread puede terminar antes de que el hilo principal procese `job_finished`/`thread.finished`; esos callbacks siguen entregándose al reanudar el event loop.

### Requirement: Isolated shutdown wrapper preserves child failures

El test que lanza la suite de cierre en un proceso pytest independiente SHALL propagar exactamente el exit status y salida del hijo. El hijo SHALL ejecutar todos los casos de cierre previstos y limpiar workers, ventanas y recursos; un fallo del wrapper no se registrará como causa raíz independiente cuando sólo replique el mismo node ID fallido.

#### Scenario: Full async shutdown file passes in a fresh process
- **WHEN** el proceso hijo ejecuta los casos de cancelación y cierre con workers retenidos
- **THEN** cada worker se completa y el proceso devuelve cero; cualquier aserción del hijo permanece visible en el test padre.

### Requirement: CI outcomes remain complete and fail-closed

La colecta del manifiesto y los ocho shards SHALL asignar/ejecutar cada ID exactamente una vez. JUnit/log/timings/report SHALL conservar el estado de los fallos reales. No se permite ocultar el error con skips/xfail, cambiar node IDs a ciegas, aumentar todos los límites, ni aceptar el test workflow como verde si cualquier shard falla.

#### Scenario: Candidate CI contains a failing node
- **WHEN** un shard termina con exit no cero, fallo, timeout, señal o cobertura incompleta
- **THEN** el resumen reporta failure y las pruebas siguen formando parte de la matriz.

## Invariants

- No se cambia el comportamiento del worker/GUI en producción salvo evidencia reproducible aprobada en esta campaña.
- No se omite ningún test ni se modifica el manifiesto para cambiar el resultado de una prueba fallida.
- El timeout externo de shard existente permanece en 300 s.
- Cierre de ventana empaquetada, aceptación manual y test harness Qt se informan por separado.
